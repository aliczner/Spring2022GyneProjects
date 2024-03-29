#wild spring queens 2022

#first need to separate out the right bees

allData <- read.csv("TriangulateHandheld.csv") #12245 rows
data <- subset(allData, projectID=="wild") #4727 rows

library(dplyr)
data2<-data %>% #checking for duplications, no duplicates were found
  distinct() 

#need to remove any points that are outside of the towers detection range
towers<-read.csv("TowerActualLocations.csv")

landcover<-terra::rast("SpringLandcoverRaster.tif") #Rare landcover, made in overwintering

library(sf)
towers <- st_as_sf(towers, coords = c("X", "Y"), crs = 4326)
towersproj <- st_transform(towers, crs = st_crs(landcover))
towersdf<-as.data.frame(towersproj) 

#making a 500 m buffer around towers
perim <-st_buffer(towersproj, dist=500) %>% st_union()

#making data spatial points
data.sp <- data
data.sp <- st_as_sf(data.sp, coords = c("X", "Y"), crs = 4326)
data.spj <- st_transform(data.sp, crs = st_crs(landcover))

#removing points not in the polygon, using datafrom below
# Perform the spatial intersection
result_sf <- st_intersection(data.spj, perim)

pointsIN.df<-data.frame(result_sf)
detsum<-aggregate(pointsIN.df$AnimalID, by=list(pointsIN.df$AnimalID), FUN=length)
detsum #number of observations per bee. 

write.csv(pointsIN.df, "WildSpring22Perim.csv")
tagsproj<-result_sf

library(dplyr)

nums<- pointsIN.df %>% 
  filter(!is.na(Var_X)) %>% 
  mutate(Var_X = as.numeric(Var_X)) %>%
  filter(!is.infinite(Var_X)) %>% 
  filter(!is.na(Var_Y)) %>% 
  mutate(Var_Y=as.numeric(Var_Y)) %>% 
  filter(!is.infinite(Var_Y))

nums %>% 
  group_by(AnimalID) %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

nums %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

############################################################
###               Bayes move                            ###
#############################################################

library(bayesmove)

######  Pre-specifying Breakpoints and Dealing with Zero-inflated Variables ####

library(lubridate)
library(tidyr)

#date is in the wrong format
#going to have to change date separately for towers vs. handheld gps. 

tags.df <- pointsIN.df
date_formats <- c("%Y-%m-%d", "%d_%m_%y") #the two date formats in the dataframe
parsed_dates <- parse_date_time(tags.df$GPSdate, orders = date_formats) 
output_format <- "%m-%d-%y"
formatted_dates <-format(parsed_dates, format = output_format) #change to format I want
tags.df$GPSdate <-formatted_dates #update GPSdate

#add a column that combines date and time together

tags.df$dtime <- mdy_hms(paste0(tags.df$GPSdate, tags.df$GPStime))

tagsgdate<-tags.df %>% 
  dplyr::select(longitude, latitude, dtime, AnimalID)
names(tagsgdate)[1] <- "x"
names(tagsgdate)[2] <- "y"
names(tagsgdate)[3] <- "date"
names(tagsgdate)[4] <- "id"
tagsgdate<-tagsgdate %>% arrange(id, date)

tagsNoNA <- tagsgdate %>% filter(!is.na(date))

tags.sf <- st_as_sf(tagsNoNA, coords = c("x", "y"), crs = 4326)
tags.pj <- st_transform(tags.sf, st_crs(landcover))
#here I just realized the geometry column at the end of tagsNoNA was the projected
#coordinates, so should do this earlier
tags.latlong <- st_coordinates(tags.pj)

tags.pj2 <- cbind(tags.pj, tags.latlong)
tags.pj2df <-as.data.frame(tags.pj2)
tags.dat<-tags.pj2df %>% 
  dplyr::select(date, id, X, Y)

### back to the bayes move 

#calculating step length, turning angle and time interval
tracks <-prep_data(dat=tags.dat %>% 
                     filter(id != "F062002"), coord.names=c("X","Y"), id="id") 
#can't run prep_data on any ids with just one occurrence.

head(tracks)
unique (tracks$id) #22 unique tracks

library(ggplot2)

ggplot() +
  geom_path(data = tracks, aes(x, y, color = id, group = id), size=0.75) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_d("id", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

boxplot(log10(tracks$dt))
tracks.big <- tracks %>% filter(dt < 10000 & dt > 1)
tracks.big2 <- tracks.big %>% filter(dt < 2000 & dt > 1)

#round time steps to specified interval
hist(tracks.big2$dt, breaks=100)

tracks <-round_track_time(dat=tracks, id="id", int=20, tol=3600, time.zone="UTC", units="secs")

#will break up the data into  rest and non-rest

tracks <- tracks %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
#now rest column where 1 = rest and 2 = non-rest

# Create list from data frame where each element is a different track
tracks.list <- df_to_list(dat = tracks, ind = "id")

# Filter observations to time interval
tracks_filt.list <- filter_time(dat.list = tracks.list, int = 20)

#### Discretize data streams (put into bins)

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
angle.bin.lims

dist.bin.lims=quantile(tracks[tracks$dt==20 & tracks$step != 0,]$step,
                       c(0.3334, 0.6667, 1.0), na.rm=T) 
dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list<- map(tracks_filt.list,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))
# Since 0s get lumped into bin 1 for SL, need to add an extra bin to store 0s
tracks_disc.list2 <- tracks_disc.list %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))

##pre-specify breakpoints

# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.list2<- map(tracks_disc.list2,
                   subset,
                   select = c(id, SL, TA, rest))

## Drop bees without non-rest state (i.e., no "2"s)
beesToDrop <- sapply(tracks.list2, function(x){
  ifelse(sum(x$rest) == nrow(x) *2 || sum(x$rest) == nrow(x), FALSE, TRUE)
})
tracks.list3 <- tracks.list2[beesToDrop]

tracks_disc.list3 <-tracks_disc.list2[beesToDrop]

# Pre-specify breakpoints based on 'rest'
breaks<- purrr::map(tracks.list3, ~find_breaks(dat = ., ind = "rest"))

### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(3,8,2)  #SL, TA, rest (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res<- segment_behavior(data = tracks.list3, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha, breakpt = breaks)

future::plan(future::sequential)  #return to single core

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res, type = "LML")

## Determine MAP for selecting breakpoints (maximum a posteriori)
MAP.est<- get_MAP(dat = dat.res$LML, nburn = 5000)
MAP.est

brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))

plot_breakpoints(data = tracks_disc.list2, as_date = FALSE, var_names = c("step","angle","rest"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Resting"), brkpts = brkpts)

tracks.seg<- assign_tseg(dat = tracks_disc.list3, brkpts = brkpts)
head(tracks.seg)

# Select only id, tseg, SL, TA, and rest columns
tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA","rest")]

# Summarize observations by track segment
nbins<- c(3,8,2)
obs<- summarize_tsegs(dat = tracks.seg2, nbins = nbins)
head(obs)

#### Running the clustering model

set.seed(1)

# Prepare for Gibbs sampler
ngibbs <- 20000  #number of MCMC iterations for Gibbs sampler
nburn <- ngibbs/2  #number of iterations for burn-in
nmaxclust <- max(nbins) - 1  #one fewer than max number of bins used for data streams
ndata.types <- length(nbins)  #number of data types

# Priors
gamma1 <- 0.1
alpha <- 0.1

# Run LDA model
res <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                        ngibbs=ngibbs, nmaxclust=nmaxclust,
                        nburn=nburn, ndata.types=ndata.types)

plot(res$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")

### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim <- extract_prop(res = res, ngibbs = ngibbs, nburn = nburn, nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means<- round(colMeans(theta.estim), digits = 3))

# Calculate cumulative sum
cumsum(theta.means)
#first 2 states comprise 94.6% of all observations

# Convert to data frame for ggplot2
theta.estim_df<- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", values_to = "prop") %>%
  modify_at("behaviour", factor)
levels(theta.estim_df$behaviour)<- 1:nmaxclust

# Plot results
ggplot(theta.estim_df, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in the first two behaviour states, some in 3

#### Classify the states as behaviours

# Extract bin estimates from phi matrix
behav.res<- get_behav_hist(dat = res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle","Resting"))

# Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c("#21908CFF","#440154FF","#FDE725FF",
                               "grey35","grey35","grey35","grey35",
                               "grey35","grey35","grey35","grey35")
                    , guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:4) +
  facet_grid(behav ~ var, scales = "free_x")

# behaviour 1 = ARS: not resting med to large step lengths, many turn angles
# behaviour 2 = resting

theta.estim.long <- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
                                   nbehav = 2, behav.names = c("ARS", "Rest"),
                                   behav.order = c(1,2))

# Plot results
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behaviour\n") +
  scale_fill_viridis_d("Behaviour") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, nrow = 2)


##### Assign behavioural states to tracks

# Merge results with original data
tracks.out<- assign_behavior(dat.orig = tracks,
                             dat.seg.list = df_to_list(tracks.seg, "id"),
                             theta.estim.long = theta.estim.long,
                             behav.names = c("ARS", "Rest"))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = tracks.out, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = tracks.out, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free", ncol = 2)

# Map of all IDs
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = behav, group = id), size=1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_d("Behaviour", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Proportion ARS
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = ARS, group = id), size=1.25, alpha=0.7) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == 1)),
             aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == n())),
             aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nARS", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "ARS") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Proportion resting
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = Rest, group = id), size=1.25, alpha=0.7) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == 1)),
             aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == n())),
             aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nRest", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "Resting") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

## conditions edits
## Bees to drop is really be bees to keep. So we need the inverse
beesToRevise <- names(beesToDrop)[! beesToDrop]

head(tracks.list2)
tracksAll <- do.call(rbind, tracks_disc.list2)

## Revise Rest == 1 for all observations
tracks.outTempRest <- tracksAll %>% 
                  mutate(angle = ifelse(is.na(angle), 0, angle))

tracks.outTempARS <- tracksAll %>% 
  mutate(angle = ifelse(is.na(angle), 0, angle))
tracks.outTempRest[tracks.outTempRest$id %in% beesToRevise &   
                     tracks.outTempRest$rest == 1 & 
                     tracks.outTempRest$dt == 20 &
                     tracks.outTempRest$angle > -1 & tracks.outTempRest$angle < 1, "behav" ] <- "Rest"
tracks.outTempRest <- tracks.outTempRest %>% 
                filter(!is.na(behav)) %>% 
                dplyr::select(-time1, -SL, -TA)

tracks.outTempARS[!tracks.outTempARS$id %in% tracks.outTempRest$id, "behav" ] <- "ARS"
tracks.outTempARS <- tracks.outTempARS %>% 
  filter(!is.na(behav)) %>% 
  dplyr::select(-time1, -SL, -TA)

## join processed dataset
tracks.outWithMissing <- tracks.out %>% 
  dplyr::select(names(tracks.outTempRest)) %>% 
  rbind(tracks.outTempRest)  %>% 
  rbind(tracks.outTempARS)
tracks.out2<-tracks.outWithMissing

write.csv(tracks.out2, "tracksout2wild.csv")

#################################################
#####         Flight Paths         ##################
#####################################################

library(ggplot2)
library(leaflet)
library(shiny)
library(dplyr)
library(htmltools)

## set leaflet CRS
UTMtoLatLon <- function(dfToConvert){
  dfIn <- dfToConvert
  coordinates(dfIn) <- ~x+y
  proj4string(dfIn) <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs "
  dfIntLatLon <- spTransform(dfIn, CRS("+proj=longlat +datum=WGS84"))
  dfLatLon <- data.frame(dfIntLatLon)
  return(dfLatLon)
}

## CSS code to style title in leaflet map
tag.map.title <- htmltools::tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(10%,20%);
    position: fixed !important;
    left: 70%;
    bottom: 10%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 16px;
  }
"))

i="CTL-02"
j=2
allIds <- unique(tracks.out2$id)
for( i in allIds) {
  beeSingular <- tracks.out2 %>% filter(id == i) %>% UTMtoLatLon(.) 
  ## need to find out number of steps before second loop
  ## maximum 65 frames for very long steps
  nSteps <- ifelse(nrow(beeSingular) > 65, 65, nrow(beeSingular))
  for( j in 2:nSteps) { ## need to start at 2 to have a line
    beeStep<-beeSingular[1:j,]
    
    colorSelect <- colorRampPalette(c("blue","orange"))
    beeStep$colour <- colorSelect(j)
    
    ## HTML code to add title - uses CSS styles from above
    title <- htmltools::tags$div(
      tag.map.title, HTML(paste0(beeStep[j, "date"]))
    )  
    
    leafletMap <- leaflet() %>% addTiles() %>% 
      addProviderTiles('Esri.WorldImagery') %>% 
      setView(-80.352, 43.377, zoom = 16) %>% 
      addControl(title, position = "topleft", className="map-title")
    for(k in 1:nrow(beeStep)){
      lineStart <- k
      lineEnd <- k+1
      leafletMap <- addPolylines(leafletMap, 
                                 data = beeStep[lineStart:lineEnd,], 
                                 lng = ~x, lat = ~y, color = ~colour)
    }
    leafletMap
    
    ## Make directory
    beeDir <- paste0("./scratchwild/",i)
    dir.create(beeDir, showWarnings = FALSE)
    ## save file
    mapview::mapshot(leafletMap, file = paste0(beeDir,"/", i, "_", j, ".png"))
    print(paste0(i, "-", j)) ## prints iteration
  }
}

## Need to rename all the pngs to have equal padding of numbers so
## that smaller numbers dont come before biggers ones
## e.g., 3 > 21
allFiles <- list.files("scratchwild", recursive = T, full.names = T)
listLengths <- strsplit(allFiles, "_")
fileID <- do.call(rbind, listLengths)[,2]
fileID_addZeros <- ifelse(nchar(fileID) == 5, paste0("00",fileID), 
                          ifelse(nchar(fileID) == 6, paste0("0",fileID), fileID))
renameFiles <- paste0(do.call(rbind, listLengths)[,1], "_",fileID_addZeros)
file.rename(allFiles, renameFiles)

## function to save as GIF
makeGIFwithPNG <- function(directory, savePath, FPS) {
  require(magick)
  ## list file names and read in
  imgs <- list.files(directory, full.names = TRUE)
  maxImages <- ifelse(length(imgs) > 50, 50, length(imgs))
  imgs <- imgs[1:maxImages] ## cap the max number of images
  img_list <- lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = FPS)
  
  ## save to disk
  image_write(image = img_animated,
              path = savePath)
}

processedBeeMaps <- list.dirs("./scratchwild")[-1]
for(k in 13:length(processedBeeMaps)){
  makeGIFwithPNG(processedBeeMaps[k], paste0(processedBeeMaps[k], ".gif"), 2)
  print(k)
  print(Sys.time())
}

# Making flight paths for each bee
library(leaflet)

allIds <- unique(tracks.out2$id)
for( i in allIds) {
  beeSingular <- tracks.out2%>% filter(id == i)
  beePlot <- ggplot(data = beeSingular, aes(x, y, color = date, label = date)) +
    geom_path(size=0.7) +
    geom_text() +
    scale_color_distiller(palette = "Spectral")
  labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  beePlot
  ggsave( paste0("flightpaths/",i,".pdf"), width = 11, height = 9)
  print(i)
}

##############################################3
#### adding landcover for stats ###############
###########################################

##to do stats with landcover need to add landcover data for true and random points
library(amt)
library(sf)
library(tidyverse)

tracksout.sp <-tracks.out2
tracksout.sp <- st_as_sf(tracksout.sp, coords = c("x", "y"), crs="ESRI:102001")


tracksRandom <- st_sample(perim, size=8876)

rand_sl = random_numbers(make_exp_distr(), n = 1e+05)
rand_ta = random_numbers(make_unif_distr(), n = 1e+05)

getRandomCoords <- function(df, N) { ## needs dataframe and number of random steps per observation
  rand_sl =  random_numbers(fit_distr(df$step, "gamma"), N) ## set random steps based on dataframe
  rand_ta = random_numbers(fit_distr(df$angle, "vonmises"), N) # set random angles based on dataframe
  randomizedCoords <- data.frame() ## empty dataframe to fill with randomized points
  for(i in 1:nrow(df) ){
    x2 <- df[i,"x"] + (rand_sl * cos(rand_ta)) ## calculate change in x-coordinate based on origin, step, and angle
    y2 <- df[i,"y"] + (rand_sl * sin(rand_ta)) ## calculate change in y-coordinate based on origin, step, and angle
    tempDf <- data.frame(id = df[i, "id"], case = "random", ## create dataframe with original bee ID, change in points, and movement info
                         x1= df[i,"x"], y1= df[i,"y"], x2, y2, 
                         sl = rand_sl, ta = rand_ta)
    randomizedCoords <- rbind(randomizedCoords, tempDf) ## output all randomizations per observation
  }
  return(randomizedCoords)
  
}
randomCoords <- getRandomCoords(tracks.out2, 1) ## calculate all random points
plot(randomCoords$x2, randomCoords$y2)

getTrueCoords <- function(df){ ## revise the true steps to match the same as random
  revisedCoordsDF <- df %>% 
    group_by(id) %>%  ## select by bee identifier
    mutate(x2 = lead(x), y2 = lead(y)) %>%  ## create a new column for the end travel point
    filter(!is.na(x2)) %>%  ## drop all initial x-y that don't have an end point
    mutate(case = "true") %>%  ## add a column to identify true
    dplyr::select(case, x1 = x, y1 = y, x2, y2, sl = step, ta = angle) ## match same data structure
  return(revisedCoordsDF)
}
trueCoords <- getTrueCoords(tracks.out2) ## revise data to get all true points
gynetracks.tf <- rbind(randomCoords,trueCoords ) 

#landcover raster needs to be a rasterstack
library(terra)
landStack <- rast()
for(i in 1:7) { 
  tempRaster <- landcover
  tempRaster[tempRaster != i] <- 0
  tempRaster[tempRaster == i] <- 1
  names(tempRaster) <- paste0("landcover",i)
  landStack <- c(landStack, tempRaster)
}

#extract landcover data but first add treatcode
gynetracks.tf <- gynetracks.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

gynetracks.tfs <- gynetracks.tf
gynetracks.tfs<- st_as_sf(gynetracks.tfs, coords = c("x2", "y2"), crs=st_crs(landcover))

# Convert 'gynetracks.tfs' to a SpatVector object to work with terra
gynetracks.tfs <- vect(gynetracks.tfs)

landex<- terra::extract(landStack, gynetracks.tfs)

#add covariates from landStack to the dataframe for SSF
gynetracks.tfs$agriculture <- landex[,1]
gynetracks.tfs$developed <- landex[,2]
gynetracks.tfs$earlyFloral <- landex[,3]
gynetracks.tfs$forest <- landex[,4]
gynetracks.tfs$lateFloral <- landex[,5]
gynetracks.tfs$lowFloral <- landex[,6]
gynetracks.tfs$wetland<- landex[,7]
head(gynetracks.tfs)

gynetracks.tfs$presence <- ifelse(gynetracks.tfs$case == "random", FALSE, TRUE)

####################################
### Step-selection functions ########
#######################################

library(survival)
gynetracks.tfs.NONA <- data.frame(gynetracks.tfs)
gynetracks.tfs.NONA[is.na(gynetracks.tfs.NONA)] <- 0

#maybe try all queens, then break it up by foraging, nest searing and queens

gynetracks.tfs.NONA$presence[gynetracks.tfs.NONA$presence == FALSE] <- 0
gynetracks.tfs.NONA$presence[gynetracks.tfs.NONA$presence == TRUE] <- 1

forageq <- subset(gynetracks.tfs.NONA, treatcode =="F")
nestq <- subset(gynetracks.tfs.NONA, treatcode =="N")
queenq <- subset(gynetracks.tfs.NONA, treatcode =="Q")

SSF1 <- clogit(as.numeric(presence) ~ earlyFloral
               + lateFloral + lowFloral,
               method="approximate", na.action="na.fail", 
               data=gynetracks.tfs.NONA)

SSF2<-clogit(as.numeric(presence) ~ agriculture + developed + earlyFloral + forest
             + lateFloral + lowFloral + wetland,
             method="approximate", na.action="na.fail", 
             data=forageq)

SSF3<-clogit(as.numeric(presence) ~ earlyFloral
             + lateFloral + lowFloral,
             method="approximate", na.action="na.fail", 
             data=nestq)

SSF4<-clogit(as.numeric(presence) ~ agriculture + developed + earlyFloral + forest
             + lateFloral + lowFloral + wetland,
             method="approximate", na.action="na.fail", 
             data=queenq)

library(MuMIn)
D1 <-dredge(SSF1)
D2 <- dredge(SSF2)
D3 <- dredge (SSF3)
D4 <- dredge (SSF4) #not enough occurrences to work. 

SSF3.b <- clogit(as.numeric(presence) ~ agriculture 
                 + lateFloral, 
                 method="approximate", na.action="na.fail", 
                 data=nestq)

###################################
### Step Lengths and turn angles ###
###################################

#using the previous dataset which had random points, those need to be removed

#Separating dataset to make a daytime column
gynetracks.period <- separate(tracks.out2, date, into = c("day","time"), sep=" ")
gynetracks.period <-separate(gynetracks.period, time, sep=":", into=c("hour","minute","sec"))

gynetracks.period[,"daytime"] <- ifelse( as.numeric(gynetracks.period$hour) < 5
                                         | as.numeric(gynetracks.period$hour) > 20, "night", "day")

## adding landcover as categorical variable
gynetracks.ps <- gynetracks.period
gynetracks.ps<- st_as_sf(gynetracks.ps, coords = c("x", "y"), crs=st_crs(landcover))

# Convert 'gynetracks.tfs' to a SpatVector object to work with terra
gynetracks.ps <- vect(gynetracks.ps)

landex<- terra::extract(landcover, gynetracks.ps) #makes a vector of landcover numbers
gynetracks.ps$landtype <- landex$SpringLandcoverRaster
gynetracks.psf<-as.data.frame(gynetracks.ps)
landClasses <- data.frame(landtype = 1:7, landTypeNew = c("agriculture","developed",
                      "earlyFloral", "forest", "lateFloral", "lowFloral", "wetland"))
gynetracks.psf <- gynetracks.psf %>% left_join(landClasses) %>% dplyr::select(-landtype) %>% rename(landtype = landTypeNew)

gynetracks.psf<- gynetracks.psf %>% 
  filter(!is.na(behav)) %>% 
  mutate(ARS = ifelse(behav == "ARS", 1, 0),
         Rest = ifelse(behav == "Rest", 1, 0))


#adding treatcode column
gynetracks.psf2 <- gynetracks.psf %>% 
  mutate(treatcode = substr(id, 0, 1))

library(pastecs)

step.queen <- subset(gynetracks.psf2, treatcode %in% c("F", "Q"))
step.queen$treatcode[step.queen$treatcode == "Q"] <- "F"
#there is only two queens with pollen so collapse into foraging
step.nests <- (subset(gynetracks.psf2, treatcode == "N"))


stat.desc(gynetracks.psf)
stat.desc(abs(gynetracks.psf$angle))
stat.desc(step.queen)
stat.desc(abs(step.queen$angle))
stat.desc(step.nests)
stat.desc(abs(gynetracks.psf$angle))

gynetracks.psf2$treatcode[gynetracks.psf2$treatcode == "Q"] <- "F"


#is there a difference between treatments for step length
library(car)

s1 <- lm(step~treatcode, data=gynetracks.psf2)
car::Anova(s1, type =2) #yes significant

# is there a difference between treatments x time period for step length

s2 <- lm (step ~ treatcode * daytime, data = gynetracks.psf2)
car::Anova(s2, type = 3)
library(emmeans)
s2.a <-emmeans(s2, pairwise ~ treatcode *daytime, data=gynetracks.psf2)
s2.b <- data.frame(s2.a$contrasts) %>% filter(p.value < 0.05)

#step length x treatment x land cover

s3 <- lm (step ~ treatcode * landtype, data= gynetracks.psf2)
car::Anova(s3, type = 3)

#step length x treatment x land cover x time period

s4<- lm (step ~ treatcode * landtype * daytime, data= gynetracks.psf2)
car::Anova(s4, type = 3)
library(emmeans)
s4.a <-emmeans(s4, pairwise ~ treatcode * landtype *daytime, data=gynetracks.psf2)
s4.b <- data.frame(s4.a$contrasts) %>% filter(p.value < 0.05)

##NEW likely do not have enough power for the three way interaction

s5 <- lm(step ~ treatcode * landtype + daytime, data=gynetracks.psf2)
#s5 <- lm(step ~ treatcode * landtype + daytime + 
           #daytime:treatcode, data=gynetracks.psf2) not sure if I care about this
car::Anova(s5, type = 3)
s5.a <-emmeans(s5, pairwise ~ treatcode * landtype, data=gynetracks.psf2)
s5.b <- data.frame(s5.a$contrasts) %>% filter(p.value < 0.05)


# step length x landcover plot

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))
stepplotdat <-gynetracks.psf2 %>% group_by (landtype) %>% 
  summarize (avgStep = mean(step, na.rm=T), errorstep = se(step))

library(ggplot2)
ggplot(data=stepplotdat, aes(x=landtype, y=avgStep)) +
  theme_classic() +
  geom_bar(stat="identity", fill="#292866") +
  geom_errorbar(aes(x = landtype, ymin = avgStep - errorstep, ymax = avgStep + errorstep), 
                width = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 16)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=14)) +
  scale_y_continuous(name= "average step length (m)") +
  scale_x_discrete(limits=c("earlyFloral", "lateFloral", "lowFloral"),
                   labels=c("early floral", "late floral", "low floral"))
#### Turning angle
#is there a difference between treatments for turning angle
library(car)

t1 <- lm(angle~treatcode, data=gynetracks.psf2)
car::Anova(t1, type =2) #significant

# is there a difference between treatments x time period for turning angles

t2 <- lm (angle ~ treatcode * daytime, data = gynetracks.psf2)
car::Anova(t2, type = 3)

#turning angle x treatment x land cover

t3 <- lm (angle ~ treatcode * landtype, data= gynetracks.psf2)
car::Anova(t3, type = 3)
library(emmeans)
t3.a <-emmeans(t3, pairwise ~ treatcode * landtype, data=gynetracks.psf2)
t3.b <- data.frame(t3.a$contrasts) %>% filter(p.value < 0.05)
t3.c <-emmeans(t3, pairwise ~ landtype, data = gynetracks.psf2)

#turning angle x treatment x land cover x time period

t4<- lm (angle ~ treatcode * landtype * daytime, data= gynetracks.psf2)
car::Anova(t4, type = 3)
library(emmeans)
t4.a <-emmeans(t4, pairwise ~ treatcode * landtype *da, data=gynetracks.psf2)
t4.b <- data.frame(t4.a$contrasts) %>% filter(p.value < 0.05)

t5 <- lm(abs(angle) ~ treatcode * landtype + daytime, data=gynetracks.psf2)
#s5 <- lm(step ~ treatcode * landtype + daytime + 
#daytime:treatcode, data=gynetracks.psf2) not sure if I care about this
car::Anova(t5, type = 3)
t5.a <-emmeans(s5, pairwise ~ treatcode * landtype, data=gynetracks.psf2)
t5.b <- data.frame(s5.a$contrasts) %>% filter(p.value < 0.05)

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))
turnplotdat <-gynetracks.psf2 %>% 
  group_by (treatcode, landtype) %>%
  summarize (avgTurn = mean(abs(angle), na.rm=T), errorturn = se(abs(angle)))

zeroValueDF <- data.frame(treatcode = "N", landtype = "earlyFloral", avgTurn = 0, errorTurn = 0)
turnplotdat2 <- rbind(turnplotdat, zeroValueDF)

library(ggplot2)
ggplot(data=turnplotdat2, aes(x=landtype, y=avgTurn, fill = treatcode)) +
  theme_classic() +
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(x = landtype, ymin = avgTurn - errorturn, 
                    ymax = avgTurn + errorturn), 
                width = 0, position = position_dodge(0.9))  +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 16)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))+
  scale_y_continuous(name= "average turn angle") +
  scale_x_discrete(limits=c("earlyFloral", "lateFloral", "lowFloral"),
                   labels=c("early floral", "late floral", "low floral")) +
  scale_fill_manual(values = c("#8070FE", "#0D3601"),
                    labels=c("foraging", "nest-searching"))



### behaviour as a response variable  ###

b1<-glm(ARS ~ treatcode, family="binomial", data=gynetracks.psf2) 
anova(b1, test="Chisq") #significant
b1.a<- emmeans(b1, pairwise ~treatcode, data=gynetracks.psf2)

b2<-glm(ARS~treatcode*daytime, family="binomial", data=gynetracks.psf2)
anova(b2, test="Chisq") #significant
b2.a<-emmeans(b2, pairwise~treatcode *daytime, data=gynetracks.psf2)

b3<- glm(ARS~treatcode * landtype, family="binomial", data=gynetracks.psf2)
anova(b3, test="Chisq")  #significant
b3.a<-emmeans(b3, pairwise~treatcode *landtype, data=gynetracks.psf2)

b4<-glm(ARS~treatcode * landtype * daytime, family="binomial", data=gynetracks.psf2)
anova(b4, test="Chisq") #significant


b5<-glm(Rest~treatcode, family="binomial", data=gynetracks.psf2)
anova(b5, test="Chisq")
# significant

b6<-glm(Rest~treatcode*landtype, family="binomial", data=gynetracks.psf2)
anova(b6, test="Chisq")
b6.a<-emmeans(b6, pairwise~treatcode*landtype, data=gynetracks.psf2)

b7 <- glm(Rest~treatcode * daytime, family = "binomial", data=gynetracks.psf2)
anova(b7, test="Chisq") #not significant

b8 <- glm(Rest ~treatcode*landtype*daytime, family = "binomial", data=gynetracks.psf2)
anova(b8, test = "Chisq") #three way interaction not significant 

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

aggtab <- gynetracks.psf2[,c("treatcode", "behav")]
behavPlot<-aggtab %>% 
  group_by(treatcode) %>% 
  mutate(isRest = behav=="Rest", isARS = behav=="ARS") %>% 
  summarize(rest = sum(isRest)/length(isRest), ARS = sum(isARS)/length(isARS)) %>% 
  tidyr::gather(behav, prop, 2:3) %>% 
  mutate(prop = round(prop,4))

ggplot(behavPlot, aes(x=treatcode, y=prop, fill=behav)) +
  geom_bar (stat="identity", position = position_dodge()) +
  theme_classic () +
  scale_fill_manual(values = c("#59C9A5", "#1C3144"), name = "behaviour") +
  scale_x_discrete(name="", labels = c("control", "cyantraniliprole")) +
  scale_y_continuous(name = "proportion of behaviour") +
  theme (axis.text.x = element_text(size=12),
         axis.title.y  = element_text(size=14),
         axis.text.y = element_text(size =12),
         legend.title = element_text(size = 14),
         legend.text = element_text(size =12))

aggtab2 <- gynetracks.psf2[,c("treatcode", "behav", "landtype")]
behavPlot2<-aggtab2 %>% 
  group_by(treatcode, landtype) %>% 
  mutate(isRest = behav=="Rest", isARS = behav=="ARS") %>% 
  summarize(rest = sum(isRest)/length(isRest), ARS = sum(isARS)/length(isARS)) %>% 
  tidyr::gather(behav, prop, 3:4) %>% 
  mutate(prop = round(prop,4))

labels <- c(CTL = "control", CYN = "cyantraniliprole")
ggplot(behavPlot2, aes(x= landtype, y = prop, fill=behav)) + theme_classic() + 
  facet_grid(treatcode ~ ., labeller = labeller(treatcode = labels)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#59C9A5", "#1C3144"), name = "behaviour") +
  scale_x_discrete(name="", limits=c("agriculture", "forest", "highFloral", "modFloral", "lowFloral"),
                   labels=c("agriculture", "forest", "high floral", "moderate floral", "low floral")) +
  scale_y_continuous(name="proportion of behaviour") +
  theme(strip.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y  = element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size =12))

