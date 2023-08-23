# overwintering queens (aka Sabrina's bees)

#first need to separate out the right bees

allData <- read.csv("TriangulateHandheld.csv") #12245 rows
data <- subset(allData, projectID=="overwintering") #370 rows

library(dplyr)
data2<-data %>% #checking for duplications, no duplicates were found
  distinct() 

write.csv(data, "overwinteringSpring2022.csv")

landcover<- st_read("Spring2022Landcover.shp")
rasterDimX = round(extent(landcover)[2]-extent(landcover)[1], 0 ) ## Rows per raster
rasterDimY = round(extent(landcover)[4]-extent(landcover)[3], 0)  ## columns per raster
empty<-raster(nrows= rasterDimX,  ## create an empty raster with the number of rows based on extent
              ncols = rasterDimY, ## create an empty raster with the number of cols based on extent
              ext = extent(landcover), 
              crs = crs(landcover))
landcover$landclassFactor <- as.factor(landcover$landclass) ## switch classes to numbers alphabetically
rastland <- rasterize(landcover, empty, "landclassFactor") ## convert polygon to raster based on class as a number
rastland2 <- focal(rastland, matrix(1,21,21), fun=modal, NAonly=T,
                   na.rm=T)  ## Take most common land cover within 10 m of NA (for NA values only)

plot(rastland2)
writeRaster(rastland2, "SpringLandcoverRaster.tif", overwrite=T)

#need to remove any points that are outside of the towers detection range
towers<-read.csv("TowerActualLocations.csv")

landcover<-raster("SpringLandcoverRaster.tif") #Rare landcover

coordinates(towers)<-~X+Y
proj4string(towers) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
towersproj<-spTransform(towers, crs(landcover)) #reprojecting CRS to utm matching landcover
towersdf<-as.data.frame(towersproj) 

#making buffer around towers
perim<-buffer(towersproj, width=400, dissolve=T)

#making data spatial points
data.sp <- data
coordinates (data.sp) <- ~ X + Y
proj4string(data.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
data.spj <-spTransform(data.sp, crs(landcover)) #reprojecting CRS to utm

#removing points not in the polygon, using datafrom below
pointsIN<- data.spj[!is.na(over(data.spj, perim)),] #346, 24 obs were outside of tower range
#2018 points before polygon, 1998 within polygon

pointsIN.df<-data.frame(pointsIN)
detsum<-aggregate(pointsIN.df$AnimalID, by=list(pointsIN.df$AnimalID), FUN=length)
detsum #number of observations per bee. 

write.csv(pointsIN.df, "TriangulatTagsPerim.csv")
tagsproj<-pointsIN

library(dplyr)

nums<- pointsIN.df %>% 
  mutate(Var_X = as.numeric(Var_X)) %>% 
  mutate(Var_Y=as.numeric(Var_Y)) %>% 
  filter(!is.na(Var_X)) %>% 
  filter(!is.infinite(Var_X)) %>% 
  filter(!is.na(Var_Y)) %>% 
  filter(!is.infinite(Var_Y))

nums %>% 
  group_by(AnimalID) %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

nums %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

#################################################################
################## Bayes move ##################################
################################################################

library(bayesmove)

######  Pre-specifying Breakpoints and Dealing with Zero-inflated Variables ####

#date is in the wrong format

library(lubridate)
library(tidyr)

tags.df <- pointsIN.df

septags<- tags.df %>% separate(GPSdate, into = c("day", "month", "year"), "_")
# septags <- septags %>% mutate(day = ifelse(as.numeric(day)<10, paste0("0",day), day),
#                                 month = ifelse(as.numeric(month)<10, paste0("0",month), month))
septags$goodDate<-paste0("20", septags$year, "-", septags$month, "-", septags$day)

septags$dtime <- ymd_hms(paste0(septags$goodDate, septags$GPStime))
septagsNoNA <- septags %>% filter(!is.na(dtime))

tagsgdate<-septagsNoNA %>% 
  dplyr::select(X, Y, dtime, AnimalID)
names(tagsgdate)[1] <- "x"
names(tagsgdate)[2] <- "y"
names(tagsgdate)[3] <- "date"
names(tagsgdate)[4] <- "id"
tagsgdate<-tagsgdate %>% arrange(id, date)

### calculating step length, turning angle and time interval
tracks <-prep_data(dat=tagsgdate %>% filter(id != "CTL-20"), coord.names=c("x","y"), id="id") #CTL-20 
#was only observed once so it cannot make a track, removing to continue the code.

head(tracks)
unique (tracks$id) #22 unique track ids 

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

#round time steps to specified interval
hist(tracks.big$dt, breaks=200)
axis(side=1, at=seq(50, 4000, 50)) 

tracks<-round_track_time(dat=tracks, id="id", int=500, tol=3600, time.zone="UTC", units="secs")

#will break up the data into  rest and non-rest

tracks <- tracks %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
#now rest column where 1 = rest and 2 = non-rest

# Create list from data frame where each element is a different track
tracks.list<- df_to_list(dat = tracks, ind = "id")

# Filter observations to time interval
tracks_filt.list<- filter_time(dat.list = tracks.list, int = 500)

#### Discretize data streams (put into bins)

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
angle.bin.lims

dist.bin.lims=quantile(tracks[tracks$dt==500,]$step,
                       c(0, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 1.0), na.rm=T) 
dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list<- map(tracks_filt.list,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))
# Since 0s get lumped into bin 1 for SL, need to add a 8th bin to store only 0s
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

tracks_disc.list3<-tracks_disc.list2 [beesToDrop]

# Pre-specify breakpoints based on 'rest'
breaks<- purrr::map(tracks.list3, ~find_breaks(dat = ., ind = "rest"))

### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(8,8,2)  #SL, TA, rest (in the order from left to right in tracks.list2)

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
nbins<- c(8,8,2)
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
#first 2 states comprise 98% of all observations

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
#most are in the first two behaviour states, a little bit in 3

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
                               "grey35","grey35","grey35","grey35"), guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")

#first behaviour: mainly not resting, some resting. Some zero step length and many at max step length
# all turning angles except directly straight, foraging or area-restricted search
#Second behaviour: Only resting. No step length, facing straight

theta.estim.long<- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
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

## Revise Rest == 1 for all observations
tracks.outTemp <- tracks.out
tracks.outTemp[tracks.outTemp$id %in% beesToRevise &   tracks.outTemp$rest == 1 & tracks.outTemp$dt == 250 &
                 tracks.outTemp$angle > -1 & tracks.outTemp$angle < 1, "behav" ] <- "Rest"
tracks.outTemp <- tracks.outTemp %>% filter(!is.na(behav))

## join processed dataset
tracks.outWithMissing <- tracks.out %>% 
  dplyr::select(names(tracks.outTemp)) %>% 
  rbind(tracks.outTemp)  
tracks.out2<-tracks.outWithMissing

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
  proj4string(dfIn) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m"
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
    beeDir <- paste0("./scratch/",i)
    dir.create(beeDir, showWarnings = FALSE)
    ## save file
    mapview::mapshot(leafletMap, file = paste0(beeDir,"/", i, "_", j, ".png"))
    print(paste0(i, "-", j)) ## prints iteration
  }
}

## Need to rename all the pngs to have equal padding of numbers so
## that smaller numbers dont come before biggers ones
## e.g., 3 > 21
allFiles <- list.files("scratch", recursive = T, full.names = T)
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

processedBeeMaps <- list.dirs("./scratch")[-1]
for(k in 1:length(processedBeeMaps)){
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
library(raster)
library(tidyverse)

tracksout.sp <-tracks.out2
coordinates(tracksout.sp) <-~x+y
proj4string(tracksout.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs" 

tracksRandom <- spsample(tracksout.sp, 486, "random")

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
randomCoords <- getRandomCoords(tracks.out2, 10) ## calculate all random points
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
landStack <- stack()
for(i in 1:7) { 
  tempRaster <- landcover
  tempRaster[tempRaster != i] <- 0
  tempRaster[tempRaster == i] <- 1
  names(tempRaster) <- paste0("landcover",i)
  landStack <- stack(landStack, tempRaster)
}

#extract landcover data but first add treatcode
gynetracks.tf <- gynetracks.tf %>% 
  mutate(treatcode = substr(id, 0, 3))

gynetracks.tfs <- gynetracks.tf
coordinates(gynetracks.tfs) <-~x2+y2
proj4string(gynetracks.tfs) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs" 
landex<- raster::extract(landStack, gynetracks.tfs)

#add covariates from landStack to the dataframe for SSF
gynetracks.tfs$agriculture <- landex[,1]
gynetracks.tfs$developed <- landex[,2]
gynetracks.tfs$forest <- landex[,3]
gynetracks.tfs$highFloral <- landex[,4]
gynetracks.tfs$lowFloral <- landex[,5]
gynetracks.tfs$modFloral <- landex[,6]
gynetracks.tfs$wetland<- landex[,7]
head(gynetracks.tfs)

gynetracks.tfs$presence <- ifelse(gynetracks.tfs$case == "random", FALSE, TRUE)

####################################
### Step-selection functions ########
#######################################

library(survival)
gynetracks.tfs.NONA <- data.frame(gynetracks.tfs)
gynetracks.tfs.NONA[is.na(gynetracks.tfs.NONA)] <- 0

controlNONA <- subset(gynetracks.tfs.NONA, treatcode =="CTL")
cyanNONA <- subset(gynetracks.tfs.NONA, treatcode =="CYN")

SSF1<-clogit(presence ~ agriculture + developed + forest + highFloral + lowFloral 
             + modFloral + wetland,
             method="approximate", na.action="na.fail", 
             data=controlNONA)

SSF2<-clogit(presence ~ agriculture + developed + forest + highFloral + lowFloral 
             + modFloral + wetland,
             method="approximate", na.action="na.fail", 
             data=cyanNONA)

library(MuMIn)
D1 <-dredge(SSF1)
D2 <- dredge(SSF2)

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
coordinates(gynetracks.ps)<-~x+y
proj4string(gynetracks.ps) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

landex<- raster::extract(landcover, gynetracks.ps) #makes a vector of landcover numbers
gynetracks.ps$landtype <- landex
gynetracks.psf<-as.data.frame(gynetracks.ps)
landClasses <- data.frame(landtype = 1:7, landTypeNew = c("agriculture","developed","forest", "highFloral", "modFloral", "lowFloral", "wetland"))
gynetracks.psf <- gynetracks.psf %>% left_join(landClasses) %>% dplyr::select(-landtype) %>% rename(landtype = landTypeNew)

gynetracks.psf<- gynetracks.psf %>% 
  filter(!is.na(behav)) %>% 
  mutate(ARS = ifelse(behav == "ARS", 1, 0),
         Rest = ifelse(behav == "Rest", 1, 0))


#adding treatcode column
gynetracks.psf2 <- gynetracks.psf %>% 
  mutate(treatcode = substr(id, 0, 3))

library(pastecs)

step.control <- (subset(gynetracks.psf2, treatcode == "CTL"))
step.cyan <- (subset(gynetracks.psf2, treatcode == "CYN"))

stat.desc(step.control)
stat.desc(step.cyan)

#is there a difference between treatments for step length
library(car)

s1 <- lm(step~treatcode, data=gynetracks.psf2)
car::Anova(s1, type =2) #yes significant

# is there a difference between treatments x time period for step length

s2 <- lm (step ~ treatcode * daytime, data = gynetracks.psf2)
car::Anova(s2, type = 3)

#step length x treatment x land cover

s3 <- lm (step ~ treatcode * landtype, data= gynetracks.psf2)
car::Anova(s3, type = 3)

#step length x treatment x land cover x time period

s4<- lm (step ~ treatcode * landtype * daytime, data= gynetracks.psf2)
car::Anova(s4, type = 3)
library(emmeans)
s4.a <-emmeans(s4, pairwise ~ treatcode * landtype *daytime, data=gynetracks.psf2)
s4.b <- data.frame(s4.a$contrasts) %>% filter(p.value < 0.05)

# step length x landcover plot

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))
stepplotdat <-gynetracks.psf2 %>% group_by (landtype) %>% 
  summarize (avgStep = mean(step, na.rm=T), errorstep = se(step))

library(ggplot2)
ggplot(data=stepplotdat, aes(x=landtype, y=avgStep)) +
  theme_classic() +
  geom_bar(stat="identity", fill="#858786") +
  geom_errorbar(aes(x = landtype, ymin = avgStep - errorstep+1, ymax = avgStep + errorstep+1), 
                width = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 14)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) +
  scale_y_continuous(name= "average step length") +
  scale_x_discrete(limits=c("agriculture", "forest", "lowFloral", "modFloral", "highFloral"),
                   labels=c("agriculture", "forest", "low floral", "moderate floral", "high floral"))
#### Turning angle
#is there a difference between treatments for turning angle
library(car)

t1 <- lm(angle~treatcode, data=gynetracks.psf2)
car::Anova(t1, type =2) #not significant

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

###########################
### Kernel Density  #######
#############################

library(tidyr)
library(MASS)
library(raster)
library (ggplot2)
library(dplyr)

tagsproj.df <-as.data.frame(tagsproj)

tagsproj.df2 <- tagsproj.df %>% 
  mutate(Treatment = case_when(
    startsWith(AnimalID, "CTL") ~ "control",
    startsWith(AnimalID, "CYN") ~ "cyantraniliprole"))

#Separating dataset to make a daytime column
tagsproj.df3 <- tidyr::separate(tagsproj.df2, GPStime, sep=":", 
                                into=c("hour","minute","sec"))

tagsproj.df3[,"daytime"] <- ifelse( as.numeric(tagsproj.df3$hour) < 5
                           | as.numeric(tagsproj.df3$hour) > 20, "night", "day")


#landcover2 <- aggregate(landcover, fact=5, fun = modal) #reduce the resolution to increase speed
landcover2 <- aggregate(landcover, fact=5, fun = modal)
control <- subset(tagsproj.df3 , Treatment=="control")
controlDay <-subset(control, daytime == "day")
controlNight <- subset(control, daytime == "night")
cyan <- subset(tagsproj.df3, Treatment == "cyantraniliprole")
cyanDay <- subset(cyan, daytime == "day")
cyanNight <-subset(cyan, daytime=="night")

control.kd<-kde2d(x=control$X, control$Y)
controlDay.kd <- kde2d (x=controlDay$X, controlDay$Y)
controlNight.kd <-kde2d (x=controlNight$X, controlNight$Y)
cyan.kd<-kde2d(x=cyan$X, y=cyan$Y)
cyanDay.kd <- kde2d (x=cyanDay$X, cyanDay$Y)
cyanNight.kd <-kde2d (x=cyanNight$X, cyanNight$Y)

control.kdr<-raster(control.kd)
cyan.kdr <- raster(cyan.kd)
controlDay.kdr <- raster(controlDay.kd)
controlNight.kdr <- raster(controlNight.kd)
cyanDay.kdr <- raster(cyanDay.kd)
cyanNight.kdr <- raster(cyanNight.kd)

#extracting landcover variables
landcoverPoly <- rasterToPoints(landcover2)
landpoint<-as.data.frame(landcoverPoly)
coordinates(landpoint)<- ~x + y
proj4string(landpoint)<-crs(landcover)

control.Totalvalues <- raster::extract(control.kdr, landpoint)
control.Totalex<-cbind(control.Totalvalues, data.frame(landpoint))
control.Dayvalues <-raster::extract(controlDay.kdr, landpoint)
control.Dayex <- cbind(control.Dayvalues, data.frame(landpoint))
control.Nightvalues <- raster::extract(controlNight.kdr, landpoint)
control.Nightex <- cbind(control.Nightvalues, data.frame(landpoint))

cyan.Totalvalues <- raster::extract(cyan.kdr, landpoint)
cyan.Totalex<-cbind(cyan.Totalvalues, data.frame(landpoint))
cyan.Dayvalues <-raster::extract(cyanDay.kdr, landpoint)
cyan.Dayex <- cbind(cyan.Dayvalues, data.frame(landpoint))
cyan.Nightvalues <- raster::extract(cyanNight.kdr, landpoint)
cyan.Nightex <- cbind(cyan.Nightvalues, data.frame(landpoint))

treatex <-cbind(control.Totalex, cyan.Totalvalues, control.Dayex, control.Nightex,
               cyan.Dayex, cyan.Nightex)
treatdata <-treatex %>%
  tidyr::gather(treatment, kdensity, c("control.Totalvalues", "cyan.Totalvalues", 
                                "control.Dayvalues", "control.Nightvalues",
                                "cyan.Dayvalues", "cyan.Nightvalues"))
treatdata2 <- treatdata[,c(1:3, ncol(treatdata), ncol(treatdata)-1)] 

write.csv (treatdata2, "KernelDensityData.csv")

##### STOPPING HERE WRITING THE DATA AS CSV TO CONTINUE LATER
kd.data <- read.csv("KernelDensityData.csv")

kd.data2 <- kd.data %>% 
  tidyr::separate(treatment, into=c("treatment","daytime"), sep="[.]")

kd.data3 <-subset(kd.data2, daytime=="Totalvalues") #any stats without daytime

kd1<- lm(kdensity ~ treatment, data=kd.data3)
car::Anova(kd1, type=2)
emmeans(kd1, pairwise~treatment, data=kd.data3)

kd2 <- lm(kdensity ~ treatment * as.factor(SpringLandcoverRaster),
          data = kd.data3)
car::Anova(kd2, type=2)
kd2.a <- emmeans(kd2, pairwise~treatment * as.factor(SpringLandcoverRaster), 
        data=kd.data3)
k2postSig<- data.frame(kd2.a$contrasts) %>% filter(p.value < 0.05)
write.csv(k2postSig, "k2Post.csv", row.names= F)

kd.data4 <- subset(kd.data2, daytime!="Totalvalues") #use this for any daytime stats

kd3 <- lm(kdensity ~ treatment* as.factor(SpringLandcoverRaster) * daytime, 
          data = kd.data4)
car::Anova(kd3, type=2)
kd3.a <- emmeans(kd3, pairwise~treatment * as.factor(SpringLandcoverRaster) *daytime, 
                 data=kd.data4)
k3postSig <- data.frame(kd3.a$contrasts) %>% filter(p.value < 0.05)
write.csv(k3postSig, "k3Post.csv", row.names= F)

kd4 <- lm(kdensity ~ treatment* daytime, data = kd.data4)
car::Anova(kd4, type=2)
kd4.a <- emmeans(kd4, pairwise~treatment *daytime, data=kd.data4)
k4postSig <- data.frame(kd4.a$contrasts) %>% filter(p.value < 0.05)
write.csv(k4postSig, "k4Post.csv", row.names= F)


##plotting results of 3 way interaction

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

kdplotdat <-kd.data4 %>% 
  group_by(SpringLandcoverRaster, treatment, daytime) %>% 
  summarize(meanKD = mean(kdensity*1000, na.rm=T), 
            errorval = se(kdensity*1000))

Timelabels <- c(Dayvalues = "day", Nightvalues = "night")
Treatlabels <- c(control = "control", cyan = "cyantraniliprole")

ggplot(kdplotdat, aes(x = SpringLandcoverRaster, y = meanKD)) + 
  geom_bar(stat = "identity", fill ="#858786") +
  geom_errorbar(aes(x = SpringLandcoverRaster, ymin = meanKD - errorval, 
                    ymax = meanKD + errorval), width = 0) +
  theme_classic() +
  facet_grid(treatment ~ daytime, labeller = labeller (daytime = Timelabels,
                                                       treatment = Treatlabels)) +
  scale_x_discrete(limits=c("1", "2", "4", "3", "5", "6", "7"),
                   labels=c("agriculture", "developed", "forest", "early\nfloral",  
                            "late\nfloral", "low\nfloral",
                            "wetland"), name="") +
  theme(axis.text.x = element_text (size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_y_continuous(name="mean kernel density") 


#need to change crs for plotting with leaflette

control.kdr<-raster(control.kd)
cyan.kdr <- raster(cyan.kd)
controlDay.kdr <- raster(controlDay.kd)
controlNight.kdr <- raster(controlNight.kd)
cyanDay.kdr <- raster(cyanDay.kd)
cyanNight.kdr <- raster(cyanNight.kd)

newproj <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
crs(controlDay.kdr) <-crs(landcover)
crs(controlNight.kdr) <- crs(landcover)
crs(cyanDay.kdr) <-crs(landcover)
crs(cyanNight.kdr) <- crs(landcover)

controlDay.kdr.p<-projectRaster(controlDay.kdr, crs=newproj)
controlNight.kdr.p<-projectRaster(controlNight.kdr, crs=newproj)
cyanDay.kdr.p<-projectRaster(cyanDay.kdr, crs=newproj)
cyanNight.kdr.p<-projectRaster(cyanNight.kdr, crs=newproj)

library(leaflet)
library(colorspace)

#set up colour palette for control
cPal<-sequential_hcl(6, palette = "PuBu", rev=TRUE)
cval = c(0,seq(0, maxValue(controlDay.kdr.p), by = maxValue(controlDay.kdr.p)/5))
cval = cval *100000 #numbers were very small, added this to clean up the legend
cval <- round(cval, 2)
cpal = c("#FFFFFF00", cPal)

#control day
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(controlDay.kdr.p, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")
#control night
cPal<-sequential_hcl(6, palette = "PuBu", rev=TRUE)
cval = c(0,seq(0, maxValue(controlNight.kdr.p), by = maxValue(controlNight.kdr.p)/5))
cval = cval *100000 #numbers were very small, added this to clean up the legend
cval <- round(cval, 2)
cpal = c("#FFFFFF00", cPal)

leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(controlNight.kdr.p, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")
#cyantraniliprole day
cPal<-sequential_hcl(6, palette = "OrYel", rev=TRUE)
cval = c(0,seq(0, maxValue(cyanDay.kdr.p), by = maxValue(cyanDay.kdr.p)/5))
cval = cval *100000 #numbers were very small, added this to clean up the legend
cval <- round(cval, 2)
cpal = c("#FFFFFF00", cPal)

leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(cyanDay.kdr.p, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")
#cyantraniliprole night
cPal<-sequential_hcl(6, palette = "OrYel", rev=TRUE)
cval = c(0,seq(0, maxValue(cyanNight.kdr.p), by = maxValue(cyanNight.kdr.p)/5))
cval = cval *100000 #numbers were very small, added this to clean up the legend
cval <- round(cval, 2)
cpal = c("#FFFFFF00", cPal)

leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(cyanNight.kdr.p, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

####  Minimum Convex Polygon ####

library(adehabitatHR)

tagsproj.df3 <- tagsproj.df3[!is.na(tagsproj.df3$X) 
                               & !is.na(tagsproj.df3$Y),] #remove rows with NA
gynept5 <- tagsproj.df3 %>% 
  group_by(AnimalID) %>%
  mutate(nHits = length(AnimalID)) %>% 
  filter(nHits > 5) #mcp function only works with at least 5 relocations 
gynept5<-gynept5[, c("AnimalID", "X", "Y")] #mcp needs only 3 columns
coordinates(gynept5) <- ~X+Y #specifying coordinates
proj4string(gynept5) <-"+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45" #Lambert Conformal Conic

gyne.mcp <- mcp(gynept5, percent = 100, unout = "km2")
gyne.mcp
names(gyne.mcp)[1]<-"AnimalID"
gyne.mcpDF <- data.frame(gyne.mcp) #dataframe of mcp output with all data columns
str(gyne.mcpDF)

gyne.mcpDF2 <- gyne.mcpDF %>% 
  mutate(Treatment = case_when(
    startsWith(AnimalID, "CTL") ~ "control",
    startsWith(AnimalID, "CYN") ~ "cyantraniliprole"))

m1<-lm(area~Treatment, data=gyne.mcpDF2) #no difference
summary(m1)
anova(m1)

gyne.mcpDF2  %>%
  group_by(Treatment) %>%
  summarise(
    areaAvg = mean(area, na.rm = T),
    sd = sd(area, na.rm=T),
    n = sum(!is.na(area)),
    se = sd / sqrt(n)
  )
##########################################
####  Last Detection Point  ####
#############################################

###calculating the number of days flying

tagsproj.df3$GPSdate <- as.Date(tagsproj.df3$GPSdate, format = "%d_%m_%y")

daysFlown <- tagsproj.df3  %>%
  group_by(AnimalID) %>% 
  summarise(firstDay = (min(GPSdate)), lastDay = (max(GPSdate))) %>% 
  mutate(daysFlown = lastDay - firstDay)

daysFlown.df<-data.frame(daysFlown)

daysFlown.df2<- daysFlown.df %>% 
  mutate(Treatment = case_when(
    startsWith(AnimalID, "CTL") ~ "control",
    startsWith(A, "CYN") ~ "cyantraniliprole"))

daysFlown.df2  %>%
  group_by(Treatment) %>%
  summarise(
    daysAvg = mean(daysFlown, na.rm = T),
    sd = sd(daysFlown, na.rm=T),
    n = sum(!is.na(daysFlown)),
    se = sd / sqrt(n)
  )

write.csv (daysFlown.df2, "DaysFlownOverwintering.csv")
#is there a difference in treatments between days flown?

f1 <- glm(as.numeric(daysFlown) ~ Treatment, family = "poisson", 
          data = daysFlown.df2)
summary(f1) # overdispersed 
library(AER)
dispersiontest(f1) #overdispersed, trying again with negative binomial

library(MASS)
f1.b <- glm.nb(as.numeric(daysFlown) ~ Treatment,  
            data = daysFlown.df2)
summary(f1.b) #possibly still overdispersed 
library(DHARMa) #AER only for poisson so trying this one
simf1.b <- simulateResiduals(f1.b)
testOverdispersion(simf1.b) # no longer overdispersed

anova(f1.b, test = "Chisq")

#extracting the last detection point

tagsproj.df3$GPSdate <-as.POSIXct(as.character(paste(tagsproj.df3$GPSdate, paste(
                                 tagsproj.df3$hour,
                                 tagsproj.df3$min,
                                tagsproj.df3$sec,
                                sep=":"))), format="%Y-%m-%d %H:%M:%S")

lastDay<-tagsproj.df3 %>% 
  group_by(AnimalID) %>% 
  slice(which.max(GPSdate))

lastDay.sp<-lastDay
coordinates(lastDay.sp)<-~X+Y
proj4string(lastDay.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

lastDay.sp$landtype<- raster::extract(landcover, lastDay.sp)

lastDay.df <- as.data.frame(lastDay.sp)

library(ggplot2)

#plot of landcover types associated with last occurrence point
ggplot(data=lastDay.df, aes(x=as.character(landtype), fill=Treatment)) +
  geom_bar(stat="count", position = position_dodge2(preserve="single")) +
  scale_y_continuous(breaks=seq(0, 10, 2), name= "count of last occurrences") +
  scale_x_discrete(limits=c("1", "3", "4", "5", "6"),
                   labels=c("agriculture", "early floral", "forest", 
                            "late floral", "low floral")) +
  scale_fill_manual(values = c("#59C9A5", "#1C3144")) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 12)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) 

#making a map of where the last points are
lastDay.sp<-lastDay
coordinates(lastDay.sp)<-~X+Y
proj4string(lastDay.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"
lastDay.sp <- spTransform(lastDay.sp, "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +towgs84=0,0,0")

set.seed(11) ## makes the plot repeatable
randomAdjustment <- function(nObs, adjRange){
  return(runif(nObs, -0.5, 0.5)/adjRange)
}
## Then run this
lastDay.sp$xAdj <- lastDay.sp$X + randomAdjustment(23, 1000)
lastDay.sp$yAdj <- lastDay.sp$Y + randomAdjustment(23, 1000)

library(leaflet)

treatpal <- colorFactor(palette= c("#59C9A5", "#1C3144"),
                        levels = c("control", "cyantraniliprole"))

leaflet(lastDay.sp) %>% addTiles() %>%
  addProviderTiles('Esri.WorldImagery') %>%
  setView(-80.352, 43.377, zoom = 16) %>%
  addCircleMarkers(~xAdj, ~yAdj, stroke=FALSE, color = ~treatpal(Treatment),
                   fillOpacity=0.9) %>%
  addLegend("topright", pal = treatpal, values= ~Treatment, opacity=1)

### making a land cover plot possibly for the paper
landcover.df <- as.data.frame(landcover, xy=TRUE) %>% 
  na.omit()

landmap <-ggplot(data = landcover.df) +
  geom_raster(aes(x=x, y=y, fill= as.factor(SpringLandcoverRaster))) +
  theme_minimal() +
  scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7"),
                      labels = c("agriculture", "developed", "early floral",
                                 "forest", "late floral", "low floral", "wetland"),
                      values = c("#FFE381", "#C0BCB5", "#827191", "#9BC53D",
                                "#F194B4", "#EFCA08", "#00C2D1")) +
  theme (legend.title = element_blank(),
         axis.title.x = element_text(size =14),
         axis.text.x = element_text(size = 12),
         axis.title.y = element_text(size =14),
         axis.text.y = element_text(size = 12))
landmap
