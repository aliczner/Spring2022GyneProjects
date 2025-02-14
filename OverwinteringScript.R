# overwintering queens (aka Sabrina's bees)
library(sf)
library(raster)

#first need to separate out the right bees

allData <- read.csv("TriangulateHandheld.csv") #12245 rows
data <- subset(allData, projectID=="overwintering") #370 rows

library(dplyr)
data2<-data %>% #checking for duplications, no duplicates were found
  distinct() 

write.csv(data, "overwinteringSpring2022.csv")

landcover <- st_read("Spring2022Landcover.shp")
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

landcover <-raster("SpringLandcoverRaster.tif") #Rare landcover

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
tracks <-prep_data(dat=tagsgdate %>% filter(id != "CTL-20"), 
                   coord.names=c("x","y"), id="id") #CTL-20 
#was only observed once so it cannot make a track, removing to continue the code.

head(tracks)
unique (tracks$id) #22 unique track ids 

library(ggplot2)

boxplot(log10(tracks$dt))
tracks.big <- tracks %>% filter(dt < 9000 & dt > 1)

#round time steps to specified interval
hist(tracks.big$dt)

tracks <-round_track_time(dat=tracks, id="id", int=2000, tol=3600, time.zone="UTC", units="secs")
tracks[is.na(tracks)] <- 0

#will break up the data into  rest and non-rest

tracks <- tracks %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
#now rest column where 1 = rest and 2 = non-rest

# Create list from data frame where each element is a different track
tracks.list <- df_to_list(dat = tracks, ind = "id")

# Filter observations to time interval
tracks_filt.list <- filter_time(dat.list = tracks.list, int = 2000)


#### Discretize data streams (put into bins)

angle.bin.lims=round(seq(from=-pi, to=pi, by=pi/4), digits = 3)
angle.bin.lims

dist.bin.lims=quantile(tracks[tracks$dt==2000,]$step,
                       c(0.75, 0.8, 0.85, 0.9, 0.95, 1), na.rm=T) 
dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list <- map(tracks_filt.list,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))

# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.list2 <- map(tracks_disc.list,
                   subset,
                   select = c(id, SL, TA))

## Drop bees without non-rest state (i.e., no "2"s)
#beesToDrop <- sapply(tracks.list2, function(x){
# ifelse(sum(x$rest) == nrow(x) *2 || sum(x$rest) == nrow(x), FALSE, TRUE)
#})
#tracks.list3 <- tracks.list2[beesToDrop]

#tracks_disc.list3<-tracks_disc.list2 [beesToDrop]

### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(5,8)  #SL, TA (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 4)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)

future::plan(future::sequential)  #return to single core

## Determine MAP for selecting breakpoints (maximum a posteriori)
MAP.est <- get_MAP(dat = dat.res$LML, nburn = 5000)
MAP.est

brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))

tracks.seg<- assign_tseg(dat = tracks_disc.list, brkpts = brkpts)
head(tracks.seg)

# Select only id, tseg, SL, TA, and rest columns
tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]

# Summarize observations by track segment
nbins<- c(5,8)
obs<- summarize_tsegs(dat = tracks.seg2, nbins = nbins)
head(obs)

#### Running the clustering model

set.seed(1)

# Prepare for Gibbs sampler
ngibbs <- 25000  #number of MCMC iterations for Gibbs sampler
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


### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim <- extract_prop(res = res, ngibbs = ngibbs, nburn = nburn, 
                            nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means<- round(colMeans(theta.estim), digits = 3))

# Calculate cumulative sum
cumsum(theta.means)
#first 2 states comprise 99% of all observations

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
#most are in the first two behaviour states

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

#first behaviour: mainly resting. Some greater step lengths
# mainly straight but some at the extremes
# could be resting with some foraging or nest searching?
#Second behaviour: Only moving. largest step lengths tortuous movement
#transit

theta.estim.long <- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
                                   nbehav = 2, behav.names = c("ARS", "Rest"),
                                   behav.order = c(1,2))

##### Assign behavioural states to tracks

# Merge results with original data
tracks.out<- assign_behavior(dat.orig = tracks,
                             dat.seg.list = df_to_list(tracks.seg, "id"),
                             theta.estim.long = theta.estim.long,
                             behav.names = c("ARS", "Rest"))


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
allIds <- unique(tracks.out$id)
for( i in allIds) {
  beeSingular <- tracks.out %>% filter(id == i) %>% UTMtoLatLon(.) 
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
for(k in 21:length(processedBeeMaps)){
  makeGIFwithPNG(processedBeeMaps[k], paste0(processedBeeMaps[k], ".gif"), 2)
  print(k)
  print(Sys.time())
}

#making flight path simple maps without leaflet for paper

library(ggplot2)
library(dplyr)

tests <- tracks.out %>%
  group_by(id) %>%
  mutate(row_number = row_number())


allIds <- unique(tests$id)
for(i in 1:length(allIds)){
  maxRow <- max(tests[tests$id == allIds[i],"row_number"] )
  ggplot(tests %>% filter(id == allIds[i]), aes(x = x, y = y, color = row_number)) +
    geom_path(lwd = 2) +
    labs(title = paste0(allIds[i]), x = "", y = "") +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    guides(color = guide_legend(title = "Step number")) +
    scale_color_viridis_c(option = "C", limits = c(1, maxRow), breaks = round(seq(1, maxRow, length.out = 11), 0)) 
  ggsave(paste0("IndividualBeeFlights/BeePath_",allIds[i],".png"))
}

##############################################3
#### adding landcover for stats ###############
###########################################

##to do stats with landcover need to add landcover data for true and random points
library(amt)
library(raster)
library(tidyverse)

tracksout.flight <- subset(tracks.out, dt == 2000)
tracksout.sp <-tracksout.flight
coordinates(tracksout.sp) <-~x+y
proj4string(tracksout.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs" 

tracksRandom <- spsample(tracksout.sp, 118, "random")

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
randomCoords <- getRandomCoords(tracksout.flight, 1) ## calculate all random points
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
trueCoords <- getTrueCoords(tracksout.flight) ## revise data to get all true points
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
gynetracks.tfs$earlyFloral <- landex[,3]
gynetracks.tfs$forest <- landex[,4]
gynetracks.tfs$lateFloral <- landex[,5]
gynetracks.tfs$lowFloral<- landex[,6]
gynetracks.tfs$wetland<- landex[,7]
head(gynetracks.tfs)

gynetracks.tfs$presence <- ifelse(gynetracks.tfs$case == "random", FALSE, TRUE)

####################################
### Step-selection functions ########
#######################################

library(survival)
gynetracks.tfs.NONA <- data.frame(gynetracks.tfs)
gynetracks.tfs.NONA[is.na(gynetracks.tfs.NONA)] <- 0

## trying to keep the datasets together

SSF.treats <- glm(presence ~ treatcode * (agriculture + developed + forest + 
                                            earlyFloral + lateFloral + lowFloral + 
                                            wetland), 
                  family = "binomial", data = gynetracks.tfs.NONA)
summary(SSF.treats)

#removing wetland because it has 0 obs and developed which has  3 obs

SSF.treats2 <- glm(presence ~ treatcode * (agriculture + forest + earlyFloral + 
                                             lateFloral + lowFloral), 
                   family = "binomial", data = gynetracks.tfs.NONA)
summary(SSF.treats2)

#nothing significant.

controlNONA <- subset(gynetracks.tfs.NONA, treatcode =="CTL")
cyanNONA <- subset(gynetracks.tfs.NONA, treatcode =="CYN")

SSF1 <-glm(presence ~ agriculture + developed + forest + earlyFloral + lowFloral 
             + lateFloral + wetland, family = "binomial", data=controlNONA)
summary(SSF1) 
#there are no observations in wetland, and very few in developed and ag, try removing
#just those two

SSF2 <- glm(presence ~ forest + earlyFloral +lowFloral + lateFloral,
            family = "binomial", data = controlNONA)
summary(SSF2) # nothing is significant

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF2, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF2)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

ggplot(data = odds_ratio.df2, aes(x = landtype, y = oddsRatio)) +
  theme_classic() +
  geom_point(size = 4, colour = c("#235789","#235789","#235789", "#235789")) +
  geom_errorbar(aes( x = landtype , ymin = lowci, 
                     ymax = highci), width = 0, 
                colour = c("#235789", "#235789","#235789", "#235789")) +  
  theme(axis.title.x = element_text (size = 14), 
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(limits = c("earlyFloral", "forest", "lateFloral", 
                              "lowFloral"),
                   labels = c("forest", "early floral", "late floral",
                              "low floral")) +
  scale_y_continuous(name = "odds ratio") +
  geom_hline(yintercept = 1, linetype = "solid", colour = "dark grey") +
  coord_flip() +
  ggtitle("control") +
  theme(plot.title = element_text(size = 16))


SSF4 <-glm(presence ~forest + earlyFloral + lowFloral 
             + lateFloral , family = "binomial", data=cyanNONA)
summary(SSF4) # late floral not modeling due to singularities ?

SSF5 <-glm(presence ~ forest + earlyFloral + lowFloral, family = "binomial", 
           data=cyanNONA)
summary(SSF5)

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF5, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF5)))

odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF5)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

ggplot(data = odds_ratio.df2, aes(x = landtype, y = oddsRatio)) +
  theme_classic() +
  geom_point(size = 4, colour = c("#8C271E", "#8C271E","black")) +
  geom_errorbar(aes( x = landtype , ymin = lowci, 
                     ymax = highci), width = 0, colour = c("#8C271E", "#8C271E",
                                                           "black")) +  
  theme(axis.title.x = element_text (size = 14), 
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(limits = c("forest", "earlyFloral", "lowFloral"),
                   labels = c("forest", "early floral","low floral")) +
  scale_y_continuous(name = "odds ratio") +
  geom_hline(yintercept = 1, linetype = "solid", colour = "dark grey") +
  coord_flip() +
  ggtitle("cyantraniliprole") +
  theme(plot.title = element_text(size = 16))

###################################
### Step Lengths and turn angles ###
###################################

#using the previous dataset which had random points, those need to be removed

#Separating dataset to make a daytime column
gynetracks.period <- tidyr::separate(tracks.out, date, into = c("day","time"), sep=" ")
gynetracks.period <-tidyr::separate(gynetracks.period, time, sep=":", into=c("hour","minute","sec"))

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
landClasses <- data.frame(landtype = 1:7, 
                          landTypeNew = c("agriculture","developed","earlyFloral", 
                          "forest", "lateFloral", "lowFloral", "wetland"))
gynetracks.psf <- gynetracks.psf %>% 
  left_join(landClasses) %>% 
  dplyr::select(-landtype) %>% 
  rename(landtype = landTypeNew)

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
car::Anova(s1, type =2) #no not significant
pwr::pwr.f2.test(u = 1, v = NULL, f2 = 0.009 , sig.level = 0.05, power = 0.8)
#not enough power to test interactions

# possibly individual factors

s2 <- lm (step ~ treatcode + daytime + landtype, data = gynetracks.psf2)
car::Anova(s2, type = 2) #none are significant. 
#power is still low

#### Turning angle
t1 <- lm(angle~treatcode, data=gynetracks.psf2)
car::Anova(t1, type =2) #not significant
pwr::pwr.f2.test(u = 1, v = 281, f2 = 0.0009 , sig.level = 0.05, power = NULL)
#power is very low 0.08

t2 <- lm (angle ~ treatcode + daytime + landtype, data = gynetracks.psf2)
car::Anova(t2, type = 2) #no significant differnce
pwr::pwr.f2.test(u = 1, v = 276, f2 = 0.03 , sig.level = 0.05, power = NULL)
#power is good

####  Net squared displacement NSD

gynetracks.nsd <-subset(gynetracks.psf2, dt == 2000)

n1 <- lm(NSD ~ treatcode, data=gynetracks.nsd)
car::Anova(n1, type =2) #significant
pwr::pwr.f2.test(u = 1, v = 116, f2 = 0.05 , sig.level = 0.05, power = NULL)
#power isn't great 0.67

## NSX x treatment x daytime 
n2 <- lm (NSD ~ treatcode * daytime, data = gynetracks.nsd)
car::Anova(n2, type = 3)# significant
n2.a <- emmeans::emmeans(n2, pairwise ~ treatcode * daytime, data = gynetracks.nsd)
pwr::pwr.f2.test(u = 1, v = 114, f2 = 0.08 , sig.level = 0.05, power = NULL)
#power is good 0.86

# NSD x treatment x land cover

n3 <- lm (NSD ~ treatcode * landtype, data= gynetracks.nsd)
car::Anova(n3, type = 3)
n3.a <- emmeans::emmeans(n3, pairwise ~ treatcode * landtype, data = gynetracks.nsd)

#figure for NSD by treatment and daytime
se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

NSDplot1 <- gynetracks.nsd[,c("treatcode", "NSD", "daytime")] %>% 
  group_by(daytime, treatcode) %>% 
  summarise(avgNSD = mean(NSD), errorNSD = se(NSD))

ggplot(NSDplot1, aes(x=treatcode, y=(avgNSD/1000), fill=daytime)) +
  geom_bar (stat="identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = avgNSD/1000 - errorNSD/1000, ymax = avgNSD/1000 + errorNSD/1000),
                width = 0, position = position_dodge(width=0.9)) +
  theme_classic () +
  scale_fill_manual(values = c("#78C0E0", "#27187E"), name = "time of day") +
  scale_x_discrete(limits=c("CTL", "CYN"),
                   labels=c("control", "cyantraniliprole")) +
  scale_y_continuous(name = "mean squared displacement (km)", 
                     labels = scales::comma) +
  theme (axis.text.x = element_text(size=12),
         axis.title.x = element_blank(), 
         axis.title.y  = element_text(size=14),
         axis.text.y = element_text(size =12),
         legend.title = element_text(size = 14),
         legend.text = element_text(size =12))

# figure for NSD and land cover
se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

NSDplot2 <- gynetracks.nsd[,c("treatcode", "NSD", "landtype")] %>% 
  group_by(landtype, treatcode) %>% 
  summarize(avgNSD = mean(NSD), errorNSD = se(NSD))

missingData <- data.frame(landtype = c("agriculture","forest"), treatcode =c("CYN","CYN"), avgNSD=c(0,0), errorNSD=c(0,0))
allData <- rbind(NSDplot2, missingData)

ggplot(allData, aes(x=landtype, y=avgNSD/1000, fill=treatcode)) +
  geom_bar (stat="identity", position = position_dodge(preserve ="single")) +
  geom_errorbar(aes(ymin = avgNSD/1000 - errorNSD/1000, ymax = avgNSD/1000 + errorNSD/1000),
                width = 0, position = position_dodge(width=0.9)) +
  theme_classic () +
  scale_fill_manual(values = c("#20BF55", "#8D6A9F"), name = "",
                    breaks = c("CTL", "CYN"),
                    labels = c("control", "cyantraniliprole")) +
  scale_x_discrete(limits=c("agriculture", "forest", "earlyFloral", "lateFloral", 
                            "lowFloral"),
                   labels=c("agriculture", "forest", "early floral", "late floral", 
                            "low floral")) +
  scale_y_continuous(name = "mean squared displacement (km)",
                     labels = scales::comma) +
  theme (axis.text.x = element_text(size=12),
         axis.title.x = element_blank(), 
         axis.title.y  = element_text(size=14),
         axis.text.y = element_text(size =12),
         legend.title = element_text(size = 14),
         legend.text = element_text(size =12))


###########################################
### behaviour as a response variable  ###
###############################################

#updated to betaregression
b1 <- betareg::betareg(Rest ~ treatcode, data = na.omit(gynetracks.psf2))
lmtest::lrtest(b1) #significant
pwr::pwr.f2.test(u = 1, v = 5, f2 = 0.09 , sig.level = 0.05, power = NULL)

b2<-betareg::betareg(Rest ~ treatcode + daytime + landtype, data=gynetracks.psf2)
lmtest::lrtest(b2) #not significant


se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

behavplot1 <- gynetracks.psf2 %>%
  select(ARS, Rest, landtype, treatcode) %>%
  gather(behav, value, ARS:Rest) %>%
  group_by(treatcode, behav) %>%
  summarize(means = mean(value, na.rm=T), ses = se(value))


ggplot(behavplot1 , aes(x=behav, y = means, fill = treatcode)) +
  geom_bar (stat="identity", 
            position = position_dodge()) +
geom_errorbar(aes(ymin = means - ses, ymax = means + ses), 
              position = position_dodge(width=0.8), width = 0) +
  theme_classic () +
  scale_fill_manual(values = c("#20BF55", "#8D6A9F"), name = "", 
                    breaks = c("CTL", "CYN"),
                    labels = c ("control", "cyantranilprole")) +
  scale_y_continuous(name = "mean proportion of behaviour") +
  theme (axis.text.x = element_text(size=12),
         axis.title.x = element_blank(), 
         axis.title.y  = element_text(size=14),
         axis.text.y = element_text(size =12),
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

kd4 <- lm(kdensity ~ treatment* daytime, data = kd.data4)
car::Anova(kd4, type=2)
kd4.a <- emmeans(kd4, pairwise~treatment *daytime, data=kd.data4)
k4postSig <- data.frame(kd4.a$contrasts) %>% filter(p.value < 0.05)
write.csv(k4postSig, "k4Post.csv", row.names= F)


### plotting results for kernel density
# one plot for treatment x land cover
# another plot for treatment x daytime 

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

kdplotdat1 <-kd.data4 %>% 
  group_by(SpringLandcoverRaster, treatment) %>% 
  summarize(meanKD = mean(kdensity*100000, na.rm=T), 
            errorval = se(kdensity*100000))

ggplot(kdplotdat1, aes(x = SpringLandcoverRaster, y = meanKD, fill = treatment)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(x = SpringLandcoverRaster, ymin = meanKD - errorval, 
                    ymax = meanKD + errorval), position = position_dodge(0.9),
                width = 0) +
  theme_classic() +
  scale_x_discrete(limits=c("1", "2", "4", "3", "5", "6", "7"),
                   labels=c("agriculture", "developed", "forest", "early\nfloral",  
                            "late\nfloral", "low\nfloral",
                            "wetland"), name="") +
  theme(axis.text.x = element_text (size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.text = element_text(size = 12)) +
  scale_y_continuous(name="mean kernel density") +
  scale_fill_manual(values = c("#20BF55", "#8D6A9F"),
                    name = "", 
                      breaks = c("control", "cyan"), 
                      labels = c("control", "cyantraniliprole"))
# dayplot 

kdplotdat2 <-kd.data4 %>% 
  group_by(daytime, treatment) %>% 
  summarize(meanKD = mean(kdensity*100000, na.rm=T), 
            errorval = se(kdensity*100000))

ggplot(kdplotdat2, aes(x = treatment, y = meanKD, fill = daytime)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(x = treatment, ymin = meanKD - errorval, 
                    ymax = meanKD + errorval), position = position_dodge(0.9),
                width = 0) +
  theme_classic() +
  scale_x_discrete(limits=c("control", "cyan"), 
                   labels = c("control", "cyantraniliprole"), name="") +
  theme(axis.text.x = element_text (size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        legend.text = element_text(size = 12)) +
  scale_y_continuous(name="mean kernel density") +
  scale_fill_manual(values = c("#78C0E0", "#27187E"),
                    name = "", 
                    breaks = c("Dayvalues", "Nightvalues"), 
                    labels = c("day", "night"))


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
    startsWith(AnimalID, "CYN") ~ "cyantraniliprole"))

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
testDispersion(simf1.b) # still overdispersed

f1.c <- glm(as.numeric(daysFlown) ~ Treatment, family="quasipoisson", 
            data = daysFlown.df2)
summary(f1.c) ## still overdispersed

# could try this bayes model later if you want but need to install
#the package, didn't want to do it as it needed to restart. 

library(rstanarm)

# Convert response variable to integer type if it's not already
daysFlown.df2$daysFlown <- as.integer(daysFlown.df2$daysFlown)

# Fit a Bayesian Poisson regression model
bayesian_model <- stan_glm(daysFlown ~ Treatment, 
                           family = poisson(), 
                           data = daysFlown.df2,
                           chains = 4, # Number of Markov chains
                           iter = 2000) # Number of iterations

# Summary of the Bayesian model
summary(bayesian_model)




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


