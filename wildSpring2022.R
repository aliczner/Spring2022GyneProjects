#wild spring queens 2022

#first need to separate out the right bees

allData <- read.csv("TriangulateHandheld.csv") #12245 rows
data <- subset(allData, projectID=="wild") #4727 rows

library(dplyr)
data2<-data %>% #checking for duplications, no duplicates were found
  distinct() 

write.csv(data, "wildSpring2022.csv")
data<-read.csv("wildSpring2022.csv")

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
tags.pj <- st_transform(tags.sf, crs = 3347)
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

tracks<-round_track_time(dat=tracks, id="id", int=20, tol=500, time.zone="UTC", units="secs")

#will break up the data into  rest and non-rest

tracks <- tracks %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
#now rest column where 1 = rest and 2 = non-rest

# Create list from data frame where each element is a different track
tracks.list<- df_to_list(dat = tracks, ind = "id")

# Filter observations to time interval
tracks_filt.list<- filter_time(dat.list = tracks.list, int = 20)

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

# behaviour 1 = resting: only rest,only zero step lengths, and straight turn angles
# behaviour 2 = ARS: not resting med to large step lengths, many turn angles

theta.estim.long <- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
                                   nbehav = 2, behav.names = c("Rest", "ARS"),
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
                             behav.names = c("Rest", "ARS"))

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
