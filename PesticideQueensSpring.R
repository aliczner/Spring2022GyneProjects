## 2022 spring pesticide experiment with cyan, flupyridifurone, imidacloprid

#first need to separate out the right bees

allData <- read.csv("TriangulateHandheld.csv") #12245 rows
data <- subset(allData, projectID =="pesticide") #7144 rows
data.nomotus <-subset(data, dataSource != "motus")

library(dplyr)
data2 <-data %>% #checking for duplications, no duplicates were found
  distinct() 

#need to remove any points that are outside of the towers detection range
towers<-read.csv("TowerActualLocations.csv")

landcover <-terra::rast("SpringLandcoverRaster.tif") #Rare landcover, made in overwintering

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

pointsIN.df <-data.frame(result_sf) #970
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
library(lubridate)
library(tidyr)
library(sf)

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
tags.latlong <- st_coordinates(tags.pj)

tags.pj2 <- cbind(tags.pj, tags.latlong)
tags.pj2df <-as.data.frame(tags.pj2)
tags.dat<-tags.pj2df %>% 
  dplyr::select(date, id, X, Y)

# there are multiple id's with only one observation. that won't work with bayesmove needs to be
# removed. 

tags.dat2 <- tags.dat %>% 
  group_by(id) %>% 
  filter(n()>1) %>% 
  ungroup
# removed ids = control:1, 10, 11, 14, imid: 1, 3

### back to the bayes move 
#calculating step length, turning angle and time interval
tracks <- prep_data(dat = tags.dat2, coord.names = c("X", "Y"),  id = "id")
tracks <- replace (tracks, is.na(tracks), 0)
#can't run prep_data on any ids with just one occurrence.

head(tracks)
unique (tracks$id) #37 unique tracks

library(ggplot2)

boxplot(log10(tracks$dt))
tracks.big <- tracks %>% filter(dt < 10000 & dt > 1)
tracks.big2 <- tracks.big %>% filter(dt < 3600 & dt > 1)
#14 at 1 hour
#31 at 24 hrs 
#36 at one week



#round time steps to the three intervals
hist(tracks.big2$dt, breaks=100)

tracks.hour <-round_track_time(dat=tracks, id="id", int=1800, tol=1800, 
                               time.zone="UTC", units="secs")

tracks.hour <- tracks.hour %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
tracks.day <-round_track_time(dat=tracks, id="id", int=43200, tol=43200, 
                               time.zone="UTC", units="secs")

tracks.day <- tracks.day %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
tracks.week <-round_track_time(dat=tracks, id="id", int=302400, tol=302400, 
                              time.zone="UTC", units="secs")

tracks.week <- tracks.week %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
# Create list from data frame where each element is a different track
tracks.list.hour <- df_to_list(dat = tracks.hour, ind = "id")
tracks.list.day <- df_to_list(dat = tracks.day, ind = "id")
tracks.list.week <- df_to_list(dat = tracks.week, ind = "id")

# Filter observations to time interval
tracks_filt.list.hour <- filter_time(dat.list = tracks.list.hour, int = 1800)
tracks_filt.list.day <- filter_time(dat.list = tracks.list.day, int = 43200)
tracks_filt.list.week <- filter_time(dat.list = tracks.list.week, int = 302400)

#### Discretize data streams (put into bins)
angle.bin.lims=c(-3.143, -2.356, -1.571, -0.785,  0, 0.785,  1.571,  2.356,  3.143)
angle.bin.lims

tracks.hour_norest <- subset(tracks.hour, rest == 2)
dist.bin.lims=quantile(tracks.hour_norest[tracks.hour_norest$dt==1800,]$step,
                       c(0, 0.1, 0.25, 0.30, 0.50, 0.75, 0.90, 1), na.rm=T) #7 bins

tracks.day_norest <- subset(tracks.day, rest == 2)
dist.bin.lims=quantile(tracks.day_norest[tracks.day_norest$dt==43200,]$step,
                       c(0, 0.1, 0.25, 0.30, 0.50, 0.75, 0.90, 1), na.rm=T) #7 bins

tracks.week_norest <- subset(tracks.week, rest == 2)
dist.bin.lims=quantile(tracks.week_norest[tracks.week_norest$dt==302400,]$step,
                       c(0, 0.1, 0.25, 0.30, 0.50, 0.75, 0.90, 1), na.rm=T) #7 bins

dist.bin.lims
dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list.hour <- map(tracks_filt.list.hour,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))

tracks_disc.list.day <- map(tracks_filt.list.day,
                            discrete_move_var,
                            lims = list(dist.bin.lims, angle.bin.lims),
                            varIn = c("step", "angle"),
                            varOut = c("SL", "TA"))

tracks_disc.list.week <- map(tracks_filt.list.week,
                           discrete_move_var,
                           lims = list(dist.bin.lims, angle.bin.lims),
                           varIn = c("step", "angle"),
                           varOut = c("SL", "TA"))

# Since 0s get lumped into bin 1 for SL, need to add an extra bin to store 0s
tracks_disc.list.hour2 <- tracks_disc.list.hour %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))

tracks_disc.list.day2 <- tracks_disc.list.day %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))

tracks_disc.list.week2 <- tracks_disc.list.week %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))
# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.hour.list2 <- map(tracks_disc.list.hour2,
                   subset,
                   select = c(id, SL, TA))

tracks.day.list2 <- map(tracks_disc.list.day2,
                        subset,
                        select = c(id, SL, TA))

tracks.week.list2 <- map(tracks_disc.list.week2,
                        subset,
                        select = c(id, SL, TA))

### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(8,8)  #SL, TA, rest (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.hour.res <- segment_behavior(data = tracks.hour.list2, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)

dat.day.res <- segment_behavior(data = tracks.day.list2, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha)


dat.week.res <- segment_behavior(data = tracks.week.list2, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha)

future::plan(future::sequential)  #return to single core

## Determine MAP for selecting breakpoints (maximum a posteriori)
get_MAP <- function(dat, nburn) {
  
  # Subset only numeric columns after burn-in and convert to matrix
  tmp <- as.matrix(dat[, (nburn + 2):ncol(dat)])  # Subset columns after burn-in and ensure matrix format
  
  # Find max LML per row (ID) after burn-in
  MAP.est <- as.integer(apply(tmp, 1, function(x) which.max(x)) + nburn)
  
  return(MAP.est)
}

MAP.est.hour <- get_MAP(dat = dat.hour.res$LML[,-1], nburn = 5000)
MAP.est.hour
idColumn <- dat.hour.res$LML[,1]

MAP.est.day <- get_MAP(dat = dat.day.res$LML[,-1], nburn = 5000)
MAP.est.day
idColumn <- dat.day.res$LML[,1]

MAP.est.week <- get_MAP(dat = dat.week.res$LML[,-1], nburn = 5000)
MAP.est.week
idColumn <- dat.week.res$LML[,1]

brkpts.hour <- get_breakpts(dat = dat.hour.res$brkpts, MAP.est = MAP.est.hour)
brkpts.day <- get_breakpts(dat = dat.day.res$brkpts, MAP.est = MAP.est.day)
brkpts.week <- get_breakpts(dat = dat.week.res$brkpts, MAP.est = MAP.est.week)

# How many breakpoints estimated per ID?
apply(brkpts.hour[,-1], 1, function(x) length(purrr::discard(x, is.na)))
apply(brkpts.day[,-1], 1, function(x) length(purrr::discard(x, is.na)))
apply(brkpts.week[,-1], 1, function(x) length(purrr::discard(x, is.na)))

tracks.seg.hour <- assign_tseg(dat = tracks_disc.list.hour2, 
                              brkpts = brkpts.hour)
head(tracks.seg.hour)

tracks.seg.day <- assign_tseg(dat = tracks_disc.list.day2, 
                              brkpts = brkpts.day)
head(tracks.seg.day)

tracks.seg.week <- assign_tseg(dat = tracks_disc.list.week2, 
                              brkpts = brkpts.week)
head(tracks.seg.week)

# Select only id, tseg, SL, TA,
tracks.seg.hour2 <- tracks.seg.hour[,c("id","tseg","SL","TA")]
tracks.seg.day2 <- tracks.seg.day[,c("id","tseg","SL","TA")]
tracks.seg.week2 <- tracks.seg.week[,c("id","tseg","SL","TA")]

# Summarize observations by track segment
nbins <- c(8,8)

obs <- summarize_tsegs(dat = tracks.seg.hour2, nbins = nbins)
obs <- summarize_tsegs(dat = tracks.seg.day2, nbins = nbins)
obs <- summarize_tsegs(dat = tracks.seg.week2, nbins = nbins)
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
res.hour <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                        ngibbs=ngibbs, nmaxclust=nmaxclust,
                        nburn=nburn, ndata.types=ndata.types)

res.day <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                             ngibbs=ngibbs, nmaxclust=nmaxclust,
                             nburn=nburn, ndata.types=ndata.types)

res.week <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                            ngibbs=ngibbs, nmaxclust=nmaxclust,
                            nburn=nburn, ndata.types=ndata.types)



### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim.hour <- extract_prop(res = res.hour, ngibbs = ngibbs, nburn = nburn,
                            nmaxclust = nmaxclust)

theta.estim.day <- extract_prop(res = res.day, ngibbs = ngibbs, nburn = nburn,
                                 nmaxclust = nmaxclust)

theta.estim.week <- extract_prop(res = res.week, ngibbs = ngibbs, nburn = nburn,
                                nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means.hour <- round(colMeans(theta.estim.hour), digits = 3))
(theta.means.day <- round(colMeans(theta.estim.day), digits = 3))
(theta.means.week <- round(colMeans(theta.estim.week), digits = 3))

# Calculate cumulative sum
cumsum(theta.means.hour)
cumsum(theta.means.day)
cumsum(theta.means.week)
#first 2 states comprise 98.7% of all observations

# Convert to data frame for ggplot2
theta.estim_df.hour<- theta.estim.hour %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", 
               values_to = "prop") %>%
  modify_at("behaviour", factor)
levels(theta.estim_df.hour$behaviour)<- 1:nmaxclust

theta.estim_df.day<- theta.estim.day %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", 
               values_to = "prop") %>%
  modify_at("behaviour", factor)
levels(theta.estim_df.day$behaviour)<- 1:nmaxclust

theta.estim_df.week<- theta.estim.week %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", 
               values_to = "prop") %>%
  modify_at("behaviour", factor)
levels(theta.estim_df.week$behaviour)<- 1:nmaxclust


# Plot results
ggplot(theta.estim_df.hour, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in the first two behaviour states

ggplot(theta.estim_df.day, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in 2 behaviour states, a little in a third

ggplot(theta.estim_df.week, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in 2 behaviour states, very few in a third

#### Classify the states as behaviours

# Extract bin estimates from phi matrix
behav.res.hour <- get_behav_hist(dat = res.hour, nburn = nburn, 
                          ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle"))
behav.res.day <- get_behav_hist(dat = res.day, nburn = nburn, 
                                 ngibbs = ngibbs, nmaxclust = nmaxclust,
                                 var.names = c("Step Length","Turning Angle"))
behav.res.week <- get_behav_hist(dat = res.week, nburn = nburn, 
                                ngibbs = ngibbs, nmaxclust = nmaxclust,
                                var.names = c("Step Length","Turning Angle"))


# Plot histograms of proportion data
ggplot(behav.res.hour, aes(x = bin, y = prop, fill = as.factor(behav))) +
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
  scale_x_continuous(breaks = 1:7) +
  facet_grid(behav ~ var, scales = "free_x")

# behaviour 1 = ARS or foraging. mainly zero SL mainly straight but has 
                #values for both SL and TA that span the gradient
# behaviour 2 = exploratory/searching. big SL and tortuous
# behaviour 3 = rest. SL = 0 TA = straight.

theta.estim.long.hour <- expand_behavior(dat = tracks.seg.hour, theta.estim = theta.estim.hour, obs = obs,
                                    nbehav = 3, behav.names = c("ARS", "Explore", "Rest"),
                                    behav.order = c(1,2,3))

ggplot(behav.res.day, aes(x = bin, y = prop, fill = as.factor(behav))) +
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
  scale_x_continuous(breaks = 1:7) +
  facet_grid(behav ~ var, scales = "free_x")
#1st behav: mainly 0 SL but some larger ones, mainly straight, some other TA
#2nd behav: mainly large SL, mainly tortuous TA
#3rd behav: 0 SL and straight TA

#behav 1 = localized use of resourced-focused station keeping. animal is primarly staying
#in a favourable area with minimal exploration outside of it. different from ARS in that
#ARS is more systematic searching behaviour. so it might show more tortuous turns, and
#fewer zero steps as the animal actively moves around looking for resources. 
#behav 2 = broad-scale ARS because of larger SL, resources may be patchy or spread out
#behav 3 = rest

theta.estim.long.day <- expand_behavior(dat = tracks.seg.day, theta.estim = theta.estim.day, 
                                         obs = obs,
                                         nbehav = 3, 
                                         behav.names = c("StationKeeping", "BroadARS", "Rest"),
                                         behav.order = c(1,2,3))

# week 
ggplot(behav.res.week, aes(x = bin, y = prop, fill = as.factor(behav))) +
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
  scale_x_continuous(breaks = 1:7) +
  facet_grid(behav ~ var, scales = "free_x")
#1st behav: all different SL and TA
#2nd behav: zero SL and striaght TA
#3rd is very similar to 1, can't really tell much difference, perhaps slightly more tortuous
#and maybe slightly more often to have shorter SL but very slight
# will only keep the first two behaviours, the 3rd is too similar and <2% of obs


#behav 1 = exploratory
#behav 2 = rest


theta.estim.long.week <- expand_behavior(dat = tracks.seg.week, theta.estim = theta.estim.week, 
                                        obs = obs,
                                        nbehav = 2, 
                                        behav.names = c("Explore", "Rest"),
                                        behav.order = c(1,2))
##### Assign behavioural states to tracks

# Merge results with original data
tracks.out.hour <- assign_behavior(dat.orig = tracks.hour,
                             dat.seg.list = df_to_list(tracks.seg.hour, "id"),
                             theta.estim.long = theta.estim.long.hour,
                             behav.names = c("ARS", "Explore", "Rest"))

write.csv(tracks.out.hour, "tracksout_pesticidehour.csv")

tracks.out.day <- assign_behavior(dat.orig = tracks.day,
                                  dat.seg.list = df_to_list(tracks.seg.day, "id"),
                                  theta.estim.long = theta.estim.long.day,
                                  behav.names = c("StationKeeping", "BroadARS", "Rest"))

write.csv(tracks.out.day, "tracksout_pesticideDay.csv")

tracks.out.week <- assign_behavior(dat.orig = tracks.week,
                                  dat.seg.list = df_to_list(tracks.seg.week, "id"),
                                  theta.estim.long = theta.estim.long.week,
                                  behav.names = c("Explore", "Rest"))

write.csv(tracks.out.week, "tracksout_pesticideWeek.csv")


#### adding landcover for stats 


##to do stats with landcover need to add landcover data for true and random points
library(amt)
library(sf)
library(tidyverse)

landcover <-terra::rast("SpringLandcoverRaster.tif") #Rare landcover, made in overwintering

tracksoutdthour <- subset(tracks.out.hour, dt == 1800) #only tracks at sampling interval
tracksoutdthour.sp <-st_as_sf(tracksoutdthour, coords = c("x", "y"), 
                          crs = st_crs(landcover))
tracksoutdtday <- subset(tracks.out.day, dt == 43200) #only tracks at sampling interval
tracksoutdtday.sp <-st_as_sf(tracksoutdtday, coords = c("x", "y"), 
                              crs = st_crs(landcover))
tracksoutdtweek <- subset(tracks.out.week, dt == 302400) #only tracks at sampling interval
tracksoutdtweek.sp <-st_as_sf(tracksoutdtweek, coords = c("x", "y"), 
                             crs = st_crs(landcover))

tracksouthour.sp <-tracks.out.hour #all tracks, for when sampling interval doesnt matter
tracksouthour.sp <- st_as_sf(tracksouthour.sp, coords = c("x", "y"), 
                         crs = st_crs(landcover))
tracksoutday.sp <-tracks.out.day #all tracks, for when sampling interval doesnt matter
tracksoutday.sp <- st_as_sf(tracksoutday.sp, coords = c("x", "y"), 
                             crs = st_crs(landcover))
tracksoutweek.sp <-tracks.out.week #all tracks, for when sampling interval doesnt matter
tracksoutweek.sp <- st_as_sf(tracksoutweek.sp, coords = c("x", "y"), 
                             crs = st_crs(landcover))



tracksRandom <- st_sample(perim, size=9630)

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

randomCoords <- getRandomCoords(tracks.out.hour, 1) ## calculate all random points
randomCoordsdt <- getRandomCoords(tracksoutdthour, 1)

randomCoords <- getRandomCoords(tracks.out.day, 1) ## calculate all random points
randomCoordsdt <- getRandomCoords(tracksoutdtday, 1)

randomCoords <- getRandomCoords(tracks.out.week, 1) ## calculate all random points
randomCoordsdt <- getRandomCoords(tracksoutdtweek, 1)


getTrueCoords <- function(df){ ## revise the true steps to match the same as random
  revisedCoordsDF <- df %>% 
    group_by(id) %>%  ## select by bee identifier
    mutate(x2 = lead(x), y2 = lead(y)) %>%  ## create a new column for the end travel point
    filter(!is.na(x2)) %>%  ## drop all initial x-y that don't have an end point
    mutate(case = "true") %>%  ## add a column to identify true
    dplyr::select(case, x1 = x, y1 = y, x2, y2, sl = step, ta = angle)
  return(revisedCoordsDF)
}

trueCoordshour <- getTrueCoords(tracks.out.hour) ## revise data to get all true points
trueCoordsdthour <-getTrueCoords(tracksoutdthour)

trueCoordsday <- getTrueCoords(tracks.out.day) ## revise data to get all true points
trueCoordsdtday <- getTrueCoords(tracksoutdtday)

trueCoordsweek <- getTrueCoords(tracks.out.week) ## revise data to get all true points
trueCoordsdtweek <- getTrueCoords(tracksoutdtweek)

hourtracks.tf <- rbind(randomCoords,trueCoordshour ) 
hourtracksdt.tf <- rbind(randomCoordsdt, trueCoordsdthour)

daytracks.tf <- rbind(randomCoords,trueCoordsday ) 
daytracksdt.tf <- rbind(randomCoordsdt, trueCoordsdtday)

weektracks.tf <- rbind(randomCoords,trueCoordsweek ) 
weektracksdt.tf <- rbind(randomCoordsdt, trueCoordsdtweek)

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
hourtracks.tf <- hourtracks.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

hourtracks.tfs <- hourtracks.tf
hourtracks.tfs <- st_as_sf(hourtracks.tfs, coords = c("x2", "y2"),
                          crs=st_crs(landcover))

hourtracksdt.tf <- hourtracksdt.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

hourtracksdt.tfs <- hourtracksdt.tf
hourtracksdt.tfs <- st_as_sf(hourtracksdt.tfs, coords = c("x2", "y2"),
                           crs=st_crs(landcover))

daytracks.tf <- daytracks.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

daytracks.tfs <- daytracks.tf
daytracks.tfs <- st_as_sf(daytracks.tfs, coords = c("x2", "y2"),
                           crs=st_crs(landcover))

daytracksdt.tf <- daytracksdt.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

daytracksdt.tfs <- daytracksdt.tf
daytracksdt.tfs <- st_as_sf(daytracksdt.tfs, coords = c("x2", "y2"),
                             crs=st_crs(landcover))

weektracks.tf <- weektracks.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

weektracks.tfs <- weektracks.tf
weektracks.tfs <- st_as_sf(weektracks.tfs, coords = c("x2", "y2"),
                          crs=st_crs(landcover))

weektracksdt.tf <- weektracksdt.tf %>% 
  mutate(treatcode = substr(id, 0, 1))

weektracksdt.tfs <- weektracksdt.tf
weektracksdt.tfs <- st_as_sf(weektracksdt.tfs, coords = c("x2", "y2"),
                            crs=st_crs(landcover))

# Convert ' to a SpatVector object to work with terra
hourtracks.tfs <- vect(hourtracks.tfs)
landex <- terra::extract(landStack, hourtracks.tfs)

hourtracksdt.tfs <- vect(hourtracksdt.tfs)
landex <- terra::extract(landStack, hourtracksdt.tfs)

daytracks.tfs <- vect(daytracks.tfs)
landex <- terra::extract(landStack, daytracks.tfs)

daytracksdt.tfs <- vect(daytracksdt.tfs)
landex <- terra::extract(landStack, daytracksdt.tfs)

weektracks.tfs <- vect(weektracks.tfs)
landex <- terra::extract(landStack, weektracks.tfs)

weektracksdt.tfs <- vect(weektracksdt.tfs)
landex <- terra::extract(landStack, weektracksdt.tfs)

#add covariates from landStack to the dataframe for SSF
hourtracks.tfs$agriculture <- landex[,2]
hourtracks.tfs$developed <- landex[,3]
hourtracks.tfs$earlyFloral <- landex[,4]
hourtracks.tfs$forest <- landex[,5]
hourtracks.tfs$lateFloral <- landex[,6]
hourtracks.tfs$lowFloral <- landex[,7]
hourtracks.tfs$wetland<- landex[,8]
head(hourtracks.tfs)

hourtracks.tfs$presence <- ifelse(hourtracks.tfs$case == "random", 
                                  FALSE, TRUE)

hourtracksdt.tfs$agriculture <- landex[,2]
hourtracksdt.tfs$developed <- landex[,3]
hourtracksdt.tfs$earlyFloral <- landex[,4]
hourtracksdt.tfs$forest <- landex[,5]
hourtracksdt.tfs$lateFloral <- landex[,6]
hourtracksdt.tfs$lowFloral <- landex[,7]
hourtracksdt.tfs$wetland<- landex[,8]
head(hourtracksdt.tfs)

hourtracksdt.tfs$presence <- ifelse(hourtracksdt.tfs$case == "random", 
                                  FALSE, TRUE)

daytracks.tfs$agriculture <- landex[,2]
daytracks.tfs$developed <- landex[,3]
daytracks.tfs$earlyFloral <- landex[,4]
daytracks.tfs$forest <- landex[,5]
daytracks.tfs$lateFloral <- landex[,6]
daytracks.tfs$lowFloral <- landex[,7]
daytracks.tfs$wetland<- landex[,8]
head(daytracks.tfs)

daytracks.tfs$presence <- ifelse(daytracks.tfs$case == "random", 
                                  FALSE, TRUE)

daytracksdt.tfs$agriculture <- landex[,2]
daytracksdt.tfs$developed <- landex[,3]
daytracksdt.tfs$earlyFloral <- landex[,4]
daytracksdt.tfs$forest <- landex[,5]
daytracksdt.tfs$lateFloral <- landex[,6]
daytracksdt.tfs$lowFloral <- landex[,7]
daytracksdt.tfs$wetland<- landex[,8]
head(daytracksdt.tfs)

daytracksdt.tfs$presence <- ifelse(daytracksdt.tfs$case == "random", 
                                    FALSE, TRUE)

weektracks.tfs$agriculture <- landex[,2]
weektracks.tfs$developed <- landex[,3]
weektracks.tfs$earlyFloral <- landex[,4]
weektracks.tfs$forest <- landex[,5]
weektracks.tfs$lateFloral <- landex[,6]
weektracks.tfs$lowFloral <- landex[,7]
weektracks.tfs$wetland<- landex[,8]
head(weektracks.tfs)

weektracks.tfs$presence <- ifelse(weektracks.tfs$case == "random", 
                                 FALSE, TRUE)


weektracksdt.tfs$agriculture <- landex[,2]
weektracksdt.tfs$developed <- landex[,3]
weektracksdt.tfs$earlyFloral <- landex[,4]
weektracksdt.tfs$forest <- landex[,5]
weektracksdt.tfs$lateFloral <- landex[,6]
weektracksdt.tfs$lowFloral <- landex[,7]
weektracksdt.tfs$wetland<- landex[,8]
head(weektracksdt.tfs)

weektracksdt.tfs$presence <- ifelse(weektracksdt.tfs$case == "random", 
                                   FALSE, TRUE)

####################################
### Step-selection functions ########
#####################################

library(survival)
hourtracksdt.tfs.NONA <- data.frame(hourtracksdt.tfs)
hourtracksdt.tfs.NONA[is.na(hourtracksdt.tfs.NONA)] <- 0

SSF.hour.C <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral, 
                  family = "binomial", 
                data = subset(hourtracksdt.tfs.NONA, treatcode == "C"))
summary(SSF.hour.C)
anova(SSF.hour.C, test = "Chisq") #nothing sig

SSF.hour.I <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral, 
                     family = "binomial", 
                     data = subset(hourtracksdt.tfs.NONA, treatcode == "I"))
summary(SSF.hour.I)
anova(SSF.hour.I, test = "Chisq")

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.hour.I, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.hour.I)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")


ggplot(data = odds_ratio.df2, aes(x = landtype, y = oddsRatio)) +
  theme_classic() +
  geom_point(size = 4, colour = c("red","red","red", "white")) +
  geom_errorbar(aes( x = landtype , ymin = lowci, 
                     ymax = highci), width = 0, 
                colour = c("#235789", "#235789","#235789", "#235789")) +  
  theme(axis.title.x = element_text (size = 14), 
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(limits = c("earlyFloral", "forest", "lateFloral"),
                   labels = c("forest", "early floral", "late floral")) +
  scale_y_continuous(name = "odds ratio") +
  geom_hline(yintercept = 1, linetype = "solid", colour = "dark grey") +
  coord_flip() +
  ggtitle("imidacloprid hour") +
  theme(plot.title = element_text(size = 16))


SSF.hour.Y <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral, 
                  family = "binomial", 
                  data = subset(hourtracksdt.tfs.NONA, treatcode == "Y"))
summary(SSF.hour.Y)
anova(SSF.hour.Y, test = "Chisq")
 #wasn't enough points to really run an SSF. 

# day
daytracksdt.tfs.NONA <- data.frame(daytracksdt.tfs)
daytracksdt.tfs.NONA[is.na(daytracksdt.tfs.NONA)] <- 0

SSF.day.C <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral,
                  family = "binomial", 
                  data = subset(daytracksdt.tfs.NONA, treatcode == "C"))
summary(SSF.day.C)
anova(SSF.day.C, test = "Chisq") #forest, early floral and late floral sig

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.day.C, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.day.C)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

SSF.day.I <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral,
                 family = "binomial", 
                 data = subset(daytracksdt.tfs.NONA, treatcode == "I"))
summary(SSF.day.I)
anova(SSF.day.I, test = "Chisq") #forest, early floral and late floral sig

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.day.I, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.day.I)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

SSF.day.Y <- glm(presence ~ forest + earlyFloral + lateFloral + lowFloral,
                 family = "binomial", 
                 data = subset(daytracksdt.tfs.NONA, treatcode == "Y"))
summary(SSF.day.Y)
anova(SSF.day.Y, test = "Chisq") #forest, early floral and late floral sig

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.day.Y, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.day.Y)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

##  week
weektracksdt.tfs.NONA %>% 
  filter(presence == TRUE) %>% 
  group_by(treatcode) %>% 
  reframe(across(agriculture:wetland, ~sum(.x, na.rm = TRUE))) %>% 
  t() %>% 
  data.frame()

weektracksdt.tfs.NONA <- data.frame(weektracksdt.tfs)
weektracksdt.tfs.NONA[is.na(weektracksdt.tfs.NONA)] <- 0

SSF.week.C <- glm(presence ~ earlyFloral + lateFloral + lowFloral, 
                  family = "binomial", 
                  data = subset(weektracksdt.tfs.NONA, treatcode == "C"))
summary(SSF.week.C)
anova(SSF.week.C, test = "Chisq")

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.week.C, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.week.C)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

SSF.week.I <- glm(presence ~ earlyFloral + lateFloral + lowFloral, 
                  family = "binomial", 
                  data = subset(weektracksdt.tfs.NONA, treatcode == "I"))
summary(SSF.week.I)
anova(SSF.week.I, test = "Chisq")

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.week.I, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.week.I)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")


SSF.week.Y <- glm(presence ~ earlyFloral + lateFloral + lowFloral, 
                  family = "binomial", 
                  data = subset(weektracksdt.tfs.NONA, treatcode == "Y"))
summary(SSF.week.Y)
anova(SSF.week.Y, test = "Chisq")

# Get confidence intervals for coefficients
coeff_ci <- confint(SSF.week.Y, level = 0.95)

# Exponentiate the confidence intervals for odds ratios
odds_ratios_ci <- data.frame(exp(coeff_ci))
odds_ratios_ci$oddsRatio <- as.numeric(exp(coefficients(SSF.week.Y)))

names(odds_ratios_ci)[1] <-"lowci"
names(odds_ratios_ci)[2] <-"highci"

odds_ratios.df <- odds_ratios_ci[row.names(odds_ratios_ci) != "(Intercept)", ]
odds_ratio.df2 <- odds_ratios.df %>% 
  tibble::rownames_to_column("landtype")

#############################################
### Step Lengths and turn angles and NSD ###
#############################################

## start over again with tracksoutdthour

# needs treatcode, landcover, and daytime columns

#Separating dataset to make a daytime column

make_daytime <- function (data){
  #separate the data into day and time
  data <-tidyr::separate(data, date, into = c("day", "time"), sep = " ")
  #separate the time into hour, minute and seconds
  data <- tidyr::separate(data, time, sep = ":", into = c("hour", "minute", "sec"))
  # Convert hour to numeric to ensure comparisons work correctly
  data$hour <- as.numeric(data$hour)
  # Add a column to classify as day or night
  data$daytime <- ifelse(data$hour < 5 | data$hour > 20, "night", "day")
  return(data)
}


hourdt <- make_daytime(tracksoutdthour)
daydt <- make_daytime(tracksoutdtday)
weekdt <- make_daytime(tracksoutdtweek)

  
## adding landcover as categorical variable
library(sf)
library(dplyr)
library(terra)
library(pastecs)

landcover <-terra::rast("SpringLandcoverRaster.tif")

add_land_treatcode <- function(data, landcover, landClasses){
  
  # Step 1: Convert to sf object and reproject to match landcover CRS
  data_sf <- st_as_sf(data, coords = c("x", "y"), crs = st_crs(landcover))
  
  # Step 2: Convert to a SpatVector for terra compatibility
  data_vect <- vect(data_sf)
  
  # Step 3: Extract landcover values
  landex <- terra::extract(landcover, data_vect)
  data_vect$landtype <- landex[, 2] # Extract relevant column (e.g., landcover values)
  
  # Step 4: Convert back to a data frame and join with land cover classes
  data_df <- as.data.frame(data_vect)
  data_df <- data_df %>%
    left_join(landClasses, by = "landtype") %>%
    dplyr::select(-landtype) %>%
    rename(landtype = landTypeNew)
  
  # Step 5: Add a "treatcode" column from the first letter of "id"
  data_df <- data_df %>%
    mutate(treatcode = substr(id, 1, 1))
  
  # Step 6: Add behavioural columns and filter out rows with missing behaviour data
  data_df <- data_df %>%
    filter(!is.na(behav)) %>%
    mutate(
      ars = ifelse(behav == "ARS", 1, 0),
      Explore = ifelse(behav == "Explore", 1, 0),
      Rest = ifelse(behav == "Rest", 1, 0),
      StationKeep = ifelse(behav == "StationKeeping", 1, 0),
      BroadARS = ifelse(behav == "BroadARS", 1, 0)
    )

  return(data_df)
}

landClasses <- data.frame(
  landtype = 1:7,
  landTypeNew = c("agriculture", "developed", "earlyFloral", "forest", "lateFloral", "lowFloral", "wetland")
)

queentracks.hour2 <- add_land_treatcode(hourdt, landcover, landClasses)
queentracks.day2 <- add_land_treatcode(daydt, landcover, landClasses)
queentracks.week2 <- add_land_treatcode(weekdt, landcover, landClasses)

##  getting descriptive stats per treatcode

library(dplyr)

describe_steps <- function(data, treatcodes) {
  # Filter data by the specified treatcodes
  data_filtered <- data %>%
    filter(treatcode %in% treatcodes)  # Only include rows with specified treatcodes
  
  # Select only the columns: step, angle, and obs
  selected_columns <- data_filtered %>%
    select(step, angle, obs)  # Only select these specific columns
  
  # Calculate descriptive statistics by treatcode for selected columns
  stats_by_treatcode <- data_filtered %>%
    group_by(treatcode) %>%
    summarise(
      across(
        .cols = all_of(names(selected_columns)),  # Apply to the selected columns only
        list(
          mean = ~abs(mean(. , na.rm = TRUE)),
          se = ~sd(. , na.rm = TRUE) / sqrt(sum(!is.na(.))),
          min = ~min(. , na.rm = TRUE),
          max = ~max(. , na.rm = TRUE),
          median = ~median(. , na.rm = TRUE),
          IQR = ~IQR(. , na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  return(stats_by_treatcode)
}

treatcodes <- c("C", "I", "Y")
stats_hour <- describe_steps(queentracks.hour2, treatcodes)
stats_day <- describe_steps(queentracks.day2, treatcodes)
stats_week <- describe_steps(queentracks.week2, treatcodes)

write.csv(stats_hour, "stats_hour.csv")
write.csv(stats_day, "stats_day.csv")
write.csv(stats_week, "stats_week.csv")


#is there a difference between treatments for step length
library(car)
library(emmeans)

s1 <- lm(step~treatcode*daytime*landtype, data=queentracks.hour2)
car::Anova(s1, type =2) #significant for 3 way interaction
summary(s1)
pwr::pwr.f2.test(u = 2, v = NULL, f2 = 0.3507 , sig.level = 0.05, power = 0.8)
#there is enough power to test three way interactions
s1.a <- emmeans(s1, pairwise ~ treatcode * daytime * landtype, data = queentracks.hour2)
s1.a

s2 <- lm(step~treatcode*daytime*landtype, data=queentracks.day2)
car::Anova(s2, type =2) #significant for 3 way interaction
summary(s2)
pwr::pwr.f2.test(u = 2, v = NULL, f2 = 0.223 , sig.level = 0.05, power = 0.8)
#there is enough power to test three way interactions
s2.a <- emmeans(s2, pairwise ~ treatcode * daytime * landtype, data = queentracks.day2)
s2.a

s3 <- lm(step ~ treatcode *daytime*landtype, data=queentracks.week2)
car::Anova(s3, type =2) #significant for 3 way interaction
summary(s3)
pwr::pwr.f2.test(u = 2, v = NULL, f2 = 0.2603 , sig.level = 0.05, power = 0.8)
#there is enough power to test three way interactions
s3.a <- emmeans(s3, pairwise ~ treatcode * daytime * landtype, data = queentracks.week2)
s3.a

library(dplyr)
library(tidyr)
library(ggplot2)

queentracks.hour2 <- queentracks.hour2 %>%
  mutate(dataset = "hour")
queentracks.day2 <- queentracks.day2 %>%
  mutate(dataset = "day")
queentracks.week2 <- queentracks.week2 %>%
  mutate(dataset = "week")

# Combine datasets
queentracks_combined <- bind_rows(queentracks.hour2, queentracks.day2, queentracks.week2) %>%
  mutate(dataset = factor(dataset, levels = c("hour", "day", "week")),  # Set row order
         daytime = recode(daytime, "day" = "day time", "night" = "night time"))  # Rename columns

# Summarize and complete data
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))

stepplot <- queentracks_combined %>%
  group_by(dataset, treatcode, daytime, landtype) %>%
  summarize(avgStep = mean(step, na.rm = TRUE),
            errorstep = se(step),
            .groups = "drop_last") %>%
  ungroup() %>%
  tidyr::complete(dataset, treatcode, daytime, landtype, fill = list(avgStep = 0, errorstep = 0))

# Plot
ggplot(data = stepplot, aes(x = landtype, y = avgStep + 1, fill = treatcode)) +
  theme_classic() +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    breaks = c("C", "I", "Y"), 
                    labels = c("control", "imidacloprid", "cyantraniliprole")) +
  geom_errorbar(aes(x = landtype, ymin = avgStep - errorstep, ymax = avgStep + errorstep), 
                width = 0, position = position_dodge(width = 0.8)) +
  facet_grid(dataset ~ daytime) +  # Facet by dataset (rows) and daytime (columns)
  theme(strip.text.x = element_text(size = 13),  # Customize column headers
        strip.text.y = element_text(size = 13)) +  # Customize row headers
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  scale_y_continuous(name = "average step length (m)") +
  scale_x_discrete(limits = c("earlyFloral", "lateFloral", "lowFloral"),
                   labels = c("early floral", "late floral", "low floral")) +
  theme(legend.title = element_blank())


### Turning angle stats
library(car)
library(emmeans)

t1 <- lm(abs(angle) ~ treatcode*daytime*landtype, data=queentracks.hour2)
car::Anova(t1, type =2) #significant for 3 way interaction
summary(t1)
pwr::pwr.f2.test(u = 2, v = NULL, f2 = 0.3723 , sig.level = 0.05, power = 0.8)
#there is enough power to test three way interactions n = 589
t1.a <- emmeans(t1, pairwise ~ treatcode * daytime * landtype, data = queentracks.hour2)
t1.a

t2 <- lm(abs(angle) ~ treatcode * daytime*landtype, data=queentracks.day2)
car::Anova(t2, type =2) #not sig for 3 way interaction
summary(t2)
pwr::pwr.f2.test(u = 1, v = NULL, f2 = 0.3354 , sig.level = 0.05, power = 0.8)
#there is enough power to test two way interactions
t2.a <- emmeans(t2, pairwise ~ treatcode * daytime, data = queentracks.day2)
t2.a
t2.b <- emmeans(t2, pairwise ~ treatcode * landtype, data = queentracks.day2)
t2.b

t3 <- lm(abs(angle) ~ treatcode * daytime *landtype, data=queentracks.week2)
car::Anova(t3, type =2) #significant for 3 way interaction
summary(t3)
pwr::pwr.f2.test(u = 2, v = NULL, f2 = 0.3489 , sig.level = 0.05, power = 0.8)
#there is enough power to test three way interactions
t3.a <- emmeans(t3, pairwise ~ treatcode * daytime * landtype, data = queentracks.week2)
t3.a

# Turning Angles Plot
library(dplyr)
library(tidyr)
library(ggplot2)

# Summarize and complete data
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))

angleplot <- queentracks_combined %>%
  group_by(dataset, treatcode, daytime, landtype) %>%
  summarize(avgAngle = mean(abs(angle), na.rm = TRUE),
            errorangle = se(abs(angle)),
            .groups = "drop_last") %>%
  ungroup() %>%
  tidyr::complete(dataset, treatcode, daytime, landtype, fill = list(avgAngle = 0, errorangle = 0))

ggplot(data = angleplot, aes(x = landtype, y = avgAngle, fill = treatcode)) +
  theme_classic() +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    breaks = c("C", "I", "Y"), 
                    labels = c("control", "imidacloprid", "cyantraniliprole")) +
  geom_errorbar(aes(x = landtype, ymin = avgAngle - errorangle, ymax = avgAngle + errorangle), 
                width = 0, position = position_dodge(width = 0.8)) +
  facet_grid(dataset ~ daytime) +  # Facet by dataset (rows) and daytime (columns)
  theme(strip.text.x = element_text(size = 13),  # Customize column headers
        strip.text.y = element_text(size = 13)) +  # Customize row headers
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  scale_y_continuous(name = "average turning angle (rad)") +
  scale_x_discrete(limits = c("earlyFloral", "lateFloral", "lowFloral"),
                   labels = c("early floral", "late floral", "low floral")) +
  theme(legend.title = element_blank())


### total flight length

tf.dat <- queentracks_combined %>% 
  group_by(id, treatcode, dataset) %>%
  summarise(total_steps = sum(step, na.rm = TRUE)) %>%
  ungroup()

#hour
f1 <- lm(total_steps ~ treatcode, data = subset(tf.dat, dataset =="hour"))
summary(f1)
anova(f1) # no difference

f2 <- lm(total_steps ~ treatcode, data = subset(tf.dat, dataset =="day"))
summary(f2)
anova(f2) # no difference

f3 <- lm(total_steps ~ treatcode, data = subset(tf.dat, dataset =="week"))
summary(f3)
anova(f3) #no difference 


####   behaviour as a response variable  ####
library(car)
library(emmeans)

##glm power analysis function
library(pwr)

calculate_model_metrics <- function(model, sig.level = 0.05, power = NULL) {
  # Extract deviance
  residual_deviance <- model$deviance
  null_deviance <- model$null.deviance
  
  # Calculate R^2
  r2 <- 1 - (residual_deviance / null_deviance)
  
  # Calculate effect size f^2
  f2 <- r2 / (1 - r2)
  
  # Degrees of freedom
  df_null <- model$df.null
  df_residual <- model$df.residual
  u <- df_null - df_residual
  v <- df_residual
  
  # Power analysis
  power_analysis <- pwr.f2.test(u = u, v = v, f2 = f2, sig.level = sig.level, power = power)
  
  # Return results as a list
  return(list(
    R2 = r2,
    f2 = f2,
    df_null = df_null,
    df_residual = df_residual,
    u = u,
    v = v,
    power_analysis = power_analysis
  ))
}

# hour
# behaviours hour = ARS, explore, rest
b1.hour <- glm(ARS ~ treatcode * daytime * landtype, family="binomial", 
               data=subset(queentracks_combined, dataset == "hour")) 
summary(b1.hour)
results <- calculate_model_metrics(b1.hour)
results #have just enough power for interactions n = 589

anova(b1.hour, test = "Chisq")

b1.hour.a <- emmeans::emmeans (b1.hour, 
                               pairwise ~ treatcode * landtype, 
                               data = subset(queentracks_combined, 
                                             dataset == "hour"))

#note for explore the model would not coverge with daytime
b2.hour <- glm(Explore ~ treatcode  * landtype *daytime, family="binomial", 
               data=subset(queentracks_combined, dataset == "hour"))
summary(b2.hour)
results <- calculate_model_metrics(b2.hour)
results #have just enough power for interactions
anova(b2.hour, test="Chisq")

b2.hour.a <- emmeans::emmeans (b2.hour, 
                               pairwise ~ treatcode * landtype, 
                               data = subset(queentracks_combined, 
                                             dataset == "hour"))

b3.hour <- glm(Rest ~ treatcode  * daytime * landtype, family="binomial", 
               data=subset(queentracks_combined, dataset == "hour"))
summary(b3.hour)
results <- calculate_model_metrics(b3.hour)
results #have just enough power for interactions

anova(b3.hour, test = "Chisq")

b3.hour.a <- emmeans::emmeans (b3.hour, 
                               pairwise ~ treatcode , 
                               data = subset(queentracks_combined, 
                                             dataset == "hour"),
                               adjust = "none")

##  day
# response variables:"StationKeeping", "BroadARS", "Rest"

# station keeping
b1.day <- glm(StationKeeping ~ treatcode * daytime * landtype, family="binomial", 
               data=subset(queentracks_combined, dataset == "day")) 
str(subset(queentracks_combined, dataset == "day")) # 849 obs
summary(b1.day)
results <- calculate_model_metrics(b1.day)
results #there is enough power
anova(b1.day, test ="Chisq")

b1.day.a <- emmeans::emmeans (b1.day, 
                               pairwise ~ treatcode * landtype, 
                               data = subset(queentracks_combined, 
                                             dataset == "day"))
b1.day.a


#BroadARS
b2.day <- glm(BroadARS ~ treatcode * daytime * landtype, family="binomial", 
              data=subset(queentracks_combined, dataset == "day")) 
summary(b2.day)
results <- calculate_model_metrics(b2.day)
results #there is enough power
anova(b2.day, test = "Chisq")

b2.day.a <- emmeans::emmeans (b2.day, 
                               pairwise ~ treatcode * daytime * landtype, 
                               data = subset(queentracks_combined, 
                                             dataset == "day"))
b2.day.a


#Rest

b3.day <- glm(Rest ~ treatcode * daytime * landtype, family="binomial", 
              data=subset(queentracks_combined, dataset == "day")) 
summary(b3.day)
results <- calculate_model_metrics(b3.day)
results # there is enough power

anova(b3.day, test = "Chisq")

b3.day.a <- emmeans::emmeans (b3.day, 
                               pairwise ~ treatcode, 
                               data = subset(queentracks_combined, 
                                             dataset == "day"))
b3.day.a

##  week

# Explore
b1.week <- glm(Explore ~ treatcode * daytime * landtype, family="binomial", 
              data=subset(queentracks_combined, dataset == "week")) 
summary(b1.week)
str(subset(queentracks_combined, dataset == "week")) #n = 934
results <- calculate_model_metrics(b1.week)
results #there is enough power
anova(b1.week, test = "Chisq")

b1.week.a <- emmeans::emmeans (b1.week, 
                               pairwise ~ treatcode* daytime * landtype, 
                               data = subset(queentracks_combined, 
                                             dataset == "week"))
b1.week.a



# behaviour plots

library(ggplot2)
library(dplyr)
library(tidyr)

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

#hour
#ARS and Explore, treatcode x landcover is sig. for Rest just treatcode sig for figs

behav.hour <- subset(queentracks_combined, dataset == "hour", 
                     select = c("treatcode", "ARS", "Explore", "Rest", "landtype"))


behav.hourPlot <- behav.hour %>% 
  group_by(treatcode, landtype) %>% 
  summarize(
    avgARS = mean(ARS, na.rm = TRUE), 
    errorARS = se(ARS),
    avgExplore = mean(Explore, na.rm = TRUE),
    errorExplore = se(Explore), 
    avgRest = mean(Rest, na.rm = TRUE),
    errorRest = se(Rest),
    .groups = "drop_last"  # Specify .groups only once here
  ) %>%
  ungroup() %>%
  tidyr::complete(
    treatcode, landtype, 
    fill = list(
      avgARS = 0, errorARS = 0,
      avgExplore = 0, errorExplore = 0,
      avgRest = 0, errorRest = 0
    )
  )

behav.hourPlot2 <- behav.hourPlot %>%
  pivot_longer(
    cols = starts_with("avg"),  # Columns containing response variables
    names_to = "response_var",  # New column for response variable names
    values_to = "response_value"  # New column for response values
  ) %>%
  mutate(
    error = case_when(
      response_var == "avgARS" ~ errorARS,
      response_var == "avgExplore" ~ errorExplore,
      response_var == "avgRest" ~ errorRest
    ),
    response_var = factor(response_var, levels = c("avgARS", "avgExplore", "avgRest"), 
                          labels = c("ARS", "Explore", "Rest"))
  )


ggplot(behav.hourPlot2, aes(x = landtype, y = response_value, fill = treatcode)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = response_value - error, ymax = response_value + error),
                width = 0, position = position_dodge(width = 0.8)) +
  facet_wrap(~response_var, scales = "free_y") +  # Create a facet for each response variable
  theme_classic() +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    breaks = c("C", "I", "Y"), 
                    labels = c("control", "imidacloprid", "cyantraniliprole")) +
  scale_x_discrete(limits = c("earlyFloral", "lateFloral", "lowFloral"),
                   labels = c("early floral", "late floral", "low floral"),
                   name = "") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)  # Customize facet labels
  ) +
  labs(y = "Mean proportion of flight") ## export as 12 x 6 

#day 
#station keeping treatcode x landtype is sig
#BroadARS three way interaction is significant
#Rest only treatcode is significant


behav.day <- subset(queentracks_combined, dataset == "day", 
                     select = c("treatcode", "StationKeeping", "BroadARS", "Rest", 
                                "landtype", "daytime"))

behav.dayPlot <- behav.day %>% 
  group_by(treatcode, landtype, daytime) %>% 
  summarize(
    avgSK = mean(StationKeeping, na.rm = TRUE), 
    errorSK = se(StationKeeping),
    avgBroadARS = mean(BroadARS, na.rm = TRUE),
    errorBroadARS = se(BroadARS), 
    avgRest = mean(Rest, na.rm = TRUE),
    errorRest = se(Rest),
    .groups = "drop_last"  # Specify .groups only once here
  ) %>%
  ungroup() %>%
  tidyr::complete(
    treatcode, 
    fill = list(
      avgSK = 0, errorSK = 0,
      avgBroadARS = 0, errorBroadARS = 0,
      avgRest = 0, errorRest = 0
    )
  )

behav.dayPlot2 <- behav.dayPlot %>%
  pivot_longer(
    cols = starts_with("avg"),  # Columns containing response variables
    names_to = "response_var",  # New column for response variable names
    values_to = "response_value"  # New column for response values
  ) %>%
  mutate(
    error = case_when(
      response_var == "avgSK" ~ errorSK,
      response_var == "avgBroadARS" ~ errorBroadARS,
      response_var == "avgRest" ~ errorRest
    ),
    response_var = factor(response_var, levels = c("avgSK", "avgBroadARS", "avgRest"), 
                          labels = c("StationKeeping", "BroadARS", "Rest"))
  ) %>%
  complete(landtype, treatcode, response_var, daytime, fill = list(response_value = 0))  # Fill missing combinations with 0


ggplot(behav.dayPlot2, aes(x = landtype, y = response_value, fill = treatcode)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  # Adjust width for separation
  geom_errorbar(aes(ymin = response_value - error, ymax = response_value + error),
                width = 0, position = position_dodge(width = 0.9)) +  # Same dodge width for error bars
  facet_grid(daytime ~ response_var, scales = "free_y") +  
  theme_classic() +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    breaks = c("C", "I", "Y"), 
                    labels = c("control", "imidacloprid", "cyantraniliprole")) +
  scale_x_discrete(limits = c("earlyFloral", "lateFloral", "lowFloral"),
                   labels = c("early floral", "late floral", "low floral"),
                   name = "") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "right",  # Position the legend to the right
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)  # Customize facet labels
  ) +
  labs(y = "Mean proportion of flight")


# week
# Explore, three way interaction is significant

behav.week <- subset(queentracks_combined, dataset == "week", 
                    select = c("treatcode", "Explore", "Rest", 
                               "landtype", "daytime"))

behav.weekPlot <- behav.week %>% 
  group_by(treatcode, landtype, daytime) %>% 
  summarize(
    avgExplore = mean(Explore, na.rm = TRUE), 
    errorExplore = se(Explore),
    avgRest = mean(Rest, na.rm = TRUE),
    errorRest = se(Rest), 
    .groups = "drop_last"  # Specify .groups only once here
  ) %>%
  ungroup() %>%
  tidyr::complete(
    treatcode, 
    fill = list(
      avgExplore = 0, errorExplore = 0,
      avgRest = 0, errorRest = 0
    )
  )

behav.weekPlot2 <- behav.weekPlot %>%
  pivot_longer(
    cols = starts_with("avg"),  # Columns containing response variables
    names_to = "response_var",  # New column for response variable names
    values_to = "response_value"  # New column for response values
  ) %>%
  mutate(
    error = case_when(
      response_var == "avgExplore" ~ errorExplore,
      response_var == "avgRest" ~ errorRest
    ),
    response_var = factor(response_var, levels = c("avgExplore", "avgRest"), 
                          labels = c("Explore", "Rest"))
  ) %>%
  complete(landtype, treatcode, response_var, daytime, fill = list(response_value = 0))  # Fill missing combinations with 0


ggplot(behav.weekPlot2, aes(x = landtype, y = response_value, fill = treatcode)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  # Adjust width for separation
  geom_errorbar(aes(ymin = response_value - error, ymax = response_value + error),
                width = 0, position = position_dodge(width = 0.9)) +  # Same dodge width for error bars
  facet_grid(daytime ~ response_var, scales = "free_y") +  
  theme_classic() +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    breaks = c("C", "I", "Y"), 
                    labels = c("control", "imidacloprid", "cyantraniliprole")) +
  scale_x_discrete(limits = c("earlyFloral", "lateFloral", "lowFloral"),
                   labels = c("early floral", "late floral", "low floral"),
                   name = "") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "right",  # Position the legend to the right
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)  # Customize facet labels
  ) +
  labs(y = "Mean proportion of flight")

#################################################
### Kernel Density  ##########################
############################################

library(MASS)
library(terra)
library (ggplot2)
library(dplyr)
library(tidyr)

landcover <-terra::rast("SpringLandcoverRaster.tif") #Rare landcover, made in overwintering
tracks.out.day

##  separating the dataset to make KD rasters 

process_pointsData <- function(data) {
  # add column treatcode
  data$treatcode <- substr(data$id, 1, 1)
  
  # Extract hour and classify into 'day' or 'night'
  data$daytime <- ifelse(as.numeric(format(data$date, "%H")) > 20 | as.numeric(format(data$date, "%H")) < 5, "night", "day")
  
  # Create subsets for each combination of treatcode and daytime
  subsets <- list(
    controlDay   = subset(data, treatcode == "C" & daytime == "day"),
    controlNight = subset(data, treatcode == "C" & daytime == "night"),
    imidDay      = subset(data, treatcode == "I" & daytime == "day"),
    imidNight    = subset(data, treatcode == "I" & daytime == "night"),
    cyanDay      = subset(data, treatcode == "Y" & daytime == "day"),
    cyanNight    = subset(data, treatcode == "Y" & daytime == "night")
  )
  
  return(subsets)
}

# Apply the function to each dataset
kdhour_subsets <- process_pointsData(tracks.out.hour)
kdday_subsets <- process_pointsData(tracks.out.day)
kdweek_subsets <- process_pointsData(tracks.out.week)

# Access individual datasets
#hour
controlDay_hour   <- kdhour_subsets$controlDay
controlNight_hour <- kdhour_subsets$controlNight
imidDay_hour      <- kdhour_subsets$imidDay
imidNight_hour    <- kdhour_subsets$imidNight
cyanDay_hour      <- kdhour_subsets$cyanDay
cyanNight_hour    <- kdhour_subsets$cyanNight
#day
controlDay_day   <- kdday_subsets$controlDay
controlNight_day <- kdday_subsets$controlNight
imidDay_day      <- kdday_subsets$imidDay
imidNight_day    <- kdday_subsets$imidNight
cyanDay_day      <- kdday_subsets$cyanDay
cyanNight_day    <- kdday_subsets$cyanNight
#week
controlDay_week   <- kdweek_subsets$controlDay
controlNight_week <- kdweek_subsets$controlNight
imidDay_week      <- kdweek_subsets$imidDay
imidNight_week    <- kdweek_subsets$imidNight
cyanDay_week      <- kdweek_subsets$cyanDay
cyanNight_week    <- kdweek_subsets$cyanNight

#calculating kernel density values
##bandwidth formula h = (1.06 * σ) / (n^0.2)

listToRast <- function(x, y, z) { #need this because terra doesn't do it
  resx <- ( x[length(x)] - x[1] ) / (length(x)-1)
  resy <- ( y[length(y)] - y[1] ) / (length(y)-1)
  xmn <- min(x) - 0.5 * resx
  xmx <- max(x) + 0.5 * resx
  ymn <- min(y) - 0.5 * resy
  ymx <- max(y) + 0.5 * resy
  z <- t(z)
  z <- z[nrow(z):1, ]
  r1 <- terra::rast(terra::ext(xmn, xmx, ymn, ymx), resolution=c(resx, resy))
  r1[] <- as.numeric(z)
  return(r1)
}
randomAdjustment <- function(nObs, adjRange){
  return(runif(nObs, -0.5, 0.5)/adjRange)
}

datasetList <- list(controlDay_hour, controlNight_hour, 
                    imidDay_hour, imidNight_hour, 
                    cyanDay_hour, cyanNight_hour,
                    controlDay_day, controlNight_day, 
                    imidDay_day, imidNight_day, 
                    cyanDay_day, cyanNight_day,
                    controlDay_week, controlNight_week, 
                    imidDay_week, imidNight_week, 
                    cyanDay_week, cyanNight_week)
names(datasetList) <- c("controlDay_hour", "controlNight_hour", 
                        "imidDay_hour", "imidNight_hour", 
                        "cyanDay_hour", "cyanNight_hour",
                        "controlDay_day", "controlNight_day", 
                        "imidDay_day", "imidNight_day", 
                        "cyanDay_day", "cyanNight_day",
                        "controlDay_week", "controlNight_week", 
                        "imidDay_week", "imidNight_week", 
                        "cyanDay_week", "cyanNight_week")


makeRaster <- function(x)
{
  bee.kern <- x
  bee.kern2 <-bee.kern[, c("x", "y")]
  sdX <-sd(bee.kern2$x)
  nX <- length(bee.kern2$x)
  hx = (1.06 * (sdX/(nX^0.2))) #bandwidth function for x
  sdY <-sd(bee.kern2$y)
  nY <- length(bee.kern2$y)
  hy = (1.06 *sdY/nY^0.2) #bandwidth function for y
  bee.kd<-kde2d(bee.kern2$x, bee.kern2$y, h=c(hx, hy)) # getting kernel density estimates
  beeRaster <- listToRast(bee.kd$x,  bee.kd$y,  bee.kd$z) # using raster function from above
  # filename <- paste0(names(x), ".tif")
  # writeRaster(beeRaster, filename, overwrite = TRUE)
  return(beeRaster)
}

controlDay_hour <- makeRaster(datasetList[[1]])
controlNight_hour <- makeRaster(datasetList[2])
imidDay_hour <- makeRaster(datasetList[[3]])
imidNight_hour <- makeRaster(datasetList[[4]])
cyanDay_hour <- makeRaster(datasetList[[5]])
cyanNight_hour <- makeRaster(datasetList[[6]])

controlDay_day <- makeRaster(datasetList[[7]])
controlNight_day <- makeRaster(datasetList[[8]])
imidDay_day <- makeRaster(datasetList[[9]])
imidNight_day <- makeRaster(datasetList[[10]])
cyanDay_day <- makeRaster(datasetList[[11]])
cyanNight_day <- makeRaster(datasetList[[12]])

controlDay_week <- makeRaster(datasetList[[13]])
controlNight_week <- makeRaster(datasetList[[14]])
imidDay_week <- makeRaster(datasetList[[15]])
imidNight_week <- makeRaster(datasetList[[16]])
cyanDay_week <- makeRaster(datasetList[[17]])
cyanNight_week <- makeRaster(datasetList[[18]])

#for some reason controlNight_hour would not work with the function so running separately
bee.kern <- controlNight_hour
bee.kern2 <-bee.kern[, c("x", "y")]
sdX <-sd(bee.kern2$x)
nX <- length(bee.kern2$x)
hx = (1.06 * (sdX/(nX^0.2))) #bandwidth function for x
sdY <-sd(bee.kern2$y)
nY <- length(bee.kern2$y)
hy = (1.06 *sdY/nY^0.2) #bandwidth function for y
bee.kd<-kde2d(bee.kern2$x, bee.kern2$y, h=c(hx, hy)) # getting kernel density estimates
controlNight_hour <- listToRast(bee.kd$x,  bee.kd$y,  bee.kd$z)


# plotting kernel density values
#need to change crs for plotting with leaflet
library(terra)
library(sf)
library(leaflet)
library(colorspace)

#setting up hour
crs(controlDay_hour) <-"ESRI:102002"
CDayHour <- terra::project(controlDay_hour, "epsg:4326")

crs(controlNight_hour) <-"ESRI:102002"
CNightHour <- terra::project(controlNight_hour, "epsg:4326")

crs(cyanDay_hour) <-"ESRI:102002"
YDayHour <- terra::project(cyanDay_hour, "epsg:4326")

crs(cyanNight_hour) <-"ESRI:102002"
YNightHour <- terra::project(cyanNight_hour, "epsg:4326")

crs(imidDay_hour) <-"ESRI:102002"
IDayHour <- terra::project(imidDay_hour, "epsg:4326")

crs(imidNight_hour) <-"ESRI:102002"
INightHour <- terra::project(imidNight_hour, "epsg:4326")

#setting up day
crs(controlDay_day) <-"ESRI:102002"
CDayDay <- terra::project(controlDay_day, "epsg:4326")

crs(controlNight_day) <-"ESRI:102002"
CNightDay <- terra::project(controlNight_day, "epsg:4326")

crs(cyanDay_day) <-"ESRI:102002"
YDayDay <- terra::project(cyanDay_day, "epsg:4326")

crs(cyanNight_day) <-"ESRI:102002"
YNightDay <- terra::project(cyanNight_day, "epsg:4326")

crs(imidDay_day) <-"ESRI:102002"
IDayDay <- terra::project(imidDay_day, "epsg:4326")

crs(imidNight_day) <-"ESRI:102002"
INightDay <- terra::project(imidNight_day, "epsg:4326")

#setting up week
crs(controlDay_week) <-"ESRI:102002"
CDayWeek <- terra::project(controlDay_week, "epsg:4326")

crs(controlNight_week) <-"ESRI:102002"
CNightWeek <- terra::project(controlNight_week, "epsg:4326")

crs(cyanDay_week) <-"ESRI:102002"
YDayWeek <- terra::project(cyanDay_week, "epsg:4326")

crs(cyanNight_week) <-"ESRI:102002"
YNightWeek <- terra::project(cyanNight_week, "epsg:4326")

crs(imidDay_week) <-"ESRI:102002"
IDayWeek <- terra::project(imidDay_week, "epsg:4326")

crs(imidNight_week) <-"ESRI:102002"
INightWeek <- terra::project(imidNight_week, "epsg:4326")

#set up colour palette for control
cPal<-sequential_hcl(6, palette = "Purples 3", rev=TRUE)

max_value <- max(values(CDayHour), na.rm = T)
cval = c(0,seq(0, max_value , by = max_value/5))

cpal = c("#FFFFFF00", cPal)

#set up colour palette for cyantraniliprole
yPal<-sequential_hcl(6, palette = "Teal", rev=TRUE)

#set up colour palette for imidacloprid
iPal<-sequential_hcl(6, palette = "OrYel", rev=TRUE)

# control day hour
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CDayHour, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# control night hour
#set sequential values using max kd
max_value <- max(values(CNightHour), na.rm = T)
cval = c(0,seq(0, max_value, by = max_value/5))

cpal = c("#FFFFFF00", cPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CNightHour, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# cyan day hour
max_value <- max(values(YDayHour), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YDayHour, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")
# cyan night hour
max_value <- max(values(YNightHour), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YNightHour, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")
#imid day hour
max_value <- max(values(IDayHour), na.rm = T)
Ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(IDayHour, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")
#imid night hour
max_value <- max(values(INightHour), na.rm = T)
Ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(INightHour, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")

# control day day
max_value <- max(values(CDayDay), na.rm = T)
cval = c(0,seq(0, max_value, by = max_value/5))

leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CDayDay, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# control night day
#set sequential values using max kd
max_value <- max(values(CNightDay), na.rm = T)
cval = c(0,seq(0, max_value, by = max_value/5))

cpal = c("#FFFFFF00", cPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CNightDay, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# cyan day day
#set sequential values using max kd
max_value <- max(values(YDayDay), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YDayDay, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")

# cyan night day
max_value <- max(values(YNightDay), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YNightDay, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")

#imid day day
max_value <- max(values(IDayDay), na.rm = T)
ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(IDayDay, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")

#imid night day
max_value <- max(values(INightDay), na.rm = T)
ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(INightDay, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")

# control day week
max_value <- max(values(CDayWeek), na.rm = T)
cval = c(0,seq(0, max_value, by = max_value/5))

cpal = c("#FFFFFF00", cPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CDayWeek, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# control night week
#set sequential values using max kd
max_value <- max(values(CNightWeek), na.rm = T)
cval = c(0,seq(0, max_value, by = max_value/5))

cpal = c("#FFFFFF00", cPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(CNightWeek, colors = cpal, opacity = 0.60) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)],
            title = "kernel density")

# cyan day week
max_value <- max(values(YDayWeek), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YDayWeek, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")

# cyan night week
max_value <- max(values(YNightWeek), na.rm = T)
yval = c(0,seq(0, max_value, by = max_value/5))

ypal = c("#FFFFFF00", yPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(YNightWeek, colors = ypal, opacity = 0.65) %>% 
  addLegend(colors=ypal[-c(1,2)], labels = yval[-c(1,2)],
            title = "kernel density")

#imid day week
max_value <- max(values(IDayWeek), na.rm = T)
ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(IDayWeek, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")

#imid night week
max_value <- max(values(INightWeek), na.rm = T)
ival = c(0,seq(0, max_value, by = max_value/5))

ipal = c("#FFFFFF00", iPal)
leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(INightWeek, colors = ipal, opacity = 0.65) %>% 
  addLegend(colors=ipal[-c(1,2)], labels = ival[-c(1,2)],
            title = "kernel density")

# combining datasets of kernel density for running stats

# need to extract landcover for stats
landcover2 <- aggregate (landcover, fact = 6)
landcoverPoly <- terra::as.points(landcover2)

rasters <- list(CDayHour, CNightHour, YDayHour, YNightHour, IDayHour, INightHour,
                CDayDay, CNightDay, YDayDay, YNightDay, IDayDay, INightDay,
                CDayWeek, CNightWeek, YDayWeek, YNightWeek, IDayWeek, INightWeek)
names(rasters) <- c("CDH", "CNH", "YDH", "YNH", "IDH", "INH", 
                    "CDD", "CND", "YDD","YND", "IDD", "IND",
                    "CDW", "CNW", "YDW", "YNW", "IDW", "INW")

# Initialize an empty list to store the extracted values
values_list <- list()

# Loop through each raster and extract combined values
for (i in seq_along(rasters)) {
  # Extract raster values at points defined by landcoverPoly
  
  landcoverPolyWGS84 <- terra::project(landcoverPoly, rasters[[i]])
  extracted_raster <- terra::extract(rasters[[i]], landcoverPolyWGS84)
  extracted_raster <- extracted_raster[, -1] 
  
  # Create a data frame for landcoverPoly values
  landcover_df <- as.data.frame(landcoverPoly)
  
  # Combine extracted raster values with landcoverPoly values
  combined_values <- cbind(landcover_df, extracted_raster)
  
  # Rename the column according to raster name
  colnames(combined_values)[ncol(combined_values)] <- names(rasters)[i]
  
  # Store combined values in values_list
  values_list[[i]] <- combined_values
}
# Combine extracted values into a data frame
combined.df <- as.data.frame(values_list)

# Initialize an empty data frame with columns
final_df <- data.frame(
)

# Loop through each column in combined.df
for (col in colnames(combined.df[seq(2,20,by=2)])) {
  # Extract the values from the current column
  values <- combined.df[[col]]
  
  # Split the three-letter code into components
  treatment_code <- substr(col, 1, 1)
  daytime_code <- substr(col, 2, 2)
  type_code <- substr(col, 3, 3)
  
  # Map the codes to their respective categories
  treatment <- ifelse(treatment_code == "C", "control",
                      ifelse(treatment_code == "Y", "cyantraniliprole", "imidacloprid"))
  daytime <- ifelse(daytime_code == "D", "day", "night")
  type <- ifelse(type_code == "H", "hour", 
                 ifelse(type_code == "D" ,"day", "week"))
  
  # Create a data frame with the new columns and values
  temp_df <- data.frame(
    treatment = treatment,
    daytime = daytime,
    type = type,
    value = values,
    landcover = combined.df$SpringLandcoverRaster,
    stringsAsFactors = FALSE
  )
  
  # Append temp_df to final_df
  final_df <- rbind(final_df, temp_df)
}

##  kernel density statistics

#separating into three datasets by time period

hour.kd <- subset(final_df, type == "hour")
day.kd <- subset(final_df, type == "day")
week.kd <-subset(final_df, type == "week")

# three way interaction with 

kd1 <- glm(value ~ treatment * daytime *as.factor(landcover), 
           family = "binomial", data= hour.kd)
summary(kd1)
anova(kd1, test = "Chisq")

kd2 <- glm(value ~ treatment * daytime * as.factor(landcover), 
           family = "binomial", data= day.kd)

summary(kd2)
anova(kd2, test = "Chisq")

kd3 <- glm(value ~ treatment * as.factor(landcover), 
           family = "binomial", data= week.kd)

summary(kd3)
anova(kd3, test = "Chisq")

################################################
####  Minimum Convex Polygon ####
###############################################

library(adehabitatHR) # this package still only runs of sp, will use that for now
library(dplyr)
library(sf)
library(terra)

# finding queens with at least 5 detections for hour
hour5 <- tracks.out.hour %>% 
  group_by(id) %>%
  mutate(nHits = length(id)) %>% 
  filter(nHits > 5) #mcp function only works with at least 5 relocations 
hour5 <- hour5[, c("id", "x", "y")] #mcp needs only 3 columns
sp::coordinates(hour5) <- ~x+y #specifying coordinates
sp::proj4string(hour5) <-"+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45" #Lambert Conformal Conic

#finding queens with at least 5 detections for day
day5 <- tracks.out.day %>% 
  group_by(id) %>%
  mutate(nHits = length(id)) %>% 
  filter(nHits > 5) #mcp function only works with at least 5 relocations 
day5 <- day5[, c("id", "x", "y")] #mcp needs only 3 columns
sp::coordinates(day5) <- ~x+y #specifying coordinates
sp::proj4string(day5) <-"+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45" #Lambert Conformal Conic

#finding queens with at least 5 detections for week
week5 <- tracks.out.week %>% 
  group_by(id) %>%
  mutate(nHits = length(id)) %>% 
  filter(nHits > 5) #mcp function only works with at least 5 relocations 
week5 <- week5[, c("id", "x", "y")] #mcp needs only 3 columns
sp::coordinates(week5) <- ~x+y #specifying coordinates
sp::proj4string(week5) <-"+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45" #Lambert Conformal Conic

#calculating mcp for hour
hour.mcp <- mcp(hour5, percent = 100, unout = "km2")
hour.mcp
names(hour.mcp)[1]<-"id"
hour.mcpDF <- data.frame(hour.mcp) #dataframe of mcp output with all data columns
str(hour.mcpDF) ## 7 polygons created

#calculating mcp for day
day.mcp <- mcp(day5, percent = 100, unout = "km2")
day.mcp
names(day.mcp)[1]<-"id"
day.mcpDF <- data.frame(day.mcp) #dataframe of mcp output with all data columns
str(day.mcpDF) ## 13 polygons created

#calculating mcp for week
week.mcp <- mcp(week5, percent = 100, unout = "km2")
week.mcp
names(week.mcp)[1]<-"id"
week.mcpDF <- data.frame(week.mcp) #dataframe of mcp output with all data columns
str(week.mcpDF) ## 23 polygons created

#getting treatment names for hour
hour.mcpDF2 <- hour.mcpDF %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole")
  )

#getting treatment names for day
day.mcpDF2 <- day.mcpDF %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole")
  )

#getting treatment names for week
week.mcpDF2 <- week.mcpDF %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole")
  )

# mcp stats for hour
m1 <-lm(area~Treatment, data=hour.mcpDF2) #no difference
summary(m1)
anova(m1)

hour.mcpDF2  %>%
  group_by(Treatment) %>%
  summarise(
    areaAvg = mean(area, na.rm = T),
    sd = sd(area, na.rm=T),
    n = sum(!is.na(area)),
    se = sd / sqrt(n)
  )

# mcp stats for day
m2 <- lm(area~Treatment, data = day.mcpDF2)
summary(m2)
anova(m2) #  no difference

day.mcpDF2  %>%
  group_by(Treatment) %>%
  summarise(
    areaAvg = mean(area, na.rm = T),
    sd = sd(area, na.rm=T),
    n = sum(!is.na(area)),
    se = sd / sqrt(n)
  )

# mcp stats for week
m3 <- lm(area~Treatment, data = week.mcpDF2)
summary(m3)
anova(m3) #  no difference

week.mcpDF2  %>%
  group_by(Treatment) %>%
  summarise(
    areaAvg = mean(area, na.rm = T),
    sd = sd(area, na.rm=T),
    n = sum(!is.na(area)),
    se = sd / sqrt(n)
  )

# was going to plot mcp with leaflet but they overlap
# too much for it to be useful

#########################################################
####  Last Detection Point and number of days flying ####
#######################################################

###calculating the number of days flying
library(lubridate)
library(dplyr)

#going back to the raw tag data 

#there is one with a mistake that needs to be removed
tagsgdate2 <- tagsgdate %>%
  filter(as.Date(date) != as.Date("2021-09-07"))

daysFlown <- tagsgdate2 %>%
  group_by(id) %>%
  summarise(firstDay = min(date), lastDay = max(date)) %>%
  mutate(daysFlown = round(as.numeric(difftime(lastDay, firstDay, units = "days"))))

daysFlown.df <-data.frame(daysFlown)

daysFlown.df2<- daysFlown.df %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole"))

# summary stats

daysFlown.df2  %>%
  group_by(Treatment) %>%
  summarise(
    daysAvg = mean(daysFlown, na.rm = T),
    sd = sd(daysFlown, na.rm=T),
    n = sum(!is.na(daysFlown)),
    se = sd / sqrt(n)
  )


# is there a difference in treatments between days flown?

f1 <- glm(as.numeric(daysFlown) ~ Treatment, family = "poisson", 
          data = daysFlown.df2)
summary(f1) # overdispersed (754.99/39 = 19.36 > 1 therefore overdispersed)
library(AER)
dispersiontest(f1) #overdispersed, trying again with negative binomial

library(MASS)
f1.b <- glm.nb(as.numeric(daysFlown) ~ Treatment,  
               data = daysFlown.df2)
summary(f1.b) #dispersion is 1.3, very slighly overdispersed but should be fine
anova(f1.b, test = "Chisq") #not sig

######################################################3
###   extracting the last detection point     #########
##############################################################

landcover <-terra::rast("SpringLandcoverRaster.tif") #Rare landcover, made in overwintering

library(dplyr)
library(sf)

lastPoint <- daysFlown.df2 %>%
  left_join(tagsgdate2 %>% dplyr::select(id, x, y), by = "id")

#make last point spatial
lastPoint.sp <-st_as_sf(lastPoint, coords = c("x", "y"), crs = 4326)

#project to the same crs as landcover

lastPoint.proj <- st_transform(lastPoint.sp, crs = st_crs(landcover))

# Extract the data
extracted_data <- terra::extract(landcover, lastPoint.proj)

# Assign the second column
lastPoint.proj$landtype <- extracted_data[, 2]

lastPoint.df <- as.data.frame(lastPoint.proj)

lastPoint.df2 <-  lastPoint.df %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole"))

# is there a difference in last occurrence point landtype between treatments?
table_data <- table(as.factor(lastPoint.df2$landtype), lastPoint.df2$Treatment)

# Perform Chi-Square Test
l1 <- chisq.test(table_data)

l1.b <- fisher.test(table_data)

l1.c <- chisq.test(table_data, simulate.p.value = TRUE, B = 10000)

# View the results
print(l1)

#posthocs
residuals_table <- chisq.test(table_data)$stdres


library(ggplot2)

#plot of landcover types associated with last occurrence point
# raw data
ggplot(data=lastPoint.df2, aes(x=as.character(landtype), fill=Treatment)) +
  geom_bar(stat="count", position = position_dodge2(preserve="single")) +
  scale_y_continuous(name= "count of last occurrences", 
                     breaks = seq(0, max(table(lastPoint.df2$landtype)), 
                                  by = 50)) +
  scale_x_discrete(limits=c("3", "5", "6"),
                   labels=c( "early floral", "late floral", "low floral")) +
  scale_fill_manual(values = c("#3F612D", "#FFC800", "#49BEAA"),
                    name = "") +
  theme_classic() +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(size = 12),
        axis.title.y=element_text(size=14), 
        axis.text.y=element_text(size=12),
        legend.text=element_text(size=12)) 



#making a map of where the last points are located

library(sf)  # Ensure sf package is loaded

## Extract x and y from geometry
lastPoint.sp2 <- lastPoint.sp %>%
  mutate(x = st_coordinates(.)[,1],  # Extract X
         y = st_coordinates(.)[,2])  # Extract Y

lastPoint.sp3 <-  lastPoint.sp2 %>% 
  mutate(Treatment = case_when(
    startsWith(id, "C") ~ "control",
    startsWith(id, "I") ~ "imidacloprid",
    startsWith(id, "Y") ~ "cyantraniliprole"))


set.seed(42)  # For reproducibility
jitter_amount <- 0.0001  # Adjust this value if needed
lastPoint.sp4 <- lastPoint.sp3 %>%
  mutate(xAdj = x + runif(n(), -jitter_amount, jitter_amount),
         yAdj = y + runif(n(), -jitter_amount, jitter_amount))

# Update geometry using adjusted coordinates
# Update geometry using xAdj and yAdj
lastPoint.sp5 <- lastPoint.sp4 %>%
  mutate(geometry = st_sfc(purrr::map2(xAdj, yAdj, ~st_point(c(.x, .y))), 
                           crs = st_crs(lastPoint.sp4)))


library(leaflet)

treatpal <- colorFactor(palette= c("red", "yellow", "blue"),
                        levels = c("control", "imidacloprid", "cyantraniliprole"))

# **Update geometry with jittered coordinates**
lastPoint.sp6 <- st_as_sf(lastPoint.sp5, coords = c("xAdj", "yAdj"), 
                          crs = st_crs(lastPoint.sp5))

# Check if new geometry is applied
print(st_coordinates(lastPoint.sp6))


# Create the map with customized legend
leaflet(lastPoint.sp6) %>% 
  addTiles() %>%
  addProviderTiles('Esri.WorldImagery') %>%
  setView(lng = -80.352, lat = 43.377, zoom = 16) %>%
  addCircleMarkers(
    ~xAdj, ~yAdj, 
    stroke = FALSE, 
    color = ~treatpal(Treatment),
    fillOpacity = 0.5,
    radius = 8
  ) %>%
  addLegend(
    position = "topright",
    pal = treatpal, 
    values = ~Treatment, 
    opacity = 1,
    title = ""  # No title
  )


