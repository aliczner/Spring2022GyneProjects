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
