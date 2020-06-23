#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the initial data download from daymet
# Inputs: NPN phenology observation file obtain from 1_NPN_download in this repository
# Outputs: Daymet data for every day of every year in every location the observations come from
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#

dat.npn <- read.csv("../data_processed/Arb_NPN_data_raw.csv")
#Setting the points to download the daymet data from
path.doc <- ("../data_processed/")

dat.npn$Year <- lubridate::year(dat.npn$observation_date)

dat.npn$Site <- 'NPN'


# Creating a point list and time range that matches your MODIS dataset
# Note: This will probably change down the road
NPN.pts <- aggregate(Year~Site+latitude+longitude, data=dat.npn, 
                     FUN=min)
names(NPN.pts)[4] <- "yr.start"
NPN.pts$yr.end <- aggregate(Year~Site+latitude+longitude, data=dat.npn, 
                            FUN=max)[,4]
NPN.pts


#Writing the csv file of lat and longs because daymetr batch function needs to read a file instead of a dataframe
write.csv(NPN.pts, file.path(path.doc, "NPN_points.csv"), row.names=FALSE)


setwd(path.doc)
#Downloading all of the damet data for each point. Internal =TRUE means it creates a nested list. Set false to actually download a file
lat.list <- daymetr::download_daymet_batch(file_location = file.path(path.doc, "NPN_points.csv"),
                                           start = min(NPN.pts$yr.start),
                                           end = max(NPN.pts$yr.end),
                                           internal = T)

#removing failed downloads 
lat.list <- lat.list[sapply(lat.list, function(x) is.list(x))]

#This will only work for arb specific data. This will need to become a loop for multiple locations
lat.df <- as.data.frame(lat.list[[2]]$data)
lat.df$latitude <- lat.list[[2]]$latitude
lat.df$longitude <- lat.list[[2]]$longitude

write.csv(lat.df, file.path(path.doc, file = "Daymet_data_raw.csv"), row.names=FALSE)
