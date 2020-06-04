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

ystart <- min(dat.npn$Year)

#make sure the yend of the data matches what you enter. Sometimes daymet truncates and this varibale will become wrong later in the script
yend <- max(dat.npn$Year)

pointsfile <- "npn_points.csv"

#Subsetting to only include lat and long (and for now the first rows to make testing easier)
q.lat <- dat.npn[, c("latitude", "longitude")]

#creating a proxy "site" column because the batch function needs it
q.lat$site <- "Daymet"
q.lat <- q.lat[,c(3,1,2)]
q.lat <- unique(q.lat)

#Writing the csv file of lat and longs because batch function needs to read a file instead of a dataframe
write.csv(q.lat, file.path(path.doc, file = pointsfile), row.names=FALSE)


setwd(path.doc)
#Downloading all of the damet data for each point. Internal =TRUE means it creates a nested list. Set false to actually download a file
lat.list <- daymetr::download_daymet_batch(file_location = pointsfile,
                                           start = ystart,
                                           end = yend,
                                           internal = T)

#removing failed downloads 
lat.list <- lat.list[sapply(lat.list, function(x) is.list(x))]
