#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the initial data download from daymet
# Inputs: NPN phenology observation file obtain from 1_NPN_download in this repository
# Outputs: Daymet data for every day of every year in every location the observations come from
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Arb_Quercus_NPN_data_leaves_raw.csv"))
#Setting the points to download the daymet data from


path.hub <- "C:/Users/lucie/Documents/GitHub/"

dat.npn$Year <- lubridate::year(dat.npn$observation_date)


# Creating a point list and time range that matches your MODIS dataset
# Note: This will probably change down the road
NPN.pts <- aggregate(Year~site_id+latitude+longitude, data=dat.npn, 
                     FUN=min)
names(NPN.pts)[4] <- "yr.start"
NPN.pts$yr.end <- aggregate(Year~site_id+latitude+longitude, data=dat.npn, 
                            FUN=max)[,4]
NPN.pts


#Writing the csv file of lat and longs because daymetr batch function needs to read a file instead of a dataframe
write.csv(NPN.pts, file.path(path.doc, "NPN_points.csv"), row.names=FALSE)


#Downloading all of the damet data for each point. Internal =TRUE means it creates a nested list. Set false to actually download a file
lat.list <- daymetr::download_daymet_batch(file_location = file.path(path.doc, "NPN_points.csv"),
                                           start = min(NPN.pts$yr.start),
                                           end = max(NPN.pts$yr.end),
                                           internal = T)

#removing failed downloads 
lat.list <- lat.list[sapply(lat.list, function(x) is.list(x))]

names(lat.list) <- NPN.pts$site # Giving the different layers of the list the site names they correspond to

# Creating a new simplified list that won't make Christy cranky
list.met <- list()
for(i in seq_along(lat.list)){
  list.met[[i]] <- data.frame(site=NPN.pts$site_id[i], latitude=NPN.pts$latitude[i], longitude=NPN.pts$longitude[i], lat.list[[i]]$data)
}
names(list.met) <-  NPN.pts$site.id


#Reading in our weather calculation function
source(file.path(path.hub, "Phenology_Forecasting/scripts/weather_calc.R"))

list.met<- lapply(list.met, weather_calc)

lat.calc <- dplyr::bind_rows(list.met)

write.csv(lat.calc, file.path(path.doc, file = "Daymet_clean_data.csv"), row.names=FALSE)


