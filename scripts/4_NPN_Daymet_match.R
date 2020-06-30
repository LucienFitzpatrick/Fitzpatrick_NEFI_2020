#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script is for cleaning up the NPN data to make it easier to work with to start 
# Inputs: Clean NPN observation dataframe from 3_NPN_clean.R script
#         Daymet_clean_data from 2_Daymet_download.R
# Outputs: Dataframe that has the weather calculations for all days of interest
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), na.strings = "-9999")
lat.calc <- read.csv(file.path(path.doc, file = "Daymet_clean_data.csv"))

dat.npn <- dat.npn[dat.npn$day_of_year >= min(lat.calc$yday),]


dat.npn$CDD5.cum <- lat.calc$CDD5.cum[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$CDD0.cum <- lat.calc$CDD0.cum[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$NCD <- lat.calc$NCD[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$Precip <- lat.calc$prcp..mm.day.[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$Precip.cum <- lat.calc$Precip.cum[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$LOD <- lat.calc$dayl..s.[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$tmax <- lat.calc$tmax..deg.c.[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$tmin <- lat.calc$tmin..deg.c.[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$tmean <- lat.calc$TMEAN[match(dat.npn$observation_date, lat.calc$Date)]


# Save dat.npn 
write.csv(dat.npn, file.path(path.doc, "Fall_Phenology_data.csv"), row.names=F)

