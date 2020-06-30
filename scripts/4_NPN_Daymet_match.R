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

# Save dat.comb 
write.csv(dat.npn, file.path(path.doc, "Fall_Phenology_data.csv"), row.names=F)
dat.npn$GDD5.cum <- lat.calc$GDD5.cum[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$GDD0.cum <- lat.calc$GDD0.cum[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$NCD <- lat.calc$NCD[match(dat.npn$observation_date, lat.calc$Date)]
dat.npn$GTmean <- GTmean[match(dat.npn$observation_date, lat.calc$Date)]

