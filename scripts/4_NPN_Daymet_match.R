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

library(ggplot2)

hist(dat.npn$day_of_year)

ggplot(data = dat.npn)+
  geom_histogram(aes(x = day_of_year), stat = "count")


all.doy <- sort(unique(dat.npn$day_of_year))

missed.doy <- setdiff(214:335, all.doy)

missed.doy

#end.doy <- setdiff(340:365, all.doy)


#Looking at invididuals trends
dat.2018 <- dat.npn[dat.npn$year == 2018 ,]
ggplot(data = dat.2018) +
  facet_wrap(~individual_id)+
  geom_line(aes(x = day_of_year, y = color.clean, color = as.character(year)))+
  scale_color_discrete()


data_all = dat.npn[,c("day_of_year","color.clean")]
data_all = data_all[order(data_all$day_of_year),]
data_all$color.clean = as.numeric(as.character(data_all$color.clean))
data_all$day_of_year = as.numeric(as.character(data_all$day_of_year))

data_all_loess_10 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.10) # 10% smoothing span
data_all_loess_10 <- predict(data_all_loess_10) 
data_all_loess_30 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.30) # 10% smoothing span
data_all_loess_30 <- predict(data_all_loess_30) 
plot(data_all$day_of_year, data_all$color.clean, type="p", main="Loess Smoothing 2019 Data", xlab="Day of Year", ylab="Fall Color")
lines(data_all_loess_10, x=data_all$day_of_year, col="red", lwd = 2)
lines(data_all_loess_30, x=data_all$day_of_year, col="blue", lwd = 2)
