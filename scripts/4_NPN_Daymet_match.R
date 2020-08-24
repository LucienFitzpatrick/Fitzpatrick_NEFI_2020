#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script is for cleaning up the NPN data to make it easier to work with to start 
# Inputs: Clean NPN observation dataframe from 3_NPN_clean.R script
#         Daymet_clean_data from 2_Daymet_download.R
# Outputs: Dataframe that has the weather calculations for all days of interest
#         Visualizations of the clean data
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
library(dplyr)
path.doc <- ("../data_processed/fall/")
path.fig <- ("../data_processed/fall/figures")
if(!dir.exists(path.fig)) dir.create(path.fig, recursive=T)

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

#Looking at the frequency of observations on different days
hist(dat.npn$day_of_year)

#Looking at the porportion of yes vs. no on every given day
png(filename= file.path(path.fig, paste0("Oak_Collection_Observation_Hist.png")))
ggplot(data = dat.npn)+
  geom_histogram(aes(x = day_of_year, group = color.clean, fill = color.clean), stat = "count")
dev.off()

#Checking proportion of yes vs no over time
dat.prop <- dat.npn[,c("day_of_year", "color.clean")]

freq <- as.data.frame(table(dat.npn$day_of_year))

freq.yes <- as.data.frame(table(dat.prop))
freq.yes <- freq.yes[freq.yes$color.clean == 1,]

freq$color <- freq.yes$Freq[match(freq$Var1, freq.yes$day_of_year)]

colnames(freq) <- c("day_of_year", "nObs", "nObs_yes")

for(i in 1:nrow(freq)){
  freq[i, "prop"] <- freq[i, "nObs_yes"]/freq[i, "nObs"]
}


png(filename= file.path(path.fig, paste0("Oak_Collection_Proportions.png")))
ggplot(data = freq)+
  ggtitle("Proportion of observed Yes values of fall color")+
  geom_point(aes(x = day_of_year, y = prop, group = 1))+
  geom_line(aes(x = day_of_year, y = prop, group = 1))+
  theme(axis.text.x = element_text(angle=70, size = 7))
dev.off()

#Looking at the days with no observations
all.doy <- sort(unique(dat.npn$day_of_year))

missed.doy <- setdiff(214:335, all.doy)

missed.doy

#end.doy <- setdiff(340:365, all.doy)


#Looking at invididuals trends
png(filename= file.path(path.fig, paste0("Oak_Collection_Individuals_2018.png")))
dat.2018 <- dat.npn[dat.npn$year == 2018 ,]
ggplot(data = dat.2018) +
  facet_wrap(~individual_id)+
  geom_line(aes(x = day_of_year, y = color.clean, color = as.character(year)))+
  scale_color_discrete()
dev.off()

png(filename= file.path(path.fig, paste0("Oak_Collection_Individuals_2019.png")))
dat.2019 <- dat.npn[dat.npn$year == 2019 ,]
ggplot(data = dat.2019) +
  facet_wrap(~individual_id)+
  geom_line(aes(x = day_of_year, y = color.clean, color = as.character(year)))+
  scale_color_discrete()
dev.off()


#Setting up trend lines for graphing
data_all = dat.npn[,c("day_of_year","color.clean")]
data_all = data_all[order(data_all$day_of_year),]
data_all$color.clean = as.numeric(as.character(data_all$color.clean))
data_all$day_of_year = as.numeric(as.character(data_all$day_of_year))

data_all_loess_10 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.10) # 10% smoothing span
data_all_loess_10 <- predict(data_all_loess_10) 
data_all_loess_30 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.30) # 30% smoothing span
data_all_loess_30 <- predict(data_all_loess_30) 

png(filename= file.path(path.fig, paste0("Oak_Collection_Loess_Trend")))
plot(data_all$day_of_year, data_all$color.clean, type="p", main="Loess Smoothing 2019 Data", xlab="Day of Year", ylab="Fall Color")
lines(data_all_loess_10, x=data_all$day_of_year, col="red", lwd = 2)
lines(data_all_loess_30, x=data_all$day_of_year, col="blue", lwd = 2)
dev.off()

