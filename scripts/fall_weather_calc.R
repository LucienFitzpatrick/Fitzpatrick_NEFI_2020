#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the initial data download from daymet
# Inputs: NPN phenology observation file obtain from 1_NPN_download in this repository
# Inputs: Lubridate package
# Outputs:This function will take a data frame of daily weather data and produce the following summary statistics
#         CDD5 = Chilling degree days at 5 degrees C 
#         CDD0 = Chilling degree days at 0 degrees C
#         NCD = Number of chilling days 
# Notes: The defaults for this funcion are
#       Start of fall tracking for chilling effects (August 1st)                     f_start = 213
#
#-----------------------------------------------------------------------------------------------------------------------------------#
# Calculating the Tmean for the growing season of that year
fall_weather_calc <- function(met.all, f_start = 213){
  
  
  #Calculating mean temperature, growing degree days using 5C, and gorwing degree days using 0c
  met.all$TMEAN <- (met.all$tmax..deg.c. + met.all$tmin..deg.c.)/2
  met.all$CDD5 <- abs(ifelse(met.all$TMEAN<5, met.all$TMEAN-5, 0))
  met.all$CDD0 <- abs(ifelse(met.all$TMEAN<0, met.all$TMEAN, 0))
  
  #Creating empty columns to fill in the loop
  met.all$CDD5.cum <- NA
  met.all$CDD0.cum <- NA
  met.all$NCD <- NA
  met.all$Precip.cum <- NA
  
  #Setting the beginning and end (using julian yday) of the growing season for our Growing season mean
  #g_start <- 1
  #g_end <- 120
  met.all <- met.all[met.all$yday >= f_start,]
  
  
  # Calculate the cumulative growing degree days for each day/year
  for(YR in unique(met.all$year)){
    dat.tmp <- met.all[met.all$year==YR, ]
    #dat.gtmean <- dat.tmp[(dat.tmp$yday>=g_start & dat.tmp$yday<=g_end), ]
    #gtmean <- mean(dat.gtmean$TMEAN, na.rm = TRUE)
    dat.tmp$Date <- as.Date(paste(dat.tmp$year, dat.tmp$yday, sep="-"), format="%Y-%j")
    cdd5.cum=0; cdd0.cum=0
    d5.miss = 0; d0.miss=0
    ncd = 0
    precip.cum = 0
    for(i in 1:nrow(dat.tmp)){
      if(is.na(dat.tmp$CDD5[i]) & d5.miss<=7){ #YOU CHANGED THIS TO 7 FOR NOW BUT CHANGE BACK
        d5.miss <- d5.miss+1 # Let us miss up to 3 consecutive days
        cdd5.cum <- cdd5.cum+0
      } else {
        d5.miss = 0 # reset to 0
        cdd5.cum <- cdd5.cum+dat.tmp$CDD5[i] 
      }
      
      if(is.na(dat.tmp$CDD0[i]) & d0.miss<=7){ #YOU CHANGED THIS TO 7 FOR NOW BUT CHANGE BACK
        d0.miss <- d0.miss+1 # Let us miss up to 3 consecutive days
        cdd0.cum <- cdd0.cum+0
      } else {
        d0.miss = 0 # reset to 0
        cdd0.cum <- cdd0.cum+dat.tmp$CDD0[i] 
      }
      if(!is.na(dat.tmp$TMEAN[i]) & dat.tmp$TMEAN[i] < 0){
        ncd <- ncd + 1
      }
      precip.cum <- precip.cum + dat.tmp$prcp..mm.day.[i]
      
      dat.tmp[i,"CDD5.cum"] <- cdd5.cum
      dat.tmp[i,"CDD0.cum"] <- cdd0.cum
      dat.tmp[i, "NCD"] <- ncd
      dat.tmp[i, "Precip.cum"] <- precip.cum
      #dat.tmp[i, "GTmean"] <- gtmean
    }
    met.all[met.all$year==YR, "CDD5.cum"] <- dat.tmp$CDD5.cum
    met.all[met.all$year==YR, "CDD0.cum"] <- dat.tmp$CDD0.cum
    met.all[met.all$year==YR, "NCD"] <- dat.tmp$NCD
    met.all[met.all$year==YR, "Precip.cum"] <- dat.tmp$Precip.cum
    #met.all[met.all$year==YR, "GTmean"] <- dat.tmp$GTmean
    met.all[met.all$year==YR, "Date"] <- dat.tmp$Date
  }
  return(met.all)
}
