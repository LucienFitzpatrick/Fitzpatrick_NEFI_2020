path.hub <- "C:/Users/lucie/Documents/GitHub/NEFI/data/"

dat.npn <- read.csv(file.path(path.hub, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), na.strings = "-9999")


for(YR in unique(dat.npn$year)){
  dat.YR <- dat.npn[dat.npn$year == YR,]
  
  for(i in unique(dat.YR$individual_id)){
    dat.tmp <- dat.YR[dat.YR$individual_id == i, ]
    count <-  0
    
    for(k in 1:nrow(dat.tmp)){
      
      if (dat.tmp[k, "color.clean"] == 1){
        count <- count + 1
      }
      
      if (count > 3 & dat.tmp[k, "day_of_year"] > 250){
        dat.tmp[k, "color.clean"] <- 1
      }
    
    }
  dat.npn[(dat.npn$year == YR & dat.npn$individual_id == i), "color.full"] <- dat.tmp$color.clean
  }
}


write.csv(dat.npn, file.path(path.hub, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), row.names=FALSE)




