path.doc <- ("../data_processed/")

dat.npn <- read.csv(file.path(path.doc, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), na.strings = "-9999")

dat.npn$color.clean <- as.numeric(car::recode(dat.npn$color.clean, "'No'='0'; 'Yes'='1'; 'NA'='-1'"))

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





write.csv(dat.npn, file.path(path.doc, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), row.names=FALSE)




