#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script is for cleaning up the NPN data to make it easier to work with to start 
# Inputs: Data frame of npn data downloaded in the 1_NPN_download script
# Outputs: data frame of clean NPN observations for colored leaves and leaves
#          Visualziations of the NPN observations
# Notes: Doesn't currently properly clean "breaking leaf buds" or "falling leaves" so they are excluded in the end
#-----------------------------------------------------------------------------------------------------------------------------------#

path.doc <- ("../data_processed/fall/")

library(ggplot2)
# ------------------------------
# Cleaning and summarizing data
# ------------------------------
dat.npn <- read.csv(file.path(path.doc, "Arb_Quercus_NPN_data_leaves_raw.csv"), na.strings = "-9999")
dat.npn$genus <- as.factor(dat.npn$genus)
dat.npn$species <- as.factor(dat.npn$species)
dat.npn$phenophase_description <- factor(dat.npn$phenophase_description, levels=c("Breaking leaf buds", "Leaves", "Colored leaves", "Falling leaves"))
dat.npn$observation_date <- as.Date(dat.npn$observation_date)
dat.npn$year <- lubridate::year(dat.npn$observation_date)
dat.npn$pheno.stat <- as.factor(car::recode(dat.npn$phenophase_status, "'0'='No'; '1'='Yes'; '-1'='Unsure'"))
dat.npn$intensity_value <- dat.npn$intensity_value
summary(dat.npn)

# Checking to see if we have enough intensity data to be useful
summary(dat.npn[dat.npn$phenophase_status==1 & dat.npn$phenophase_description=="Colored leaves",]) # 176/576 NAs ==> 70% of obs have intensity... could be good.


# Doing a quick graph:
day.labels <- data.frame(Date=seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by="month"))
day.labels$yday <- lubridate::yday(day.labels$Date)
day.labels$Text <- paste(lubridate::month(day.labels$Date, label=T), lubridate::day(day.labels$Date))
summary(day.labels)

if(!dir.exists("figures")) dir.create("figures")
png("figures/LeafPhenophases.png", heigh=4, width=6, units="in", res=180)
ggplot(data=dat.npn[dat.npn$pheno.stat=="Yes",]) +
  ggtitle("The Morton Arboretum, Oak Collection") + 
  facet_grid(year ~ phenophase_description) +
  geom_histogram(aes(x=day_of_year, y=..count.., fill=phenophase_description), binwidth=7) +
  scale_x_continuous(name="Day of Year", expand=c(0,0), breaks=day.labels$yday[seq(2,12, by=3)], label=day.labels$Text[seq(2,12, by=3)]) +
  scale_y_continuous(name="# Observations", expand=c(0,0)) +
  scale_fill_manual(values=c("springgreen3", "green4", "orange3", "tan4")) +
  guides(fill=F) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=-45, hjust=0)) 
dev.off()

#Converting less than 5% to be a no
for(i in 1:nrow(dat.npn)){
  dat.npn[i, "intensity_value"] <- ifelse(is.na(dat.npn[i, "intensity_value"]), 0, dat.npn[i, "intensity_value"])
  if(dat.npn[i, "intensity_value"]== "Less than 5%"){
    dat.npn[i, "pheno.stat"] <-  "No"
  }
}


# # Turning the data into a "wide" format right now to help
# dat.wide <- reshape2::recast(data=dat.npn, phenophase_description ~ individual_id + species_id + genus + observation_date + day_of_year)
# dat.wide <- tidyr::spread(data=dat.npn, key=phenophase_description, value=pheno.stat)
# summary(dat.wide)

dat.wide <- reshape2::dcast(data=dat.npn, individual_id + species_id + site_id + genus + observation_date + year+ day_of_year ~ phenophase_description, value.var="pheno.stat")

names(dat.wide)[(ncol(dat.wide)-3):ncol(dat.wide)] <- c("buds", "leaves", "color", "falling")
dat.wide$individual_id <- as.factor(dat.wide$individual_id)
dat.wide$buds <- as.factor(dat.wide$buds)
dat.wide$leaves <- as.factor(dat.wide$leaves)
dat.wide$color <- as.factor(dat.wide$color)
dat.wide$falling <- as.factor(dat.wide$falling)
summary(dat.wide)

ggplot(data=dat.wide) +
  facet_grid(.~year) +
  geom_tile(aes(x=day_of_year, y=individual_id, fill=color)) +
  theme(axis.text.y=element_blank())


# Cleaning up some of our spotty observations
#  Only count it as colored if either the next X or previous X have been "yes"
WIN <- 3 # Window to check for QAQC
for(IND in unique(dat.wide$individual_id)){
  for(YR in unique(dat.wide$year[dat.wide$individual_id==IND])){
    row.ind <- which(dat.wide$individual_id==IND & dat.wide$year==YR)
    
    dat.ind <- dat.wide[row.ind,]
    obs.dates <- unique(dat.ind$observation_date)
    for(OBS in seq_along(obs.dates)){
      obs.leaf <- dat.ind$leaves[OBS]
      obs.color <- dat.ind$color[OBS]
      
      if(OBS>WIN & OBS<length(obs.dates)-WIN){
        leaf.clean <- all(dat.ind$leaves[OBS:(OBS+WIN)]==obs.leaf) | all(dat.ind$leaves[(OBS-WIN):(OBS)]==obs.leaf)
        color.clean <- all(dat.ind$color[OBS:(OBS+WIN)]==obs.color) | all(dat.ind$color[(OBS-WIN):(OBS)]==obs.color)
      } else if(OBS<length(obs.dates)-WIN) {
        leaf.clean <- all(dat.ind$leaves[OBS:(OBS+WIN)]==obs.leaf)
        color.clean <- all(dat.ind$color[OBS:(OBS+WIN)]==obs.color)
      } else {
        leaf.clean <- all(dat.ind$leaves[(OBS-WIN):(OBS)]==obs.leaf) | obs.leaf=="No"
        color.clean <- all(dat.ind$color[(OBS-WIN):(OBS)]==obs.color) | obs.color=="No"
      }
      
      # Enter the cleaned value -- if the observed makes sense, use it; otherwise use the opposite
      dat.ind[OBS,"leaf.clean"] <- ifelse(leaf.clean, paste(obs.leaf), ifelse(obs.leaf=="Yes", "No", "Yes"))
      dat.ind[OBS,"color.clean"] <- ifelse(color.clean, paste(obs.color), ifelse(obs.color=="Yes", "No", "Yes"))
    } # End obs
    dat.wide[row.ind,c("leaf.clean", "color.clean")] <- dat.ind[,c("leaf.clean", "color.clean")]
  } # End year

} # End individual

# ggplot(data=dat.ind) +
  # geom_tile(aes(x=day_of_year, y=individual_id, fill=leaf.clean)) +
  # theme(axis.text.y=element_blank())

png("figures/Color_Cleaned.png", heigh=6, width=4, units="in", res=180)
ggplot(data=dat.wide) +
  facet_grid(.~year) +
  geom_tile(aes(x=day_of_year, y=individual_id, fill=color.clean)) +
  scale_x_continuous(name="Day of Year", expand=c(0,0), breaks=day.labels$yday[seq(2,12, by=3)], label=day.labels$Text[seq(2,12, by=3)]) +
  theme(legend.position="top", 
        axis.text.y=element_blank())
dev.off()

png("figures/Leaf_Cleaned.png", heigh=6, width=4, units="in", res=180)
ggplot(data=dat.wide) +
  facet_grid(.~year) +
  geom_tile(aes(x=day_of_year, y=individual_id, fill=leaf.clean)) +
  scale_x_continuous(name="Day of Year", expand=c(0,0), breaks=day.labels$yday[seq(2,12, by=3)], label=day.labels$Text[seq(2,12, by=3)]) +
  theme(legend.position="top", 
        axis.text.y=element_blank())
dev.off()

# Save a version of the individual-level cleaned leaf & color data
# Note: haven't cleaned up bud break or falling leaves
write.csv(dat.wide[,names(dat.wide)[!names(dat.wide) %in% c("buds", "falling")]], file.path(path.doc, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), row.names=FALSE)
# ------------------------------------------------------------------#
#THIS SECTION WILL EVENTUALLY CHANGE
#THIS section makes it so fall color never drops which makes modeling easier for now
#-------------------------------------------------------------------#

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
