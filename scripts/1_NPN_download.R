#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the initial data download from NPN
# Inputs: N/A
# Outputs: data frame of NPN observations
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#

#This MUST BE DOWNLOADED FROM GITHUB: https://github.com/usa-npn/rnpn
#The CRAN version does not have the full capabilities needed for our project.
# install.packages("devtools")
# library('devtools')
# devtools::install_github("usa-npn/rnpn")
library(rnpn)

#setting file path
# path.doc <- "C:/Users/lucie/Documents/NPN_data/"
path.doc <- "NPN_data"
if(! dir.exists(path.doc)) dir.create(path.doc)

# ------------------------------
# 1. Retrieving npn data
# ------------------------------
#rnpn packages has tools to show the corresponding id's for these queries. Request source is your name/affiliation
#Phenophases: 
#  - breaking leaf buds = 371
#  - leaves = 483
#  - colored leaves = 498
#  - falling leaves = 471
#network 720 is The Morton Arboretum
# Station  26202 is the Oak Collection
dat.npn <- npn_download_status_data(request_source='Morton Arboretum', years=c(2018:2019), 
                                    network_ids = c(720), station_ids = c(26202), phenophase_ids =c(371, 483, 498, 471))

write.csv(dat.npn, file.path(path.doc, file = "Arb_Quercus_NPN_data_leaves_raw.csv"), row.names=FALSE)
# ------------------------------





