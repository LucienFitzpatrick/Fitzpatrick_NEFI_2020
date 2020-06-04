#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the initial data download from NPN
# Inputs: N/A
# Outputs: data frame of NPN observations
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#

#This MUST BE DOWNLOADED FROM GITHUB
#The CRAN version does not have the full capabilities needed for our project.
install.packages("devtools")
library('devtools')
devtools::install_github("usa-npn/rnpn")

library(rnpn)

#setting file path
dir.create("../data_processed/", recursive = T, showWarnings = F)


#Retrieving npn data
#rnpn packages has tools to show the corresponding id's for these queries. Request source is your name/affiliation
#Phenophase 498 is colored leaves
#network 720 is the morton arboretum
dat.npn <- npn_download_status_data(request_source='Morton Arboretum', years=c(2017:2019), 
                                    network_ids = c(720), phenophase_ids =c(498))


write.csv(dat.npn, "../data_processed/Arb_NPN_data_raw.csv", row.names=F)



