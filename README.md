# NEFI
  This repository serves to house scripts and data used during the NEFI 2020 summer course.
  The overall purpose of this project is to create a near term iterative state space model that can take historical NPN and weather data combined with weather forecasts and live updated phenology monitoring to regularly update a model predicting the intensity of fall color for collections of trees.

# Our data

  All of our data has it's own README in the data folder but the general scope of the data gathering and the institutions behind them is outline here
  
  
## National Phenology Network Observations 

  This data comes from a project that monitors oak fall color phenology within The Morton Arboretum Oak Collection using protocols from the National Phenology Network (usanpn.org). Our monitoring was started before the monitoring season of 2017 through present. For this project we are using 2017 to 2019 because 2020 data hasn’t been collected for fall yet. Observations were made by volunteers that had received training from our volunteer coordinator Brendon Reidy. Observations were made 5-9 days apart starting in the spring before bud burst and continuing through the fall until all individuals have dropped their leaves. Observation schedules are not uniform so while an individual tree is monitored every 5-9 days there are likely observations being taken on any given day. Individuals do not need to be monitored in the fall/winter after leaf drop. Species have been recorded in Nature's Notebook, the National Phenology Network database, allowing the data to be freely available for anyone who wishes to use it. More details on NPN observation protocols can be found here: https://www.usanpn.org/nn/guidelines

## Daymet
Data Overview

  This data is retrieved from the daymet project run by NASA. The Daymet dataset provides gridded estimates of daily weather parameters. Seven surface weather parameters are available at a daily time step, 1 km x 1 km spatial resolution, with a North American spatial extent. Access to the Daymet dataset is available from the ORNL DAAC through a variety of tools and formats allowing a rich resource of daily surface meteorology.


# Workflow/Script order

**1_NPN_download.R** - This will create a dataframe containing the raw data from www.usanpn.org for the specified parameters. This requires the developmental version of the "rnpn" package that must be downloaded directly from github.

**2_Daymet_download.R** - This will take the dataframe of npn observations and provide a nested list containing daily weather data for every location for the range of years.

**3_NPN_clean.R** - This will take the raw dataframe created from 1_NPN_download.R and to both organize the data and clean away potentially misentered observation data. It does this by checking all "YES" values for the presence of leaves or colored leaves and making sure they also have either a previous or future "YES" observation. This is because it wouldn't make sense for a tree to have leaves for one week and then not anymore (This is an assumption that proves accurate for our Oaks and most species but is not neccessairly universal)

# Output Files

Inside the data folder in this repository is a read me containing information on the values in these various data frames

**Arb_Quercus_NPN_data_leaves_raw.csv** - Created by 1_NPN_download.R. A raw data frame of all NPN leaf observations for all monitored Oaks in the Morton Arboretum from 2018-2019

**Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv** - Created by 3_NPN_clean.R. A data frame of all NPN observations for presence of leaves and colored leaves that have been cleaned for potential false "Yes" values.

**Daymet_data_raw.csv** - Created by 2_Daymet_download.R. A raw data frame of the daily weather metrics at the Morton Arboretum for every day from 2018-2019. 
