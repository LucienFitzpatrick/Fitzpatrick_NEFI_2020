# NEFI
  This repository serves to house scripts and data used during the NEFI 2020 summer course.
  The overall purpose of this project is to create a near term iterative state space model that can take historical NPN and weather data combined with weather forecasts and live updated phenology monitoring to regularly update a model predicting the intensity of fall color for collections of trees.

# Our data

## National Phenology Network Observations 

  This data comes from a project that monitors oak fall color phenology within The Morton Arboretum Oak Collection using protocols from the National Phenology Network (usanpn.org). Our monitoring was started before the monitoring season of 2017 through present. For this project we are using 2017 to 2019 because 2020 data hasn’t been collected for fall yet. Observations were made by volunteers that had received training from our volunteer coordinator Brendon Reidy. Observations were made 5-9 days apart starting in the spring before bud burst and continuing through the fall until all individuals have dropped their leaves. Observation schedules are not uniform so while an individual tree is monitored every 5-9 days there are likely observations being taken on any given day. Individuals do not need to be monitored in the fall/winter after leaf drop. Species have been recorded in Nature's Notebook, the National Phenology Network database, allowing the data to be freely available for anyone who wishes to use it. More details on NPN observation protocols can be found here: https://www.usanpn.org/nn/guidelines

## Daymet
Data Overview

  This data is retrieved from the daymet project run by NASA. The Daymet dataset provides gridded estimates of daily weather parameters. Seven surface weather parameters are available at a daily time step, 1 km x 1 km spatial resolution, with a North American spatial extent. Access to the Daymet dataset is available from the ORNL DAAC through a variety of tools and formats allowing a rich resource of daily surface meteorology.


# Workflow/Script order

**1_NPN_download.R** - This will create a dataframe containing the raw data from www.usanpn.org for the specified parameters. This requires the developmental version of the "rnpn" package that must be downloaded directly from github.

**2_Daymet_download.R** - This will take the dataframe of npn observations and provide a nested list containing daily weather data for every location for the range of years.

**3_NPN_clean.R** - This will take the raw dataframe created from 1_NPN_download.R and to both organize the data and clean away potentially misentered observation data. It does this by checking all "YES" values for the presence of leaves or colored leaves and making sure they also have either a previous or future "YES" observation. This is because it wouldn't make sense for a tree to have leaves for one week and then not anymore (This is an assumption that proves accurate for our Oaks and most species but is not neccessairly universal)

**4_xx_model.R** -These are the various models that we currently have. They follow a format where the title explains what the major difference is. They are all state space models and all but the starting models use covariates in their title name

**4_Starting_model.R** - Random walk state space null model. Early model no longer used

**4_Starting_model_increaseOnly.R** - Random walk state space null model that includes constraints so that they proportion can only increase. Early base model

**4_LengthOfDay.R** - State space model using day length as a covariate. Model assumes that once color occurs it doesn't go away. It uses daylength as calculated from daymet

**4_Precip_model.R** -State space model using precipitaiton as a covariate. Model assumes that once color occurs it doesn't go away. It uses daily precipitation as calculated from daymet

**4_tmin.R** - State space model using minimum temperature as a covariate. Model assumes that once color occurs it doesn't go away. It uses minimum temperature as calculated from daymet

**4_multi_covars.R** - State space model using multiple covariates. Currently non-functional

**5_forecast.R** - This will take one of the single covariate models and allow you to forecast using that data. Minor modifications will be required depending on which model is being used for the forecast

**6_iterative** - This will take one of the single covariate models and allow you to create an iterative 2018 hindcast

# Output Files

## **Arb_Quercus_NPN_data_leaves_raw.csv**

Created by 1_NPN_download.R. A raw data frame of all NPN leaf observations for all monitored Oaks in the Morton Arboretum from 2018-2019

**Column Description**

Variable | Description | Unit
-------- | ----------- | ----
observation_id | ID tag for individual observation | Integer 
update_datetime | Date and time of last update to observation | Date: YYYY-MM-DD HH:MM:SS
site_id | ID tag for site of observation | Integer
latitude | Decimal latitude of observation with datum WGS84 | Integer
longitude | Decimal longitude of observation with datum WGS84 | Integer
elevation_in_meters | Elevation of observation | meters
state | USA state where observation occurred | Character string
species_id | ID tag for species of individual observed | Integer
genus | Genus name of individual observed | Character string
species | Species name of individual observed | Character string
common_name | Common name of individual observed | Character string
kingdom | Kingdom of individual observed | Character string
individual_id | ID tag of individual observed | Integer
phenophase_id | ID tag of phenophase being observed | Integer
phenophase_description | Description of phenophase | Character string
observation_date | Date of observation | Date: YYYY-MM-DD
day_of_year | Day of year of observation | Integer
phenophase_status | Values denoting “Yes (1)” “No (0)” or “Unsure (-1)” for the phenophase being observed| 1 for Yes 0 for No -1 for Unsure
intensity_category_id | ID for category of intensity value | Integer
intensity_value | Value of intensity of phenophase | Percentage range
abundance_value | Measure used for animal abundance| Percentage range

## **Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv** 

Created by 3_NPN_clean.R. A data frame of all NPN observations for presence of leaves and colored leaves that have been cleaned for potential false "Yes" values.

**Column Description**
Variable | Description | Unit
-------- | ----------- | ----
individual_id | ID tag of individual observed | Integer
species_id | ID tag for species of individual observed | Integer
genus | Genus name of individual observed | Character string
observation_date | Date of observation | Date: YYYY-MM-DD
year | Year of measurements | Integer YYYY
day_of_year | Day of year of observation | Integer
leaves | If leaves were observed on the tree in our raw data | Binary (Yes/No)
color | If colored leaves were observed on the tree in our raw data | Binary (Yes/No)
leaf.clean | If leaves were observed on the tree AND the previous or future value is also yes | Binary (Yes/No)
color.clean | If colored leaves were observed on the tree AND the previous or future value is also yes | Binary (Yes/No)

## **Daymet_data_raw.csv** 

Created by 2_Daymet_download.R. A raw data frame of the daily weather metrics at the Morton Arboretum for every day from 2018-2019. 

Note: All Daymet years are 1 –365 days, including leap years. The Daymet database includes leap-days. Values for December31 are discarded from leap years to maintain a 365-day year.

**Column description**

Variable | Description | Unit
-------- | ----------- | ---- 
year | Year of measurements | Integer YYYY
yday | Day of year of measurements Values ranging from 1-365 | Integer
dayl..s. | Duration of the daylight period for the day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon | (s/day)
prcp..mm.day. | Daily total precipitation, sum of all forms converted to water-equivalent | (mm/day)
srad..W.m.2. | Incident shortwave radiation flux density, taken as an average over the daylight period of the day | (W/m^2)
swe..kg.m.2. | Snow water equivalent. The amount of water contained within the snowpack | kg/m^2)
tmax..deg.c. | Daily maximum 2-meter air temperature | (degrees C)
tmin..deg.c. | Daily minimum 2-meter air temperature | (degrees C)
vp..Pa. | Water Vapor Pressure. Daily average partial pressure of water vapor. | Pascals
