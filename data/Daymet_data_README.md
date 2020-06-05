# **Daymet_data_raw.csv**

Note: All Daymet years are 1 –365 days, including leap years. The Daymet database includes leap-days. Values for December31 are discarded from leap years to maintain a 365-day year.

## **Column description**

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
