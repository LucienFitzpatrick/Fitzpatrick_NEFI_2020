**Arb_Quercus_NPN_data_leaves_raw.csv**

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


**Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv**

The leaf.clean and color.clean values are calculated by checking our raw "Yes" observations and making sure that either the previous or future value also has a "Yes". This is in case of misentered values becasue if a tree was said to have leaves one observations without previous or future observations also seeing leaves it is assumed a mistake.

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
