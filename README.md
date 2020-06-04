# NEFI
Folder for scripts resulting from the NEFI summer course

Workflow/Script order

1_NPN_download.R - This will create a dataframe containing the raw data from www.usanpn.org for the specified parameters. This requires the developmental version of the "rnpn" package that must be downloaded directly from github.

2_Daymet_download.R - This will take the dataframe of npn observations and provide a nested list containing daily weather data for every location for the range of years.
