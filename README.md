# Biodiversity2018
repository of 2018 species richness analysis 

 -> Single_File_Script.R : Run's primary analyses utilized in paper submission. The script is well annotated
  1. Plot observed richness across domains
  2. Unify and merge environmental correlates across domains
  3. Train model to predict richness
  5. Visualize residuals
  
 -> dev > This folder contains additional scripts utilized in paper analyses
  a. Some preprocessing of time series of temperature and solar insolation
  b. The development and testing of variable sensitivity/importance for paper NN model
  c. Log linear var imp slope testing 
  d. Supplemental plotting
  
 __data__ > paths within the Single_File_Script.R refer to data within in this folder. Within this script path directories may need to be changed to reflect the home directory of the user. Data within this folder is either biodiversity data, primary production time series, insolation, precipitation, oxygen, or physical relief. 

ClintonExport_sept22 <- is Clinton Jenkins managed terrestrial richness data
ECMWF <- European center for moderate range weather forecasting SST/SAT
Gabriel <- Richness and environ data relative to Gabriel Reygondeau contributions
land shapefiles <- as labeled
marine diversity <- deprecated marine diversity data from Clinton Jenkins and Aquamaps
NEO data <- NDVI/CHla from NASA NEO FTP servers
terrestrial diversity <- deprecated diversity data from Clinton Jenkins

