library(raster)
library(dplyr)
library(spatial.tools)

sea_mondo = stack( "/Users/tgagne/Biodiversity2018/data/Rdata/sea_mondo_Jan15.grd" )

ECMWF_SAT <- brick("/Users/tgagne/Biodiversity2018/data/ECMWF/SAT07_17.nc") - 273.15 # K to C
SAT_mean  <- calc(ECMWF_SAT, fun = mean )
SAT_sd    <- calc(ECMWF_SAT, fun = sd )
SAT_range <- calc(ECMWF_SAT, fun = function(x){max(x, na.rm = T) - min(x, na.rm = F)})
SAT_cv    <- calc(ECMWF_SAT, fun = cv)

ECMWF_SST <- brick("/Users/tgagne/Biodiversity2018/data/ECMWF/SST07_17.nc") - 273.15
SST_mean  <- calc(ECMWF_SST, fun = mean )
SST_sd    <- calc(ECMWF_SST, fun = sd )
SST_range <- calc(ECMWF_SST, fun = function(x){max(x, na.rm = T) - min(x, na.rm = F)})
SST_cv    <- calc(ECMWF_SST, fun = cv)

SST_ECMWF <- stack(SST_mean,SST_sd,SST_range,SST_cv) %>% rotate() %>% spatial_sync_raster(reference = sea_mondo)
SAT_ECMWF <- stack(SAT_mean,SAT_sd,SAT_range,SAT_cv) %>% rotate() %>% spatial_sync_raster(reference = sea_mondo)

names(SST_ECMWF) <- c("SST_means","SST_sd","SST_range","SST_cv")
names(SAT_ECMWF) <- c("SAT_means","SAT_sd","SAT_range","SAT_cv")


par(mfrow=(c(4,2)))
plot(SST_ECMWF)

par(mfrow=(c(4,2)))
plot(SAT_ECMWF)

# writeRaster(SST_ECMWF, "/Users/tgagne/Biodiversity2018/data/ECMWF/ECMWF_SST_Jan15.grd", format="raster", overwrite=TRUE)
# writeRaster(SAT_ECMWF, "/Users/tgagne/Biodiversity2018/data/ECMWF/ECMWF_SAT_Jan15.grd", format="raster", overwrite=TRUE)








# function for creating file name with leading zeros
# makes it easier to process them sequentially
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.png', sep=''))
  }
}

for(i in 1:196){
name <- rename(i)
png(name)
ncin <- raster("/Users/tgagne/output.nc", band = i)
ncin
plot(ncin)
dev.off()
}

my_command <- 'convert *.png -delay 3 -loop 0 animation.gif'
system(my_command)










