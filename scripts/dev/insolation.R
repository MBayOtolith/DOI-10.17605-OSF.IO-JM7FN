#################################
#  Solar Insolation stack prep  #
#################################
library(WaveletComp)                 # continous wavelet power analysis 
library(lubridate)                   # date handling
library(reshape2)                    # df manipulation  
library(ggplot2)                     # plotting
library(dplyr)                       # df manipulation
library(raster)
library(viridis)
library(spatial.tools)

# KEEP in mind that within this script there is a bug that I resolve in the Single_File_Script.R later
# The seasonal wavelet intensity calculation generates rasters that are flipped 180 degrees. 
# This is resolved in the final "Single script"

setwd("~/Biodiversity-and-productivity-2017/data/NEO data/solar_insolation")  # set working directory of solar insolation data
list_r <- list.files('.') %>% lapply( FUN = raster )                          # get me the names of the files for solar insolation data, turn each of those files in to a raster
list_noMore <- do.call(what = stack, args = list_r)                           # stack all the layers (30 day composite)

# sea_mondo comes from Single_File_Script.R match to sea_mondo current res (i.e. lowest/highest common)
list_noMore <- spatial_sync_raster(list_noMore, reference = sea_mondo, method  = "bilinear", verbose = T)

# manipulate stack in to DF format that works with wavelet transformation
terre_df_1 <- data.frame( rasterToPoints( list_noMore ) )          # to data frame 
terre_df_1_cc <- terre_df_1[complete.cases(terre_df_1), ]          # remove rows containing NA, listwise deletion #

terre_df_1_cc$x <- NULL
terre_df_1_cc$y <- NULL
emp_frame_6mth <- NULL
emp_frame_12mth <- NULL

#####################
# wavelet intensity #   # only rerun if needed, raster is available for loading below
#####################
for( i in 1:nrow(terre_df_1_cc) ){
  dates_str <- colnames(terre_df_1_cc)
  cell_ts   <- terre_df_1_cc[i,] %>% as.vector() %>% as.numeric()
  bound <- data.frame(dates_str,cell_ts)
  bound$dates_str <- as.character(bound$dates_str)
  test_split <- colsplit(string = bound$dates_str, pattern = c("\\_"), names = c("1B","2B","aa1","aa2","aa3"))
  test_split <- colsplit(string = test_split$aa2, pattern = "\\.", names = c("year", "month", "trash"))
  bound <- data.frame(year = test_split$year, month = test_split$month, value = cell_ts)
  bound$date <- paste0(bound$month,"/",1,"/",bound$year) %>% mdy()
  my.w = analyze.wavelet(bound, "value",
                         loess.span = 0,
                         dt = 1, dj = 1/5,
                         lowerPeriod = 1,
                         upperPeriod = 16,
                         make.pval = T, n.sim = 10,
                         verbose = F)
  #wt.image(my.w, color.key = "quantile", n.levels = 250,legend.params = list(lab = "wavelet power levels", mar = 4.7), show.date = T)
  #reconstruct(my.w, plot.waves = F, lwd = c(1,2), legend.coords = "bottomleft")
  # pull out 6 month and 12 month +- 2months averages of wavelet power levels. This is periodicity intensity
  # i.e. low average power level for the 6 month window = low intensity 6 month periodicity
  test_wave <- as.data.frame(my.w$Power)
  test_wave$period <- my.w$Period
  test_wave <- melt(test_wave, id.var = "period", variable.name = "month")
  test_wave$month <- as.numeric(test_wave$month)
  
  #ggplot(test_wave,aes(x=month,y=period, color = value))+geom_point()+scale_color_distiller(palette = "Spectral")
  six_month_intensity <- as.numeric( test_wave %>% filter( period > 3, period < 9) %>% summarise(mean(value)) )
  emp_frame_6mth <- rbind(six_month_intensity, emp_frame_6mth)
  twelve_month_intensity <- as.numeric( test_wave %>% filter( period > 9, period < 15) %>% summarise(mean(value)) )
  emp_frame_12mth <- rbind(twelve_month_intensity, emp_frame_12mth)
  print(i)
}

# bind cell_id and x-y with emp_frame
six_month_inten <- as.data.frame(emp_frame_6mth)
six_month_inten$x <- terre_df_1[complete.cases(terre_df_1), ]$x
six_month_inten$y <- terre_df_1[complete.cases(terre_df_1), ]$y
twelve_month_inten <- as.data.frame(emp_frame_12mth)
twelve_month_inten$x <- terre_df_1[complete.cases(terre_df_1), ]$x
twelve_month_inten$y <- terre_df_1[complete.cases(terre_df_1), ]$y

ggplot(six_month_inten,aes(x,y,fill = V1))+geom_tile()+scale_fill_distiller(palette = "Spectral")+coord_fixed()
ggplot(twelve_month_inten,aes(x,y,fill = V1))+geom_tile()+scale_fill_distiller(palette = "Spectral")+coord_fixed()

xyz <- data.frame( x = six_month_inten[,"x"], y = six_month_inten[,"y"],twel_inten = twelve_month_inten[,"V1"],six_inten = six_month_inten[,"V1"])
sun_phencyc_intensity <- rasterFromXYZ(xyz, crs = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
plot(sun_phencyc_intensity)
#writeRaster(sun_phencyc_intensity, "sun_phencyc_intensity.grd", format="raster", overwrite=TRUE)




# get conventional inter-intra summary stats #
#############################################
####    INSOLATION ANNUAL SUM STATS      ####   
#############################################

# load the phenocycle intensity stack that was calculated earlier above
sun_phencyc_intensity <- stack("/Users/tgagne/Biodiversity-and-productivity-2017/data/NEO data/solar_grid_files/sun_phencyc_intensity.grd")

years <- seq(2006,2016,by=1) #sequence of years, adjustable final year to speed up or slowdown process
sun_stack <- NULL

for(x in 1:length(years)){
  all_file_names <- list.files('.')
  names_of_year_x <- grep(years[x],all_file_names)
  file_names_of_year_x <- all_file_names[names_of_year_x]
  list_r <- lapply(file_names_of_year_x, FUN = raster )        # turn each filename in to a raster a build a list
  #list_r <- lapply(list_r,                 FUN = na_remove )   # turn 99999 coded cells in to NA
  list_noMoreAN <- do.call(what = stack, args = list_r)          # stack the year
  summed <- calc(list_noMoreAN, fun = sum, na.rm = T)            # calculate the sum of that year
  #summed<- reclassify(summed,c(-Inf,0,NA))                     # reclassify below 0 to NA
  #summed <- reclassify(summed, c(quantile(summed, probs = c(0.9), type = 7), Inf, quantile(summed, probs = c(0.9), type = 7)))
  names(summed) <- years[x]
  sun_stack <- stack(summed,sun_stack)
}
sun_stack        # should be a stack of several years, depends on length of year sequence above

# sun productivity metrics with annual sums
sun_stack_mean <- calc(sun_stack, fun = mean)
sun_stack_sd <- calc(sun_stack, fun = sd)
sun_stack_cv <- calc(sun_stack, fun = cv)
sun_stack_range <- calc(sun_stack, fun = function(x){max(x) - min(x)})              
sun_metrics <- stack(sun_stack_mean,sun_stack_sd,sun_stack_cv,sun_stack_range)
sun_metrics <- spatial_sync_raster(sun_metrics, sea_mondo, method = "bilinear", verbose = T)
names(sun_metrics) <- c("ann_mean_sun","ann_std_dev_sun","ann_cv_sun","ann_range_sun")
plot(sun_metrics)

# sub-annual metrics #
sub_sd         <- calc(list_noMore, fun = sd, na.rm = T)            # calculate the total time period SD
sub_cv         <- calc(list_noMore, fun = cv, na.rm = T)            # calculate the total time period CV
sub_mean       <- calc(list_noMore, fun = mean, na.rm = T)          # calculate the total time period mean
sub_range      <- calc(list_noMore, fun = function(x){max(x, na.rm = T) - min(x, na.rm = T)})
sub_sun_stack <- stack(sub_sd,sub_cv,sub_mean,sub_range)
names(sub_sun_stack) <- c("sub_sd_sun","sub_cv_sun","sub_mean_sun","sub_range_sun")

# stack them all in to a single stack that can be pulled in to domain stack in single_file_script.R #
sun_metrics
sub_sun_stack
sun_phencyc_intensity

sun_stack <- stack(sun_metrics,sub_sun_stack,sun_phencyc_intensity)
#setwd('/Users/tgagne/Biodiversity-and-productivity-2017/data/NEO data/solar_grid_files')
#writeRaster(sun_stack, "sun_stack.grd", format="raster", overwrite=TRUE)
