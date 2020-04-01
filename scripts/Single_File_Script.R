# FULL SCRIPT - October 24, 2017

# In this script I am pulling together a number of disparate scripts that were intially
# developed in an exploratory sense. Because the project is moving forward from an initially EDA
# style of work in to a product driven stage; I feel like parsing of more extemporaneous scripts will
# help clean up what is currently a hard to follow procedure and seqeunce. 


# First lets load relavent libraries and comment on what their relavancy is
library(raster)         # Basic R raster management package
library(tidyverse)      # R packages for data science, highly capable plotting tool, data manipulation, orginization, df, etc
library(reshape2)       # useful tool to cast and melt dataframes
library(viridis)        # scale color 
library(broom)          # turn model objects in to tidy data frames
library(caret)          # ML modeling fitting environment
library(spatial.tools)  # spatial sync tool is excellent for matching projections and resampling resolutions
library(gridExtra)      # useful for displaying multiple ggplot objects in a mfrow(par()) style layout
library(pdp)            # useful machine model visualization
library(ICEbox)         # useful partial dependency variation
library(RColorBrewer)   # palette manipulation
library(spdep)          # auto-covariate for generated area weighted var to address spatial-autocorrelation
library(usdm)           # raster Variogram
library(grid)           # plot static annotation
library(scales)         # single col rescale func
library(ClusterR)       # speed up raster calculation with paralellization
library(rnaturalearth)  # contains useful shapefiles for map overlays
library(ggtern)         # ggplot ternary plot extension
library(MASS) 
library(scales)
library(plyr)   
library(gstat)          # useful geospatial analysis library with variograms and model fittings
library(rgdal)

###############################
 #     FUNCTION LIBRARY      #
###############################
rescale <-                                                             # function for rescaling raster values
  function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))}
na_remove <-                                                           # function to handle 99999 coded NA values
  function(x){      
  x[ x[] == 99999 ] <- NA 
  return(x)}
zero_recode <-                                                         # recode all persistent 0 values to NA
  function(x){          
  x[ x[] == 0 ] <- NA 
  return(x)}
rangisize <-                                                           # 0-1 standardization function for data_frames
  function(x) { 
  x <- sweep(x, 2, apply(x, 2, min)) 
  sweep(x, 2, apply(x, 2, max), "/") }
quantile_trim <-                                                       # quantile trimming function, adjust % in function at 99%/1% currently
  function(df, varname){                                               # trim quantiles function 
  
  high_tail <- 0.99
  low_tail  <- 0.01
  
  upper <- quantile(df[,varname],  c(high_tail)) %>% as.numeric()
  lower <- quantile(df[,varname],  c(low_tail) )  %>% as.numeric()
  
  #df <- transmute(df[,varname] = replace(df[,varname], df[,varname] < lower, lower))
  #df <- transmute(df[,varname] = replace(df[,varname], df[,varname] > upper, upper))
  
  df[,varname] <- ifelse(df[,varname] > lower, df[,varname], lower)
  df[,varname] <- ifelse(df[,varname] < upper, df[,varname], upper)
  #df <-  df[ (lower < df[varname] & df[varname] < upper), ]
  return(df)
  #str(df)
  }

MANY_LOESS <-          
  function(data_l,SPAN,N_size,response,predictor,pt_alpha,N_LOESS, FRAC,ymax, color, pt_color){
  LAND_LOESS_DF <- data.frame(mean = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500))
  
  for(i in 1:N_LOESS){
    # sample 1000 points
    LAND_sample <- data_l %>% sample_n(N_size)
    # fit a loess
    xx <- LAND_sample[,predictor]
    yy <- LAND_sample[,response]
    tp_est <- loess(yy ~ xx , span = SPAN) 
    # predict accross range of x using loess model
    loess_vec <- data.frame(predict(tp_est, newdata = data.frame(xx = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500)) ))
    colnames(loess_vec) <- as.character(i)
    # repeat x times
    LAND_LOESS_DF <- cbind(LAND_LOESS_DF,loess_vec)
  }
  
  LAND_long <- melt(LAND_LOESS_DF, id = "mean")
  point_sample <- sample_n(data_l,FRAC)
  xp <- point_sample[,predictor]
  yp <- point_sample[,response]
  LOW_PLOT <<-ggplot()+
    geom_point(aes(x=xp,y=yp), alpha = pt_alpha ,size = 0.02, color = pt_color) +
    geom_line(data = LAND_long,aes(x=mean,y=value,group=variable),size = 0.1, color = color)+ 
    annotation_custom( grobTree( textGrob(predictor, x = 0.1, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6 ) ) ) )+
    scale_y_continuous(limits = c(0,ymax),expand = c(0.005, 0.005))+
    scale_x_continuous(expand = c(0, 0)) +
    xlab(NULL)+
    ylab(NULL)+
    themeo + 
    theme(axis.text.x = element_blank(),
          axis.text.y=element_blank(),
          plot.margin=unit(c(0,0,0,0), "null"))
  print(LOW_PLOT)
}
MANY_LOESS_MAR <- 
  function(data_l,SPAN,N_size,response,predictor,pt_alpha,N_LOESS, FRAC,ymax, color, pt_color){
    LAND_LOESS_DF <- data.frame(mean = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500))
    
    for(i in 1:N_LOESS){
      # sample 1000 points
      LAND_sample <- data_l %>% sample_n(N_size)
      # fit a loess
      xx <- LAND_sample[,predictor]
      yy <- LAND_sample[,response]
      tp_est <- loess(yy ~ xx , span = SPAN) 
      # predict accross range of x using loess model
      loess_vec <- data.frame(predict(tp_est, newdata = data.frame(xx = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500)) ))
      colnames(loess_vec) <- as.character(i)
      # repeat x times
      LAND_LOESS_DF <- cbind(LAND_LOESS_DF,loess_vec)
    }
    
    LAND_long <- melt(LAND_LOESS_DF, id = "mean")
    point_sample <- dplyr::sample_n(data_l,FRAC)
    xp <- point_sample[,predictor]
    yp <- point_sample[,response]
    LOW_PLOT <<-ggplot()+
      geom_point(aes(x=xp,y=yp), alpha = pt_alpha ,size = 0.02, color = pt_color) +
      geom_line(data = LAND_long,aes(x=mean,y=value,group=variable),size = 0.1, color = color)+ 
      annotation_custom( grobTree( textGrob(predictor, x = 0.1, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6 ) ) ) )+
      scale_y_continuous(limits = c(0,.6),expand = c(0.005, 0.005))+
      scale_x_continuous(expand = c(0, 0)) +
      xlab(NULL)+
      ylab(NULL)+
      themeo + 
      theme(axis.text.x = element_blank(),
                     axis.text.y=element_blank(),
                     plot.margin=unit(c(0,0,0,0), "null"))
    print(LOW_PLOT)
  }
MANY_LOESS_CV <- 
  function(data_l,SPAN,N_size,response,predictor,pt_alpha,N_LOESS, FRAC,ymax, color, pt_color){
  LAND_LOESS_DF <- data.frame(mean = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500))
  
  for(i in 1:N_LOESS){
    # sample 1000 points
    LAND_sample <- data_l %>% sample_n(N_size)
    # fit a loess
    xx <- LAND_sample[,predictor]
    yy <- LAND_sample[,response]
    tp_est <- loess(yy ~ xx , span = SPAN) 
    # predict accross range of x using loess model
    loess_vec <- data.frame(predict(tp_est, newdata = data.frame(xx = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500)) ))
    colnames(loess_vec) <- as.character(i)
    # repeat x times
    LAND_LOESS_DF <- cbind(LAND_LOESS_DF,loess_vec)
  }
  
  LAND_long <- melt(LAND_LOESS_DF, id = "mean")
  
  # gather quantiles from data
  sum_data <- do.call(data.frame, aggregate(value ~ mean, data = LAND_long, 
                                            function(x) c(quantile(x, 0.975), quantile(x, 0.025), median(x))))
  # build plots
  
  point_sample <- dplyr::sample_n(data_l,FRAC)
  xp <- point_sample[,predictor]
  yp <- point_sample[,response]
  LOW_PLOT <<-ggplot()+
    #geom_point(aes(x=xp,y=yp), alpha = pt_alpha ,size = 0.02, color = pt_color) +
    geom_ribbon(data = sum_data, aes(x=mean,ymax=value.97.5.,ymin = value.2.5.),  
                size = 0, 
                #color = color, 
                fill = color, 
                alpha = .5) +
    geom_line(data = sum_data, aes(x=mean,y=value.V3),       size = 1, color = color)+ # Median
    #geom_line(data = sum_data, aes(x=mean,y=value.97.5.), size = 0.3, color = color)+ # Upper
    #geom_line(data = sum_data, aes(x=mean,y=value.2.5.),  size = 0.3, color = color)+ # Lower
    annotation_custom( grobTree( textGrob(predictor, x = 0.1, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6 ) ) ) )+
    scale_x_continuous(expand = c(0.000, 0.000), breaks = c(.25,.5,.75))+
    scale_y_continuous(expand = c(0.000, 0.000), breaks = c(0.2,0.4,0.6)) +
    coord_cartesian(ylim = c(c(0,.75)))+
    xlab(NULL)+
    ylab(NULL)+
    themeo + 
    theme(axis.text.x = element_blank(),
          axis.text.y=element_blank(),
          plot.margin=unit(c(4,4,4,4), "pt"))
  print(LOW_PLOT)
}
MANY_LOESS_CV_MAR <- 
  function(data_l,SPAN,N_size,response,predictor,pt_alpha,N_LOESS, FRAC,ymax, color, pt_color){
  LAND_LOESS_DF <- data.frame(mean = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500))
  
  for(i in 1:N_LOESS){
    # sample 1000 points
    LAND_sample <- data_l %>% sample_n(N_size)
    # fit a loess
    xx <- LAND_sample[,predictor]
    yy <- LAND_sample[,response]
    tp_est <- loess(yy ~ xx , span = SPAN) 
    # predict accross range of x using loess model
    loess_vec <- data.frame(predict(tp_est, newdata = data.frame(xx = seq(min(data_l[,predictor]), max(data_l[,predictor]), length.out = 500)) ))
    colnames(loess_vec) <- as.character(i)
    # repeat x times
    LAND_LOESS_DF <- cbind(LAND_LOESS_DF,loess_vec)
  }
  
  LAND_long <- melt(LAND_LOESS_DF, id = "mean")
  
  # gather quantiles from data
  sum_data <- do.call(data.frame, aggregate(value ~ mean, data = LAND_long, 
                                            function(x) c(quantile(x, 0.975), quantile(x, 0.025), median(x))))
  # build plots
  
  point_sample <- dplyr::sample_n(data_l,FRAC)
  xp <- point_sample[,predictor]
  yp <- point_sample[,response]
  LOW_PLOT <<-ggplot()+
    #geom_point(aes(x=xp,y=yp), alpha = pt_alpha ,size = 0.02, color = pt_color) +
    geom_ribbon(data = sum_data, aes(x=mean,ymax=value.97.5.,ymin = value.2.5.),  
                size = 0, 
               # color = color, 
                fill = color, 
                alpha = .5) +
    geom_line(data = sum_data, aes(x=mean,y=value.V3),       size = 1, color = color)+ # Median
    #geom_line(data = sum_data, aes(x=mean,y=value.97.5.), size = 0.3, color = color)+ # Upper
    #geom_line(data = sum_data, aes(x=mean,y=value.2.5.),  size = 0.3, color = color)+ # Lower
    annotation_custom( grobTree( textGrob(predictor, x = 0.1, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6 ) ) ) )+
    scale_x_continuous(expand = c(0.000, 0.000), breaks = c(.25,.5,.75))+
    scale_y_continuous(expand = c(0.000, 0.000), breaks = c(0.2,0.4,0.6)) +
    coord_cartesian(ylim = c(c(0,0.75)))+
    xlab(NULL)+
    ylab(NULL)+
    themeo + 
    theme(axis.text.x = element_blank(),
          axis.text.y=element_blank(),
          plot.margin=unit(c(4,4,4,4), "pt"))
  print(LOW_PLOT)
}

##############################
###  Popular ggPlot theme  ###
##############################
themeo <-theme_classic()+
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(margin = margin( 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(c(1, 0.2), unit = "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0)
        )

##############################
###   Useful projections   ###
##############################
# PROJECTION that MODIS (FROM NASA NEO) data came in:
mod_proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# PROJECTION that Clinton uses:
eck <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +over" 
# PROJECTION equal area Berhmann, centered on 195 deg
equalArea <- CRS('+proj=cea +lon_0=195 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

#################################
#                               #
# Pulling in BioDiversity maps  #
#                               #
#################################
# Gabriel's occurence data grid details for reference
# class       : RasterLayer 
#  dimensions  : 360, 720, 259200  (nrow, ncol, ncell)
#  resolution  : 0.5, 0.5  (x, y)
#  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#  coord. ref. : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0
#  phytoplankton, microalge, are not in this?

setwd('~/Biodiversity2018/data/marine diversity/gabriel richness/')
gabriel_richness <- read.csv('TYLER_SR_OBSERVED_05D.csv')
names(gabriel_richness) <- c("x","y","z")
marine_richness <- rasterFromXYZ(gabriel_richness, crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')) 
marine_richness <-  rescale(marine_richness,x.min = marine_richness @data @min ,x.max = marine_richness @data @max, new.min = 0, new.max =1)
names(marine_richness)[[1]] <- "marine_richness"
#plot(marine_richness, col = viridis(50))
#hist(marine_richness)

# OLD Marine Richness data, this is the former aquamaps and clinton's richness inputs
# scaled total Marine Richness from BioDiversityMapping.org
#marine_richness <-'~/Biodiversity2018/data/marine diversity/reprojected_all_spp/marine_export.tif' %>% raster()
#names(marine_richness) <- "marine_richness"
# read in a reprojection of AquaMaps richness, rescale as well
#setwd('/Users/tgagne/Biodiversity2018/data/Aquamaps')
#aquamaps_richness <- raster('am_all_rast.tif') %>% projectRaster(crs = mod_proj) 
#aquamaps_richness <- rescale(aquamaps_richness,
#                             x.min = aquamaps_richness @data @min ,x.max = aquamaps_richness @data @max, 
#                             new.min = 0, new.max =1)

#################################
#            ######             #
#    Terrestrial vertebrates    #     This is species richness data from Clinton's work. biodiversitymapping.org
#            ######             #
#################################

# Clinton's old richness data     # earlier prior recalc with new 50km res and projection
# also these appear to be counts
#terrestrial_mammal_richness <-'~/Biodiversity2018/data/terrestrial diversity/Mammals/richness_10km_all_spp_raster.tif' %>% raster()     # Mammals
#terrestrial_birds_richness  <-'~/Biodiversity2018/data/terrestrial diversity/Birds/richness_10km_all_raster.tif' %>% raster()           # Birds
# the amphib layer below is tricky, areas with no amphibs are reported as NA rather than 0, this gets propagated throughout the model
#terrestrial_amphib_richness <-'~/Biodiversity2018/data/terrestrial diversity/Amphibians/richness_10km_all_spp_raster.tif' %>% raster()  # Amhpibians

# Clinton's recalculated richness data - Mar2018
terrestrial_mammal_richness <-'~/Biodiversity2018/data/terrestrial diversity/BioD_mar22/richness_50km_Mammals_Feb20111.tif' %>% raster()     # Mammals
terrestrial_birds_richness  <-'~/Biodiversity2018/data/terrestrial diversity/BioD_mar22/richness_50km_Birds_v6_spp_e11.tif' %>% raster()           # Birds
# the amphib layer below is tricky, areas with no amphibs are reported as NA rather than 0, this gets propagated throughout the model
terrestrial_amphib_richness <-'~/Biodiversity2018/data/terrestrial diversity/BioD_mar22/richness_50km_Amphibians_Dec11.tif' %>% raster()  # Amhpibians

#
NAvalue(terrestrial_mammal_richness) <- 255
NAvalue(terrestrial_birds_richness)  <- 65536  #uncertain why?
NAvalue(terrestrial_amphib_richness) <- 255

# rename raster
names(terrestrial_mammal_richness) <- "mammals"
names(terrestrial_birds_richness) <- "birds"
names(terrestrial_amphib_richness) <- "amphibs"


plot(terrestrial_mammal_richness)
plot(terrestrial_birds_richness)
plot(terrestrial_amphib_richness)

# I still need terrestrial inverts.

######################################
# combine these rasters into a stack #
######################################
# extent is mismtached, addressed below
#xmin <- min(bbox(terrestrial_mammal_richness)[1,1], bbox(terrestrial_birds_richness)[1,1], bbox(terrestrial_amphib_richness)[1,1])
#xmax <- max(bbox(terrestrial_mammal_richness)[1,2], bbox(terrestrial_birds_richness)[1,2], bbox(terrestrial_amphib_richness)[1,2])  
#ymin <- min(bbox(terrestrial_mammal_richness)[2,1], bbox(terrestrial_birds_richness)[2,1], bbox(terrestrial_amphib_richness)[2,1])  
#ymax <- max(bbox(terrestrial_mammal_richness)[2,2], bbox(terrestrial_birds_richness)[2,2], bbox(terrestrial_amphib_richness)[2,2])  
#newextent <- c(xmin, xmax, ymin, ymax)

#terrestrial_mammal_richness = extend(terrestrial_mammal_richness, newextent)
#terrestrial_birds_richness  = extend(terrestrial_birds_richness, newextent)
#terrestrial_amphib_richness = extend(terrestrial_amphib_richness, newextent)
s123 = stack(terrestrial_mammal_richness, terrestrial_birds_richness, terrestrial_amphib_richness)

# calculate sum total of species = single rasterlayer
# & handling the tricky zero/NA amphib locations..
# s123[is.na(s123[])] <- 0              # recode all NA to 0  # this is an issue in earlier files. No 0 is coded is 0 early
spp_sum <- calc(s123, fun = sum)      # sum all the rasters
# spp_sum<-zero_recode(spp_sum)

# rescale rasters 0-1
terrestrial_richness <- rescale(spp_sum,
                                x.min = spp_sum @data @min,
                                x.max = spp_sum @data @max, 
                                new.min = 0, new.max =1)
names(terrestrial_richness) <- "terrestrial_richness"

rm(spp_sum,s123,terrestrial_amphib_richness,terrestrial_birds_richness,terrestrial_mammal_richness,xmin,xmax,ymin,ymax,newextent)

# reproject diversity layers in to common NASA NEO modis output
marine_richness      <- projectRaster(marine_richness, crs = mod_proj)
terrestrial_richness <- projectRaster(terrestrial_richness, crs = mod_proj)

###############################################
###############################################
 #### DIVERSITY LAYER INPUT COMPLETE HERE ####
###############################################
###############################################

##########################################
#                                        #
#    Productivity metrics generation     #
#                                        #
##########################################

###########################################################################
  ####    TERRESTRIAL PRODUCTIVITY ANNUAL SUM METRIC DEVELOPMENT    ####     i.e. yearly stacks for interannual metric dev
###########################################################################
setwd("~/Biodiversity2018/data/NEO data/Terrestrial")
years <- seq(2003,2016,by=1) #sequence of years, adjustable final year to speed up or slowdown process
land_stack <- NULL

for(x in 1:length(years)){
      all_file_names <- list.files('.')
      names_of_year_x <- grep(years[x],all_file_names)
      file_names_of_year_x <- all_file_names[names_of_year_x]
      list_r <- lapply(file_names_of_year_x, FUN = raster )        # turn each filename in to a raster a build a list
      list_r <- lapply(list_r,                 FUN = na_remove )   # turn 99999 coded cells in to NA
      list_noMore <- do.call(what = stack, args = list_r)          # stack the year
      summed <- calc(list_noMore, fun = sum, na.rm = T)            # calculate the sum of that year
      summed<- reclassify(summed,c(-Inf,0,NA))                     # reclassify below 0 to NA
      summed <- reclassify(summed, c(quantile(summed,              # tail values above 99% quantile are truncated at 99% values
                                              probs = c(0.99), type = 7), Inf, quantile(summed, probs = c(0.99), type = 7)))
      names(summed) <- years[x]
      land_stack <- stack(summed,land_stack)
    }
land_stack        # this is a stack of yearly sums of NDVI

##############################################################
     ####    MARINE PRODUCTIVITY ANNUAL SUM DEV      ####          i.e. yearly stacks for interannual metric dev
##############################################################
setwd("~/Biodiversity2018/data/NEO data/Marine")
sea_stack <- NULL

for(x in 1:length(years)){
      all_file_names <- list.files('.')
      names_of_year_x <- grep(years[x],all_file_names)
      file_names_of_year_x <- all_file_names[names_of_year_x]
      list_r <- lapply(file_names_of_year_x, FUN = raster )        # turn each filename in to a raster a build a list
      list_r <- lapply(list_r,                 FUN = na_remove )   # turn 99999 coded cells in to NA
      list_noMore <- do.call(what = stack, args = list_r)          # stack the year
      summed <- calc(list_noMore, fun = sum, na.rm = T)            # calculate the sum of that year
      summed<- reclassify(summed,c(-Inf,0,NA))                     # reclassify below 0 to NA
      summed <- reclassify(summed, c(quantile(summed,              # tail values above 99% quantile are truncated at 99% values
                                              probs = c(0.99), type = 7), Inf, quantile(summed, probs = c(0.99), type = 7)))
      names(summed) <- years[x]
      sea_stack <- stack(summed,sea_stack)
}

rm(all_file_names,file_names_of_year_x,list_noMore,list_r,summed,x,years,names_of_year_x)
sea_stack       # this is a stack of yearly sums of Chl-A


###################################################################
##     Developing ANNUAL metrics of productivity and stacking    ##
###################################################################
####  This is where we should consider quantile trimming before metric calculation...
####  Or not, considering those cells just become NA anyway and are not in the model regardless of where they are lost?

# Ocean productivity metrics with annual sums
    sea_stack_mean <- calc(sea_stack, fun = mean)            
    sea_stack_sd <- calc(sea_stack, fun = sd)
    sea_stack_cv <- calc(sea_stack, fun = cv)
sea_stack_range <- calc(sea_stack, fun = function(x){max(x) - min(x)})               
    sea_metrics <- stack(sea_stack_mean,sea_stack_sd,sea_stack_cv,sea_stack_range)
    names(sea_metrics) <- c("ann_mean_pp","ann_std_dev_pp","ann_cv_pp","ann_range_pp")

# Land productivity metrics with annual sums
    land_stack_mean <- calc(land_stack, fun = mean)
    land_stack_sd <- calc(land_stack, fun = sd)
    land_stack_cv <- calc(land_stack, fun = cv)
land_stack_range <- calc(land_stack, fun = function(x){max(x) - min(x)})             
    land_metrics <- stack(land_stack_mean,land_stack_sd,land_stack_cv,land_stack_range)
    names(land_metrics) <- c("ann_mean_pp","ann_std_dev_pp","ann_cv_pp","ann_range_pp")

# clean up
rm(sea_stack_mean,sea_stack_sd,sea_stack_cv,land_stack_mean,land_stack_sd,land_stack_cv,land_stack_range,sea_stack_range)

##########################################
#                                        #
# generate within in year metrics for PP #                                # this the dev of intra-annual summary metrics
#                                        #
########################################## 

# GENERATE SUB-ANNUAL METRICS OF VARIABILITY FOR PP
# TERRESTRIAL
setwd("~/Biodiversity2018/data/NEO data/Terrestrial")
      all_file_names <- list.files('.')
      list_r <- lapply(all_file_names, FUN = raster )                     # turn each filename in to a raster a build a list
      list_r <- lapply(list_r,                 FUN = na_remove )          # turn 99999 coded cells in to NA
      list_noMore <- do.call(what = stack, args = list_r)                 # stack all the layers (14 day composite? double check)
      list_noMore <- list_noMore + .1                                     # NOTE: because NDVI is measured on -0.1 to .9 scale it handles CV with NULL values, 
                                                                          # I am bumping the whole set up .1 units. This didn't show up earlier because the multi-month 
                                                                          # sums ending up generating min values greater than 0 
     # quant <- quantile(list_noMore, probs = c(0.999), type = 7)
     # list_noMore <- reclassify(list_noMore, c(quant, Inf, quant))
      
      sub_sd         <- calc(list_noMore, fun = sd, na.rm = T)            # calculate the total time period SD
      sub_cv         <- calc(list_noMore, fun = cv, na.rm = T)            # calculate the total time period CV
      sub_mean       <- calc(list_noMore, fun = mean, na.rm = T)          # calculate the total time period mean
      sub_range      <- calc(list_noMore, fun = function(x){max(x, na.rm = T) - min(x, na.rm = T)})
   
      sub_land_stack <- stack(sub_sd,sub_cv,sub_mean,sub_range)
      names(sub_land_stack) <- c("sub_sd_pp","sub_cv_pp","sub_mean_pp","sub_range_pp")

# MARINE
setwd("~/Biodiversity2018/data/NEO data/Marine")
      all_file_names <- list.files('.')
      list_r <- lapply(all_file_names, FUN = raster )                     # turn each filename in to a raster a build a list
      list_r <- lapply(list_r,                 FUN = na_remove )          # turn 99999 coded cells in to NA
      list_noMore <- do.call(what = stack, args = list_r)                 # stack all the layers (14 day composite? double check)
      
      # quant <- quantile(list_noMore, probs = c(0.999), type = 7)
      #list_noMore <- reclassify(list_noMore, c(quant, Inf, quant))
      
      sub_sd        <- calc(list_noMore, fun = sd, na.rm = T)             # calculate the total time period SD
      sub_cv        <- calc(list_noMore, fun = cv, na.rm = T)             # calculate the total time period CV
      sub_mean      <- calc(list_noMore, fun = mean, na.rm = T)           # calculate the total time period mean
      sub_range     <- calc(list_noMore, fun = function(x){max(x, na.rm = T) - min(x, na.rm = T)})
      
      sub_sea_stack <- stack(sub_sd,sub_cv,sub_mean,sub_range)
      names(sub_sea_stack) <- c("sub_sd_pp","sub_cv_pp","sub_mean_pp","sub_range_pp")

# environment clean-up
rm(list_r,sub_cv,sub_mean, sub_range, list_noMore, all_file_names, sub_sd)

#####################################
##                                 ##
###                               ###
#####                            ####
###### What have we got so far ######
#####                            ####
###                               ###
##                                 ##
#####################################
# Diversity
#aquamaps_richness
marine_richness
terrestrial_richness

# Annual Metrics
sea_metrics
land_metrics

# sub_Annual Metrics
sub_sea_stack
sub_land_stack

##################
# RASTER SYNCING #
##################
# TERRESTRIAL
land_annual_metrics_synced  <- spatial_sync_raster( land_metrics,      terrestrial_richness, method = "ngb", size_only = FALSE, verbose = T)
land_sub_metrics_synced     <- spatial_sync_raster( sub_land_stack,    terrestrial_richness, method = "ngb", size_only = FALSE, verbose = T)

# MARINE
sea_annual_metrics_synced  <- spatial_sync_raster( sea_metrics,       marine_richness, method = "ngb", size_only = FALSE, verbose = T)
sea_sub_metrics_synced     <- spatial_sync_raster( sub_sea_stack,     marine_richness, method = "ngb", size_only = FALSE, verbose = T)
#sea_aquamaps_synced        <- spatial_sync_raster( aquamaps_richness, marine_richness, method = "ngb", size_only = FALSE, verbose = T)

###################
# RASTER STACKING #
###################
land_mondo <- stack(terrestrial_richness
                    ,land_annual_metrics_synced
                    ,land_sub_metrics_synced) 

sea_mondo  <- stack(marine_richness
                    #, sea_aquamaps_synced
                    ,sea_annual_metrics_synced
                    ,sea_sub_metrics_synced)

# more clean-up
rm(land_metrics,sea_stack,land_stack,
   sub_land_stack,sub_sea_stack,sea_metrics,
   aquamaps_richness,terrestrial_richness,
   land_annual_metrics_synced,land_sub_metrics_synced,
   marine_richness, sea_aquamaps_synced,
   sea_annual_metrics_synced,sea_sub_metrics_synced)


#############################################################################
#############################################################################
#############################################################################
###  DONT RUN THIS AGAIN UNLESS NEEDED TO READJUST DATABASE DEVELOPMENT   ###
###  MOVING FORWARD IS DATA VISUALIZATION/EXPLORATION/MODELING            ###
#############################################################################
#############################################################################
#############################################################################
# setwd("/Users/tgagne/Biodiversity2018/data/Rdata/")
# writeRaster(land_mondo, "land_mondo_Mar22.grd", format="raster", overwrite=TRUE)
# writeRaster(sea_mondo,  "sea_mondo_Mar22.grd" , format="raster", overwrite=TRUE)
# land_mondo = stack( "land_mondo_Mar22.grd" )
# sea_mondo  = stack( "sea_mondo_Mar22.grd" )

# Abiotic inputs should be synced and stacked here to the mondo stacks #
# options to look at:
# elevation (bathy/topography) - We have this in GEBCO and SRTM
# net solar radiation          - We can get this in NASA CERES, net radiation, diff between how much enters and how much escapes
# solar insolation             - CEREs also, differs in that insolulation = how much reaches surface f(cloud, incident angle, atmosphere, aerosols, etc.)
# SST                          - average, range, 12 months of the year available, 2003 - 2016, via ERRSST

# see merge further down
SAT <- stack( "~/Biodiversity2018/data/ECMWF/ECMWF_SAT_Jan15.grd")
SST <- stack( "~/Biodiversity2018/data/ECMWF/ECMWF_SST_Jan15.grd")


   #############
  ##           ##
#### elevation ####
  ##           ##
   #############
bathy <- raster('~/Biodiversity2018/data/NEO data/Topo/GEBCO_BATHY_2002.FLOAT.TIFF') %>% 
  na_remove() %>% spatial_sync_raster(sea_mondo, method = "ngb", size_only = FALSE, verbose = T)
names(bathy) <- "elev"
topo  <- raster('~/Biodiversity2018/data/NEO data/Topo/SRTM_RAMP2_TOPO_2000-02-11_rgb_3600x1800.FLOAT.TIFF') %>% 
  na_remove() %>% spatial_sync_raster(land_mondo, method = "ngb", size_only = FALSE, verbose = T)
names(topo) <- "elev"

### Aspect and Slope calc ###
### Neighborhood = 8 is computation of slope and aspect are computed according to Horn 1981
#slope_marine   <- terrain(x = bathy, opt = 'slope', unit = 'degrees', neighbors = 8)
#aspect_marine  <- terrain(x = bathy, opt = 'aspect', unit = 'degrees', neighbors = 8)
#slope_land     <- terrain(x = topo, opt = 'slope', unit = 'degrees', neighbors = 8)
#aspect_land    <- terrain(x = topo, opt = 'aspect', unit = 'degrees', neighbors = 8)

### ### ### ### ### ### ####
# solar insolation summary # # cv,sd,mean,12 inten, etc.
### ### ### ### ### ### ####
sun_stack <- stack('~/Biodiversity2018/data/NEO data/solar_grid_files/sun_stack.grd') %>% spatial_sync_raster(sea_mondo, method = "bilinear", verbose = T)

# The solar insolation intensity approximations coordinate system 
# are rotated incorrectly 180 degrees. This is resolved below:
sun_inten_fix <- data.frame(rasterToPoints(sun_stack))[,c("x","y","six_inten","twel_inten")]
sun_inten_fix[,c("x","y")] <- sun_inten_fix[,c("x","y")] * -1 # rotate the whole plane 180 degrees by * -1
sun_inten_fix <- rasterFromXYZ(sun_inten_fix, crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% spatial_sync_raster(sea_mondo, method = "bilinear", verbose = T)
sun_stack <- dropLayer(sun_stack,c("six_inten","twel_inten"))
sun_stack <- stack(sun_stack, sun_inten_fix)

sun_stack_land <- spatial_sync_raster(sun_stack,land_mondo, verbose =T)


###################################################################
###       THIS IS DATA VISUALIZATION/EXPLORATION/MODELING       ###
###################################################################
# stack in the additional inputs #
land_mondo <- stack(land_mondo, topo,  sun_stack_land) #, SAT )#SAT_means, SAT_sd, SAT_range)
sea_mondo  <- stack(sea_mondo,  bathy, sun_stack) #,      SST)# SST_means,SST_sd,SST_range,)
###################################################
### DOWNSAMPLE HERE for working the steps below ###
###         Also match resolutions              ###
###################################################
# match land_mondo to sea_mondo res
# downsample below: match land_mondo to sea_mondo res
land_mondo <- raster::resample(x = land_mondo,y = sea_mondo, method = "bilinear" )

#stack in new ECMWF SST and SAT
land_mondo <- stack(land_mondo, SAT)
sea_mondo  <- stack(sea_mondo, SST)

###################################################
       # Second round of finalized script #
###################################################

setwd("~/Biodiversity2018/data/Rdata/")
# writeRaster(land_mondo, "land_mondo_post_Mar22.grd", format="raster", overwrite=TRUE)
# writeRaster(sea_mondo,  "sea_mondo_post_Mar22.grd" , format="raster", overwrite=TRUE)
land_mondo = stack( "land_mondo_post_Mar22.grd" )
sea_mondo  = stack( "sea_mondo_post_Mar22.grd" )

# richness without inverts for ancillary analysis post-review

gabriel_richness <- read.csv('~/Biodiversity2018/data/marine diversity/gabriel richness/Marine_Vertebrate_SR_NBOBS_05D_GRID.csv')[,1:3] %>% 
  rasterFromXYZ(crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') )

gabriel_richness <- zero_recode(gabriel_richness)
gabriel_richness <- spatial_sync_raster(unsynced = gabriel_richness, reference = sea_mondo)
gabriel_richness <- reclassify(gabriel_richness, c(quantile(gabriel_richness,probs = c(0.9999), type = 7), Inf, quantile(gabriel_richness, probs = c(0.9999), type = 7)))
gabriel_richness <-  rescale(gabriel_richness,x.min = gabriel_richness @data @min ,x.max = gabriel_richness @data @max, new.min = 0, new.max =1)

par(mfrow = c(2,2))
hist(gabriel_richness)
plot(gabriel_richness)
hist(sea_mondo$marine_richness)
plot(sea_mondo$marine_richness)

sea_mondo$marine_richness <- gabriel_richness
plot(sea_mondo$marine_richness)

## Adding a few additional covariates
#  ocean o2 concentration
setwd('~/Biodiversity2018/data/Gabriel/')
gabriel_OXY <- read.csv('OXYGEN_DATA.csv')
names(gabriel_OXY)[1:2] <- c("x","y")
gabriel_OXY <- rasterFromXYZ(gabriel_OXY, crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')) 

# terrestrial precipitation
gabriel_PRECIP<- read.csv('Precip_DATA.csv')
names(gabriel_PRECIP)[1:2] <- c("x","y")
gabriel_PRECIP <- rasterFromXYZ(gabriel_PRECIP, crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')) 

# NDVI SWap Analysis - seasonal wavelet intensity
gabriel_NDVI_swap<- read.csv('NDVI_SWAP_Jan30.csv')
names(gabriel_NDVI_swap)[1:2] <- c("x","y")
gabriel_NDVI_swap <- rasterFromXYZ(gabriel_NDVI_swap, crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')) %>% spatial_sync_raster(reference = land_mondo)

# ChlA SWap Analysis - seasonal wavelet intensity
gabriel_ChlA_swap<- read.csv('SWAP_CHLA_Mar22.csv')
names(gabriel_ChlA_swap)[1:2] <- c("x","y")
gabriel_ChlA_swap <- rasterFromXYZ(gabriel_ChlA_swap, crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')) %>% spatial_sync_raster(reference = land_mondo)

# stack them
sea_mondo <- stack(sea_mondo, gabriel_OXY, gabriel_ChlA_swap)
land_mondo<- stack(land_mondo, gabriel_PRECIP, gabriel_NDVI_swap)

# clean up
rm(gabriel_OXY, gabriel_PRECIP,gabriel_NDVI_swap,gabriel_ChlA_swap,SAT_cv, SST_cv,eck, zero_recode)

#### #### ##### ##### #### #### ##### #### ##### #### ##### #### #### #### ##
###   RESAMPLE rasters down to lower resolutions here for additional runs   ###
##### #### ##### ##### #### #### ##### ######### #### ##### #### #### #### ##
# starting point is 0.5 degree
#land_mondo <- aggregate(land_mondo, fact = 10 ) # .5 deg start: fact 2: 2 deg, fact 6: 3 deg, fact 10: 5 deg
#sea_mondo  <- aggregate(sea_mondo,  fact = 10 ) #

####? ####? ####? ####? ####? ####? ####? ####? ####
##### Reproject in to Behrmann equal area here? #### non-finite point loss is a function of edge and corner point loss
####? ####? ####? ####? ####? ####? ####? ####? ####
land_mondo <- projectRaster(land_mondo, crs = equalArea)
sea_mondo  <- projectRaster(sea_mondo,  crs = equalArea)

##############################################################################################
## Generate a basic plot of richness  and various interesting inputs ## i.e. Fig 1    ########
##############################################################################################
SEA_map <- data.frame( rasterToPoints( sea_mondo ) )   # to data frame #
LAND_map <- data.frame( rasterToPoints( land_mondo ) )

#bring in landmass shapefiles
# generate a reference geotiff from which to use in ArcMap reprojection as reference proj
#setwd('/Users/tgagne/Desktop')
#writeRaster(sea_mondo[[1]], "ref_rast", format = "GTiff")

# Equal Area projection with Pacific Ocean centered meridian
countries_df <- readOGR(dsn = "~/Biodiversity2018/data/land shapefiles/country_repro1", layer = "ocean_50m") %>% tidy() 

# latlon projection convential from Natural Earth Data
#countries_df <- readOGR(dsn = "/Users/tgagne/Biodiversity2018/data/land shapefiles/ne_50m_ocean", layer = "ne_50m_ocean") %>% tidy() 

###############################
# convert stacks to dataframe #
###############################
SEA <- data.frame( rasterToPoints( sea_mondo ) )   # to data frame #
LAND <- data.frame( rasterToPoints( land_mondo ) )

#SEA$elev <- SEA$elev * -1 # decreasing depth = increasing depth
 
SEA_noNA <- SEA[complete.cases(SEA), ]             # remove rows of mismatched NA, listwise deletion #
LAND_noNA <- LAND[complete.cases(LAND), ]

rm( LAND_map, SEA,LAND, SAT_means, 
   SAT_range, SAT_sd, SST_means, SST_range, SST_sd, bathy, topo)

## Window subsets for speedy model testing
#SEA_noNA <- SEA_noNA %>% filter(x > -130,
#                                x < -120,
#                                y >  25,
#                                y <  45)
#LAND_noNA <- LAND_noNA %>% filter(x > -130,
#                                  x < -120,
#                                  y >  25,
#                                  y <  45)


# Quick spatial coverage plot check
# where is the missing vals coming from, means of antartica and greenland
#melt(SEA_noNA, id.vars = c("x","y"))  %>% ggplot(aes(x = x, y = y, fill = value))+facet_wrap(~variable, scales = "free")+geom_raster(show.legend = F)+coord_fixed()
#melt(LAND_noNA, id.vars = c("x","y")) %>% ggplot(aes(x = x, y = y, fill = value))+facet_wrap(~variable, scales = "free")+geom_raster(show.legend = F)+coord_fixed()

############################
# species richness overlay #  Fig. 1
############################
dev.new()
richness<-ggplot()+
  geom_raster(data = SEA_map, aes(x,y, fill = log(marine_richness+1))) +
  geom_raster(data = LAND_noNA, aes(x,y, fill = log(terrestrial_richness+1)) ) +
  geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .25)+
  themeo+
  
  scale_x_continuous(expand = c(-0.005,-0.005)) +
  scale_y_continuous(expand = c(-0.01,-0.01)) +
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent") +
  #theme(axis.text = element_blank(),
  #    axis.title = element_blank())+
  ggtitle("A. global terrestrial and marine richness")+
  coord_fixed()
richness
dev.off()

########################
########################
#   Quantile trimming  #
########################
########################
# generate vector of variable names
# note these indexes, adjust as needed when including new pred/response vars
Land_varnames <- colnames(LAND_noNA)[4:length(LAND_noNA)] 
sea_varnames  <- colnames(SEA_noNA)[4:length(SEA_noNA)]

##############
##  Marine  ##
##############
prePlot <- melt(SEA_noNA, id.vars = c("x","y"))
ggplot(prePlot,aes(x = value, fill = variable))+   # Pre trim histogram plotting check MARINE
  geom_histogram(show.legend = F)+
  facet_wrap(~variable, scales = "free")+
  themeo+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(length(SEA_noNA)))

# loop through varnames # SEA loop
submission <- SEA_noNA 

# loop
for( f in 1:length(sea_varnames)) {
  submission <- quantile_trim(
    df = submission, 
    varname = sea_varnames[f])
}

postPlot <- melt(submission, id.vars = c("x","y"))
ggplot(postPlot,aes(x = value, fill = variable))+    # Post trim histogram plot
  geom_histogram(show.legend = F)+
  facet_wrap(~variable, scales = "free")+
  themeo+
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8, "Dark2"))(length(SEA_noNA))))
SEA_submission <- submission


###################
##  Terrestrial  ##
###################
prePlot <- melt(LAND_noNA, id.vars = c("x","y"))
ggplot(prePlot,aes(x = value, fill = variable))+   # Pre trim histogram plotting check MARINE
  geom_histogram(show.legend = F)+
  facet_wrap(~variable, scales = "free")+
  themeo+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(length(LAND_noNA)))

# loop through varnames # LAND loop
submission <- LAND_noNA 

# loop
for( f in 1:length(Land_varnames)) {
  submission <- quantile_trim(
    df = submission, 
    varname = Land_varnames[f])
}

postPlot <- melt(submission, id.vars = c("x","y"))
ggplot(postPlot,aes(x = value, fill = variable))+    # Post trim histogram plot
  geom_histogram(show.legend = F)+
  facet_wrap(~variable, scales = "free")+
  themeo+
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8, "Dark2"))(length(SEA_noNA))))
LAND_submission <- submission

rm(submission, postPlot, prePlot,LAND_noNA, SEA_noNA)

############################
#  Column standardization  #  # 0-1 for all input metrics
############################
LAND_submission[,Land_varnames] <- rangisize(LAND_submission[,Land_varnames])
SEA_submission[,sea_varnames]   <- rangisize(SEA_submission[,sea_varnames])

# Lets look at what we've got here now #
str(LAND_submission)
str(SEA_submission)
summary(LAND_submission)
summary(SEA_submission)

####################################################
####### Spatial Autocorrelation Assessment ######### There is obvious SAC in the response variable, though we
#################################################### only are concerned about SAC in the residuals of the model see the check of this later

# variogram SEA
#SEA_points <- SEA_submission
#coordinates(SEA_points) = ~ x + y
#EAVariogram=variogram(marine_richness~1, data=SEA_points)
#f_sea <- fit.variogram(SEAVariogram,vgm(c("Exp", "Mat", "Sph")))
#f_sea
#plot(SEAVariogram, model=f_sea) 

# variogram LAND
#LAND_points <- LAND_submission
#coordinates(LAND_points) = ~ x + y
#LANDVariogram=variogram(terrestrial_richness~1, data=LAND_points)
#f_land <- fit.variogram(LANDVariogram,vgm(c("Exp", "Mat", "Sph")))
#f_land
#plot(LANDVariogram, model=f_land) 


# redo varnames after developing covariate
Land_varnames <- colnames(LAND_submission)[4:length(LAND_submission)] # note these indexes, adjust as needed when including new pred/response vars
sea_varnames  <- colnames(SEA_submission)[4:length(SEA_submission)]


#####################
#  Exploratory Viz  #
#####################

#############################################
#       MANY_LOESS function parameters      #
#############################################

#         data      =  ex: LAND_noNA
#         SPAN      =  ex: 1                                                  # LOESS Span parameter
#         N_size    =  ex: 100 , number of points to take to fit loess to     # Size of sample from which to make a single LOESS fit
#         response  =  ex: terrestrial richness 
#         predictor =  ex: mean 
#         pt_alpha  =  ex: 0.2 transparency of background points
#         N_LOESS   =  ex: 200, number of LOESS lines to draw                 # Number of LOESS models to show
#         FRAC      =  ex: number of points to plot
#         ymax      =  ex: y-scale cutoff                                     # y scale max break
#         color     =  line color
#         pt_color  =  point color 

########################################
# Multiple LOOPS for full grid product #  # MANY_LOESS
########################################

#  PARAMS
#SPAN = .95
#N_size = 10000 # SOLID RUNNER, THOUGH SLOW
#N_LOESS = 100
#FRAC = 10000

SPAN = .95
N_size = 500 # QUICK RUNNER, THOUGH SPOTTY
N_LOESS = 100
FRAC = 0
#N_size = 300 # 

par(mfrow = c(2,2))

# transformation ideas, exploration
hist(SEA_submission[,"marine_richness"])
hist(log(SEA_submission[,"marine_richness"]+1))   # this looks like the go to...

hist(LAND_submission[,"terrestrial_richness"])
hist(log(LAND_submission[,"terrestrial_richness"]+1))   # this looks like the go to...

## Log SEA_submission input/output
SEA_submission[,"marine_richness"] <- log(SEA_submission[,"marine_richness"]+1)
LAND_submission[,"terrestrial_richness"] <- log(LAND_submission[,"terrestrial_richness"]+1)

#SEA_submission[,sea_varnames] <- log(SEA_submission[,sea_varnames]+1)
#LAND_submission[,Land_varnames] <- log(LAND_submission[,Land_varnames]+1)

## Rescale back to 0-1
SEA_submission[,"marine_richness"] <- scales::rescale(x = SEA_submission[,"marine_richness"], to = c(0,1))
SEA_submission[,sea_varnames] <- rangisize(SEA_submission[,sea_varnames])
LAND_submission[,"terrestrial_richness"] <- scales::rescale(x = LAND_submission[,"terrestrial_richness"], to = c(0,1))
LAND_submission[,Land_varnames] <- rangisize(LAND_submission[,Land_varnames])

# check ranges
summary(SEA_submission)
summary(LAND_submission)


# LOESS CONFIDENCE INTERVAL RUNS #
### LAND LOOP
pred_names = Land_varnames
land_grp = list()
for( p in 1:length(pred_names)){
  
  MANY_LOESS_CV(data_l = LAND_submission,
             response = "terrestrial_richness",
             predictor = pred_names[p],
             
             SPAN = SPAN,
             N_size = N_size,
             N_LOESS = N_LOESS,
             FRAC = FRAC,
             
             color = "#1a9850", ymax = 1,
             pt_color = "gray",pt_alpha = 1)
  
  land_grp[[p]]  <- LOW_PLOT
}
### Aquamaps model / Marine richness
pred_names = sea_varnames
ocean_AM_grp = list()
for( p in 1:length(pred_names)){
  
  MANY_LOESS_CV_MAR(
    data_l = SEA_submission,
    response = "marine_richness",
    predictor = pred_names[p],
    
    SPAN = SPAN,
    N_size = N_size,
    N_LOESS = N_LOESS,
    FRAC = FRAC,
    
    color = "#4575b4", ymax = 1,
    pt_color = "gray",pt_alpha = 1)
  
  ocean_AM_grp[[p]]  <- LOW_PLOT
}

# TERNARY PLOTS 
# These plots are useful but we need to be careful about how to use them..
# To me they represent mixture composition proxies
# i.e. what is the characterisitc mixes tied to particular richness values..
# we also need to decide about presenting raw points vs an interpolated surface
# simple point plots (high memory with high point count, sample_n() helps)
#####

# hex binning 2d for examples 
#ggplot(sample_n(LAND_submission,10000),aes(x = ann_mean_pp, y = SAT_means))+
#  stat_summary_hex(aes(z = terrestrial_richness),fun = 'median', bins = 10)+
#  geom_point(aes(color = terrestrial_richness))+
#  scale_fill_distiller(palette = "Spectral")+
#  scale_color_distiller(palette = "Spectral")

#ggplot(sample_n(SEA_submission,10000),aes(x = ann_mean_pp, y = SST_means))+
#  stat_summary_hex(aes(z = marine_richness), bins = 10, fun = 'median')+
#  geom_point(aes(color = marine_richness))+
#  scale_color_distiller(palette = "Spectral")+
#  scale_fill_distiller(palette = "Spectral")

################################
# BINNED MEDIANS TERNARY PLOTS #
################################
BINS <- 20

#### FIRST SET - contains SUN, PP, SAT/SST ####
ocean_tern <- ggtern(SEA_submission, aes(x=ann_mean_sun, y =ann_mean_pp , z = SST_means)) +
 geom_tri_tern(bins = BINS,fun = median, aes( value = marine_richness, fill = ..stat..), show.legend = F) +
  theme_arrownormal()+
  scale_fill_distiller(palette = "Spectral")
b<-ggtern(SEA_submission,aes(x=ann_mean_sun, y =ann_mean_pp , z = SST_means, color = marine_richness))+
  geom_point(alpha = .2)+
  theme_arrownormal()+
  scale_color_distiller(palette = "Spectral")
grid.arrange(ocean_tern,b)


land_tern <- ggtern(LAND_submission, aes(x=ann_mean_sun, y =ann_mean_pp , z = SAT_means)) +
  geom_tri_tern(bins = BINS,fun = median, aes( value = terrestrial_richness, fill = ..stat..), show.legend = F) +
  theme_arrownormal()+
  scale_fill_distiller(palette = "Spectral")
b<- ggtern(LAND_submission,aes(x=ann_mean_sun, y = ann_mean_pp , z = SAT_means, color = terrestrial_richness))+
  geom_point(alpha = .2)+scale_color_distiller(palette = "Spectral")+
  theme_arrownormal()
grid.arrange(land_tern,b)


#### SECOND SET - contains XXX #### 
b<-ggtern(SEA_submission,aes(y = sub_cv_pp ,x = sub_cv_sun,  z = SST_cv, color = marine_richness))+
  geom_point(alpha = .2, show.legend = F)+
  theme_arrownormal()+
  scale_color_distiller(palette = "Spectral")
ocean_tern2<-ggtern(SEA_submission, aes(y = sub_cv_pp ,x = sub_cv_sun,  z = SST_cv)) +
  geom_tri_tern(bins = BINS,fun = median, aes( value = marine_richness, fill = ..stat..), show.legend = F) +
  theme_arrownormal()+
  scale_fill_distiller(palette = "Spectral")
grid.arrange(ocean_tern2,b)


b<- ggtern(LAND_submission,aes(y = sub_cv_pp ,x = sub_cv_sun, z = SAT_cv, color = terrestrial_richness))+
  geom_point(alpha = .2, show.legend = F)+
  scale_color_distiller(palette = "Spectral")+
  theme_arrownormal()
land_tern2 <-ggtern(LAND_submission, aes(y = sub_cv_pp ,x = sub_cv_sun,  z = SAT_cv)) +
  geom_tri_tern(bins = BINS,fun = median, aes( value = terrestrial_richness, fill = ..stat..), show.legend = F) +
  theme_arrownormal()+
  scale_fill_distiller(palette = "Spectral")
grid.arrange(land_tern2,b)

##### Grid arrangement of LOESS and ternary plots ####
Fig2layout <- 
  # loess fits
  rbind(
    # LOESS PLOTS
    c(1, 2, 31, 32),
    c(3, 4, 33, 34),
    c(5, 6, 35, 36),
    c(7, 8, 37, 38),
    c(9, 10,39, 40),
    c(11,12,41, 42),
    c(13,14,43, 44),
    c(15,16,45, 46),
    c(17,18,47, 48),
    c(19,20,49, 50),
    c(21,22,51, 52),
    c(23,24,53, 54),
    c(25,26,55, 56),
    c(27,28,57, 58),
    c(29,30,59, 60),
    
    # ternary plots
    c(61,61,62, 62),
    c(61,61,62, 62),
    c(61,61,62, 62),
    c(61,61,62, 62),
    c(63,63,64, 64),
    c(63,63,64, 64),
    c(63,63,64, 64),
    c(63,63,64, 64))
  
# megagrid list
mega_list <- c(land_grp,ocean_AM_grp) 
mega_list[[61]] <- land_tern
mega_list[[62]] <- ocean_tern
mega_list[[63]] <- land_tern2
mega_list[[64]] <- ocean_tern2

# generate grid arrange device plot
# dev.new()
# pdf('/Users/tgagne/Desktop/testFig2.pdf', width = 3, height = 18)
# grid.arrange(grobs = mega_list, layout_matrix = Fig2layout)
# dev.off()

rm(f, FRAC, land_grp, LOW_PLOT, N_LOESS, p, SPAN, ocean_AM_grp, N_size)

grid.arrange(grobs = mega_list[1:60])
grid.arrange(grobs = mega_list[c(61,62,63,64)])


########################################
########################################
# BASIC MODEL DEVELOPMENT AND TESTING  #
########################################
########################################

# were not interested in keeping the combined SWap 1Y/6M amplitude as its already parsed seperately
LAND_submission$AMP_1Y6M <- NULL
SEA_submission$AMP_1Y6M  <- NULL

# large map set for supplement
test_mega_df <- merge(x = LAND_submission, y = SEA_submission, all.x = T, all.y = T)
str(test_mega_df)
mappo <- melt(test_mega_df, id.vars = c("x","y"))
ggplot(mappo,aes(x = x, y = y, fill = value))+
  facet_wrap(~variable, scales = "free")+
  geom_raster(show.legend = F)+
  #coord_fixed()+
  themeo+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_distiller(palette = "Spectral")+
  scale_x_continuous(expand = c(-0.005,-0.005)) +
  scale_y_continuous(expand = c(-0.01,-0.01))

# redo varnames after dropping covariate combined SWap
Land_varnames <- colnames(LAND_submission)[4:length(LAND_submission)] # note these indexes, adjust as needed when including new pred/response vars
sea_varnames  <- colnames(SEA_submission)[4:length(SEA_submission)]

# developing basic model on small subset of full dataset and predicting output

# Look quickly at the datasets
str(LAND_submission)
str(SEA_submission)

# sample to speed workup time
LAND_noNA  <- LAND_submission #%>% sample_n(10000, replace = T)
OCEAN_noNA <- SEA_submission  #%>% sample_n(10000, replace = T)

##
#####
# At this point you can run the sensitivity analysis for variable importance figure.
#####
##

# split in to test and training sets
# create data partition for training set 
LandInTrain  <- createDataPartition(LAND_noNA$terrestrial_richness, p=.95,list=FALSE)
OceanInTrain <- createDataPartition(OCEAN_noNA$marine_richness,     p=.95,list=FALSE)

#manip
LandInTrain  <- createDataPartition(LAND_noNA$terrestrial_richness, p=.90,list=FALSE)
OceanInTrain <- createDataPartition(OCEAN_noNA$marine_richness,     p=.90,list=FALSE)

# partion out training set
LAND_trainer<-LAND_noNA[LandInTrain,]
OCEAN_trainer<-OCEAN_noNA[OceanInTrain,]

# partition out test set
LAND_test<-LAND_noNA[-LandInTrain,]
OCEAN_test<-OCEAN_noNA[-OceanInTrain,]

rm(LandInTrain,OceanInTrain)


# caret foundation setup
numFolds <- trainControl(method = 'cv', number = 5, verboseIter = TRUE)  # establish training regim, 5-fold cross validation
numFolds <- trainControl(method = 'none', verboseIter = TRUE)            # speeding up the fit of a final model based on k-fold stability
method_mod <- 'mxnet'
tuneGrid <- expand.grid(
   layer1 = 10
  ,layer2 = 10
  ,layer3 = 10
  ,learning.rate = 0.02
  #, learning.rate = 2e-6 # stock
  ,momentum = 0.9
  ,dropout = .2
  ,activation = "relu"
)

library(mxnet)

# land tune
terre <- train(y = LAND_trainer[,"terrestrial_richness"],
               x = LAND_trainer[,Land_varnames],
               method     = method_mod, 
               preProcess = c('center', 'scale'), 
               trControl  = numFolds, 
               tuneGrid   = tuneGrid
               ,array.batch.size = 128
               ,num.round = 10
               
)
terre

terre$results

# sea tune
oceaus <-train(y = OCEAN_trainer[,"marine_richness"],
               x = OCEAN_trainer[,sea_varnames],
               method     = method_mod, 
               preProcess = c('center', 'scale'), 
               trControl  = numFolds, 
               tuneGrid   = tuneGrid
               ,array.batch.size = 128
               ,num.round = 10
)
oceaus

oceaus$results

##### Quick testset evaluation #####

#LAND
# RMSE
caret::RMSE(predict(object = terre, newdata = LAND_trainer),LAND_trainer[,"terrestrial_richness"])
caret::RMSE(predict(object = terre, newdata = LAND_test),LAND_test[,"terrestrial_richness"])
# R2
1 - (sum((LAND_trainer[,"terrestrial_richness"]-predict(object = terre, newdata = LAND_trainer) )^2)/sum((LAND_trainer[,"terrestrial_richness"]-mean(LAND_trainer[,"terrestrial_richness"]))^2))
1 - (sum((LAND_test[,"terrestrial_richness"]-predict(object = terre, newdata = LAND_test) )^2)/sum((LAND_test[,"terrestrial_richness"]-mean(LAND_test[,"terrestrial_richness"]))^2))

# OCEAN
# RMSE
caret::RMSE(predict(object = oceaus, newdata = OCEAN_trainer),OCEAN_trainer[,"marine_richness"])
caret::RMSE(predict(object = oceaus, newdata = OCEAN_test),OCEAN_test[,"marine_richness"])
# R2
1 - (sum((OCEAN_trainer[,"marine_richness"]-predict(object = oceaus, newdata = OCEAN_trainer) )^2)/sum((OCEAN_trainer[,"marine_richness"]-mean(OCEAN_trainer[,"marine_richness"]))^2))
1 - (sum((OCEAN_test[,"marine_richness"]-predict(object = oceaus, newdata = OCEAN_test) )^2)/sum((OCEAN_test[,"marine_richness"]-mean(OCEAN_test[,"marine_richness"]))^2))


#######################
# Visualize the model #
#######################

library(doParallel)  # load the parallel backend
cl <- makeCluster(6)  # use 4 workers
registerDoParallel(cl)  # register the parallel backend

##### quick 2D partial dependencies #####
# these are the examples use in the manuscript 4 panel example

#D2pdp_oceaus <- partial(oceaus, pred.var = c("elev", "SST_means"), grid.resolution = 50, progress = "text")
#D2pdp_terre <- partial(terre, pred.var = c("elev", "SAT_means"), grid.resolution = 50,  progress = "text")
D2pdp_oceaus <- partial(oceaus, pred.var = c("elev", "SST_means"), grid.resolution = 20, progress = "text", parallel = T)
D2pdp_terre <- partial(terre, pred.var = c("elev", "SAT_means"), grid.resolution = 20,  progress = "text", parallel = T)

# other vars
D2pdp_oceaus2 <- partial(oceaus, pred.var = c("elev", "AMP_1Y"), grid.resolution = 20, progress = "text", parallel = T)
D2pdp_terre2 <- partial(terre, pred.var = c("elev", "AMP_1Y"), grid.resolution = 20,  progress = "text", parallel = T)

# last round of vars
D2pdp_oceaus3 <- partial(oceaus, pred.var = c("elev", "OXY_CV"), grid.resolution = 20, progress = "text", parallel = T)
D2pdp_terre3 <- partial(terre, pred.var = c("elev", "P_CV"), grid.resolution = 20,  progress = "text", parallel = T)

D2pdp_oceaus <- as.data.frame(D2pdp_oceaus)
D2pdp_terre <- as.data.frame(D2pdp_terre)
D2pdp_oceaus2 <- as.data.frame(D2pdp_oceaus2)
D2pdp_terre2 <- as.data.frame(D2pdp_terre2)
D2pdp_oceaus3 <- as.data.frame(D2pdp_oceaus3)
D2pdp_terre3 <- as.data.frame(D2pdp_terre3)


# surfaces
land_pdep <-  ggplot()+
  geom_tile(data = D2pdp_terre, aes(elev,SAT_means,fill = yhat))+
  geom_rug(data = sample_n(LAND_trainer,1000),aes(elev,SAT_means), alpha = .09, position = "jitter")+
  scale_fill_distiller(palette = "RdYlBu")+
  coord_fixed()+
  labs(x = "elevation", y = "mean surface air temperature", title = "Terrestrial")+
  themeo
ocean_pdep <- ggplot()+
  geom_tile(data = D2pdp_oceaus, aes(elev*-1,SST_means,fill = yhat))+
  geom_rug(data = sample_n(OCEAN_trainer,1000),aes(elev*-1,SST_means), alpha = .09, position = "jitter")+
  scale_fill_distiller(palette = "RdYlBu")+
  coord_fixed()+
  labs(x = "elevation", y = "mean sea surface temperature", title = "Marine")+
  themeo
land_pdep2 <- ggplot()+
  geom_tile(data = D2pdp_terre2, aes(elev,AMP_1Y,fill = yhat))+
  geom_rug(data = sample_n(LAND_trainer,1000),aes(elev,AMP_1Y), alpha = .09, position = "jitter")+
  scale_fill_distiller(palette = "PRGn")+
  coord_fixed()+
  labs(x = "elevation", y = "annual PP seasonality", title = "Terrestrial")+
  themeo
ocean_pdep2 <-ggplot()+
  geom_tile(data = D2pdp_oceaus2, aes(elev*-1,AMP_1Y,fill = yhat))+
  geom_rug(data = sample_n(OCEAN_trainer, 1000),aes(elev*-1,AMP_1Y), alpha = .09, position = "jitter")+
  scale_fill_distiller(palette = "PRGn")+
  coord_fixed()+
  labs(x = "elevation", y = "annual PP seasonality", title = "Marine", legend = "species richness")+
  themeo
land_pdep3 <- ggplot()+
  geom_tile(data = D2pdp_terre3, aes(elev,P_CV,fill = yhat))+
  geom_rug(data = sample_n(LAND_trainer,1000),aes(elev,AMP_1Y), alpha = .09, position = "jitter")+
  scale_fill_viridis()+
  coord_fixed()+
  labs(x = "elevation", y = "precipitation", title = "Terrestrial")+
  themeo
ocean_pdep3 <-ggplot()+
  geom_tile(data = D2pdp_oceaus3, aes(elev*-1,OXY_CV,fill = yhat))+
  geom_rug(data = sample_n(OCEAN_trainer, 1000),aes(elev*-1,AMP_1Y), alpha = .09, position = "jitter")+
  scale_fill_viridis()+
  coord_fixed()+
  labs(x = "elevation", y = "oxygen", title = "Marine", legend = "species richness")+
  themeo

grid.arrange(land_pdep,ocean_pdep,land_pdep2,ocean_pdep2,land_pdep3,ocean_pdep3, nrow = 3)

# perspectives
land_pdep <- ggplot(D2pdp_terre,aes(elev,yhat,color = SAT_means))+
  scale_color_distiller(palette = "RdYlBu")+
  geom_point(size = 4,show.legend = F#, shape = 22, color = "white", alpha = .2
             )+themeo

ocean_pdep <- ggplot(D2pdp_oceaus,aes(elev*-1,yhat, color = SST_means))+
  scale_color_distiller(palette = "RdYlBu")+
  geom_point(size = 4,show.legend = F#, alpha = .2#, shape = 22, color = "white"
             )+themeo

land_pdep2 <- ggplot(D2pdp_terre2,aes(elev,yhat,color = AMP_1Y))+
  scale_color_distiller(palette = "PRGn", direction = 1)+
  geom_point(size = 4,show.legend = F#, shape = 22, color = "white", alpha = .2
  )+themeo

ocean_pdep2 <-ggplot(D2pdp_oceaus2,aes(elev*-1,yhat, color = AMP_1Y))+
  scale_color_distiller(palette = "PRGn", direction = 1)+
  geom_point(size = 4,show.legend = F#, alpha = .2#, shape = 22, color = "white"
  )+themeo

land_pdep3 <- ggplot(D2pdp_terre3,aes(elev,yhat,color = P_CV))+
  scale_color_viridis()+
  geom_point(size = 4,show.legend = F#, shape = 22, color = "white", alpha = .2
  )+themeo

ocean_pdep3 <-ggplot(D2pdp_oceaus3,aes(elev*-1,yhat, color = OXY_CV))+
  scale_color_viridis()+
  geom_point(size = 4,show.legend = F#, alpha = .2#, shape = 22, color = "white"
  )+themeo

grid.arrange(land_pdep,ocean_pdep,land_pdep2,ocean_pdep2,land_pdep3,ocean_pdep3, nrow = 3)

stopCluster(cl)


# grid of all variables as bivariate partial dependency plots below
# plot one on top of another

# centered on 0
par(mfrow = c(5,6))
for(i in 1:length(Land_varnames)) {
  
  Land_varname <- Land_varnames[i]
  pd_land <- partial(terre, pred.var = Land_varname, parallel = T)
  pd_land$yhat <- scale(pd_land$yhat, center = T, scale = F)
  plot(pd_land[,1],pd_land[,2],type = 'l', xlab = colnames(pd_land)[1],ylab = 'tp', col = "#1a9850",ylim = c(-0.1,0.1))

  sea_varname <- sea_varnames[i]
  pd_sea <- partial(oceaus, pred.var = sea_varname, parallel = T)
  pd_sea$yhat <- scale(pd_sea$yhat, center = T, scale = F)
  lines(pd_sea[,1],pd_sea[,2],type = 'l', xlab = colnames(pd_sea)[1],ylab = 'tp', col = "#4575b4")
}

# fixed first value to zero
par(mfrow = c(5,6))
for(i in 1:length(Land_varnames)) {
  
  Land_varname <- Land_varnames[i]
  pd_land <- partial(terre, pred.var = Land_varname)
  pd_land$yhat <- pd_land$yhat - pd_land$yhat[1]
  plot(pd_land[,1],pd_land[,2],type = 'l', xlab = colnames(pd_land)[1],ylab = 'tp', col = "#1a9850",ylim = c(-0.1,0.1))
  
  sea_varname <- sea_varnames[i]
  pd_sea <- partial(oceaus, pred.var = sea_varname)
  pd_sea$yhat <- pd_sea$yhat - pd_sea$yhat[1]
  lines(pd_sea[,1],pd_sea[,2],type = 'l', xlab = colnames(pd_sea)[1],ylab = 'tp', col = "#4575b4")
}
  
  
# Locating interaction potential with individual conditional expectation plots
# classic ice plot vis
par(mfrow = c(5,6))
for(i in 1:length(Land_varnames)){
  pd_ice <- ice(object = terre, X = LAND_trainer[,Land_varnames], y = LAND_trainer[,"terrestrial_richness"], predictor = Land_varnames[i], 
                frac_to_build = .01,
                verbose = T) %>% 
    plot(x_quantile = F, plot_pdp = TRUE, frac_to_plot = 1, centered = T, pts_preds_size = .5 )}

par(mfrow = c(5,6))
for(i in 1:length(sea_varnames)){
  pd_ice <- ice(object = oceaus, X = OCEAN_trainer[,sea_varnames], y = OCEAN_trainer[,"marine_richness"], predictor = sea_varnames[i], 
                frac_to_build = .01,
                verbose = T) %>% 
    plot(x_quantile = F, plot_pdp = TRUE, frac_to_plot = 1, centered = T, pts_preds_size = .5 )}

# alternatives
# LAND run with pdp
QuantilesTorF = F
full_df_land1 <- NULL
for( i in 1:length(Land_varnames)){
  pd1 <- partial(terre, pred.var = Land_varnames[i], ice = T, quantiles = QuantilesTorF
                 #, probs = c(seq(from = 0, to = 1, length = 100))
                 )
  pd1$pred_var <- Land_varnames[i]
  colnames(pd1) <- c("x","yhat","yhat.id","predvar")
  full_df_land1<-rbind(full_df_land1,pd1)}

str(full_df_land1)

#############################################################################
# Sample subset of full pdep data set to avoid overplotting and speed it up #
#############################################################################

# generate random sample of yhat.id
full_df_land1$predvar <- as.factor(full_df_land1$predvar)
test <- full_df_land1 %>% filter(yhat.id == c(seq(1,11, 1)))
ggplot(test, aes(x, yhat, group = yhat.id))+geom_line(col = "black", size = 0.1)+facet_wrap(~as.factor(predvar))

# OCEAN run with pdp
QuantilesTorF = F
full_df_sea1 <- NULL
for( i in 1:length(sea_varnames)){
  pd1 <- partial(oceaus, pred.var = sea_varnames[i], ice = T, quantiles = QuantilesTorF
                 #, probs = c(seq(from = 0, to = 1, length = 100))
  )
  pd1$pred_var <- sea_varnames[i]
  colnames(pd1) <- c("x","yhat","yhat.id","predvar")
  full_df_sea1<-rbind(full_df_sea1,pd1)}

str(full_df_sea1)

#############################################################################
# Sample subset of full pdep data set to avoid overplotting and speed it up #
#############################################################################
# generate random sample of yhat.id
full_df_sea1$predvar <- as.factor(full_df_sea1$predvar)
test <- full_df_sea1 %>% filter(yhat.id == c(seq(1,11, 1)))
ggplot(test, aes(x, yhat, group = yhat.id))+geom_line(col = "black", size = 0.1)+facet_wrap(~as.factor(predvar))


#############################################################################
#############################################################################

#### #### #### #### #### #### #### #### 
#### spatially plotting model error ####
#### #### #### #### #### #### #### #### 

# observed vs modeled
modeledL <- predict(object = terre, newdata = LAND_test)
LAND_submission$modeled <- predict(object = terre, newdata = LAND_submission)
LAND_submission$mod_error <- LAND_submission$terrestrial_richness - LAND_submission$modeled


# test set land
L_test <- ggplot()+
  geom_point(aes(x = LAND_test$terrestrial_richness, y = modeledL ), size = 1)+
  geom_line(aes(seq(from = 0,to = 1,length = 2), seq(from = 0,to = 1, length = 2)),color = "red", size = 1)+
  scale_y_continuous(expand = c(0,0))+ 
  scale_x_continuous(expand = c(0,0))+
  labs(x = "Observed", y = "Predicted")+
  themeo

# full dataset land
# size = .01
L_train <- ggplot()+
  geom_point(aes(x = LAND_submission$terrestrial_richness, y = LAND_submission$modeled, fill =  LAND_submission$mod_error), show.legend = F,size = 2.5, shape = 21)+
  geom_line(aes(seq(from = 0,to = 1,length = 2), seq(from = 0,to = 1, length = 2)),color = "red", size = 1)+
  scale_y_continuous(expand = c(0,0))+ 
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradientn(colours = rev(c("#2166ac","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#b2182b")))+
  labs(x = "Observed", y = "Predicted")+
  themeo

# ocean
modeledO <- predict(object = oceaus, newdata = OCEAN_test)
SEA_submission$modeled <- predict(object = oceaus, newdata = SEA_submission)
SEA_submission$mod_error <- SEA_submission$marine_richness - SEA_submission$modeled

# look at inputs
#melt(SEA_submission, id.vars = c("x","y"))  %>% ggplot(aes(x = x, y = y, fill = value))+facet_wrap(~variable, scales = "free")+geom_raster(show.legend = F)+coord_fixed()


# test set ocean
O_test <- ggplot()+
  geom_point(aes(x = OCEAN_test$marine_richness, y = modeledO ), size = 1)+
  geom_line(aes(seq(from = 0,to = 1,length = 2), seq(from = 0,to = 1, length = 2)),color = "red", size = 1)+
  scale_y_continuous(expand = c(0,0))+ 
  scale_x_continuous(expand = c(0,0))+
  labs(x = "Observed", y = "Predicted")+
  themeo

# full dataset ocean
O_train <- ggplot()+
  geom_point(aes(x = SEA_submission$marine_richness, y = SEA_submission$modeled, fill =  SEA_submission$mod_error), show.legend = F,size = 2.5, shape = 21)+
  geom_line(aes(seq(from = 0,to = 1,length = 2), seq(from = 0,to = 1, length = 2)),color = "red", size = 1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradientn(colours = rev(c("#2166ac","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#b2182b")))+
  labs(x = "Observed", y = "Predicted")+
  themeo

grid.arrange(L_train,L_test,O_train,O_test) 

# rescale both domain errors.
SEA_submission[,"mod_error"] <- scales::rescale(x = SEA_submission[,"mod_error"], to = c(0,1))
LAND_submission[,"mod_error"] <- scales::rescale(x = LAND_submission[,"mod_error"], to = c(0,1))

# overlayed, red means the model predicted more richness than observed, blue means the model under predicted observed richness
ggplot(SEA_submission)+
  geom_raster(data = SEA_submission, aes(x,y, fill = mod_error))+coord_fixed()+
  geom_raster(data = LAND_submission, aes(x,y, fill = mod_error))+coord_fixed()+
  geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .25)+
  scale_fill_gradientn(colours = rev(c("#154268","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#6b0f1f")))+
  
  #scale_fill_gradientn(colours = rev(c("#154268","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#6b0f1f")))+
  #scale_fill_gradientn(colours = rev(c("#2166ac","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#b2182b")))+
  #scale_fill_gradientn(colours = rev(c("#2166ac","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#b2182b")))+
  #scale_fill_gradientn(colours = rev(c("#2166ac","#2166ac","#4575b4","#74add1","#d1e5f0","white","white","white","#fddbc7","#f46d43","#d73027","#b2182b","#b2182b")))+
  #scale_fill_distiller(palette = "RdBu", direction = 1)+
  
  themeo+
  scale_x_continuous(expand = c(-0.005,-0.005)) +
  scale_y_continuous(expand = c(-0.01,-0.01)) +
  #theme(axis.text = element_blank(),
  #    axis.title = element_blank())+
  ggtitle("A. model residuals w/o SAC covariate w error rescaled")+
  coord_fixed()


#########################################################################################
###### This is where I am checking the model residuals for spatial autocorrelation ######
#########################################################################################

# Based on the degree of spatial autocorrelation observed and the sill of the semi-variogram
# I will select an appropriate neighborhood, from there I will calculate a spatial autocovariate
# from which I will add to the full dataframe and refit the model
# Then I can recalculate model error and plot the residuals.

####################################################
####### Spatial Autocorrelation Assessment ######### 
#################################################### 


# variogram SEA
SEA_points <- SEA_submission
coordinates(SEA_points) = ~ x + y
SEAVariogram=variogram(mod_error~1, data=SEA_points)
f_sea <- fit.variogram(SEAVariogram,vgm(c("Exp", "Mat", "Sph")))
f_sea
plot(SEAVariogram, model=f_sea) 

# variogram LAND
LAND_points <- LAND_submission
coordinates(LAND_points) = ~ x + y
LANDVariogram=variogram(mod_error~1, data=LAND_points)
f_land <- fit.variogram(LANDVariogram,vgm(c("Exp", "Mat", "Sph")))
f_land
plot(LANDVariogram, model=f_land) 


#######################################################################
#  Generate an auto covariate term to handle spatial autocorrelation  # # 
#######################################################################
# option-B is a revised selection from historical option-W that Doorman intially utilized in refuting autocov terms
# B is sutiable by assessment of Guillera et al 2015, ensures symettry of the neighborhood matrix (Bardos et al 2015)
# nbs is a function of the range approximated by the variogram model, though its still computationally extremely demanding
LAND_submission$auto_cov <- autocov_dist(
  z = LAND_submission$mod_error , 
  xy = SpatialPoints(data.frame(x = LAND_submission$x, y = LAND_submission$y)), 
  nbs = floor(f_land$range[2]), 
  type = "inverse", 
  zero.policy = NULL, 
  style = "B")

SEA_submission$auto_cov <- autocov_dist(
  z = SEA_submission$mod_error , 
  xy = SpatialPoints(data.frame(x = SEA_submission$x,y = SEA_submission$y)), 
  nbs = floor(f_sea$range[2]), 
  type = "inverse", 
  zero.policy = NULL, 
  style = "B")


# checking nature of correlation between auto_cov and dependent var
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(64)
k <- kde2d(LAND_submission$auto_cov, LAND_submission$terrestrial_richness, n=200)
cor(LAND_submission$auto_cov, LAND_submission$terrestrial_richness)
image(k, col=r)
k <- kde2d(SEA_submission$auto_cov, SEA_submission$marine_richness, n=200)
cor(SEA_submission$auto_cov, SEA_submission$marine_richness)
image(k, col=r)

# redo varnames after developing covariate
Land_varnames <- colnames(LAND_submission)[4:length(LAND_submission)] # note these indexes, adjust as needed when including new pred/response vars
sea_varnames  <- colnames(SEA_submission)[4:length(SEA_submission)]


#autocov plotting
ggplot(SEA_submission)+
  geom_raster(data = SEA_submission, aes(x,y, fill = auto_cov))+coord_fixed()+
  geom_raster(data = LAND_submission, aes(x,y, fill = auto_cov))+coord_fixed()+
  #geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .25)+
  scale_fill_distiller(palette = 'Spectral')+
  themeo+
  scale_x_continuous(expand = c(-0.005,-0.005)) +
  scale_y_continuous(expand = c(-0.01,-0.01)) +
  #theme(axis.text = element_blank(),
  #    axis.title = element_blank())+
  ggtitle("A.residual autocovariate values based on nbs as ridge of semi-variogram")+
  coord_fixed()



# model richness based on learned relationships
ggplot()+
  geom_raster(data = SEA_submission, aes(x,y, fill = modeled))+coord_fixed()+
  geom_raster(data = LAND_submission, aes(x,y, fill = modeled))+coord_fixed()+
 # geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .5)+
  scale_fill_distiller(palette = "Spectral", direction = -1)
ggplot()+
  geom_raster(data = SEA_submission, aes(x,y, fill = marine_richness))+coord_fixed()+
  geom_raster(data = LAND_submission, aes(x,y, fill = terrestrial_richness))+coord_fixed()+
  # geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .5)+
  scale_fill_distiller(palette = "Spectral", direction = -1)


 # write.csv(SEA_submission,file = "/Users/tgagne/Biodiversity2018/data/Gabriel/SEA_submission.csv")
 # write.csv(LAND_submission,file = "/Users/tgagne/Biodiversity2018/data/Gabriel/LAND_submission.csv")
 
    
    
    
    