#### code/map/interpolate_sulfate.R ####
# Benjamin D. Peterson

#### Get cleaned up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(dismo) # Species Distribution Modeling
library(gstat)
library(mapStats)
library(maptools)
library(raster)
library(rgdal)
library(rgeos)
library(spdep)
library(patchwork)
library(readxl)
library(tidyverse)
library(viridisLite)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Read in WCA boundaries ####
wca_boundaries <- readOGR(dsn = "dataEdited/gisData/shape_files_needed",
                          layer = "wca_bnds")
crs(wca_boundaries)



#### Read in geoTIFF data from Carl Fitz ####
# Read in data and check out projection
GDALinfo("dataEdited/sulfateMap/USGS_Base.MeanPOS.SO4concSfAvg20001223.tif")
sulfate.map <- raster("dataEdited/sulfateMap/USGS_Base.MeanPOS.SO4concSfAvg20001223.tif")
crs(sulfate.map)

# Made WCA boundaries projections match sulfate map
crs(wca_boundaries)
wca_boundaries <- spTransform(wca_boundaries,
                              crs(sulfate.map))
crs(wca_boundaries)


#### Trim up the sulfate map to area that we want ####
plot(sulfate.map)
sulfate.map.cropped <- crop(sulfate.map,
                            wca_boundaries)
plot(sulfate.map.cropped)


#### Replace NA's with 0, so we're not left with the open white spots ####
sulfate.map.cropped@data@values[which(sulfate.map.cropped@data@values == "NaN")] <- 1
plot(sulfate.map.cropped)
sulfate.map.cropped.masked <- mask(sulfate.map.cropped,
                                   wca_boundaries)
plot(sulfate.map.cropped.masked)


#### Transform data ####
# sulfate.map.cropped.masked@data@values <- log(sulfate.map.cropped.masked@data@values, 10)


#### Set up the break points ####
max.sulfate <- max(sulfate.map.cropped.masked@data@values,
                   na.rm = TRUE)
min.sulfate <- min(sulfate.map.cropped.masked@data@values,
                   na.rm = TRUE)
# Never gets below 1, must be the way the model is set up (that's
# why I went back and converted the NAs to 1 instead of 0).
cuts  <- seq(min.sulfate,
             max.sulfate,
             by = 1) #set breaks


#### Check out the map ####
plot(sulfate.map.cropped.masked,
     breaks = cuts,
     col = viridis(length(cuts))[length(cuts):1])


#### Set colors for when it's brought into qGIS ####
sulfate.map.cropped.masked@legend@values <- cuts
sulfate.map.cropped.masked@legend@colortable <- viridis(length(cuts))[length(cuts):1]


#### Save out raster ####
writeRaster(sulfate.map.cropped.masked,
            'dataEdited/sulfateMap/sulfate_map_edited.tif',
            options = c('TFW=YES'),
            overwrite = TRUE)
