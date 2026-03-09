##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SBE/Danum LiDAR Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load libraries
library(terra)
library(sf)
library(MultiscaleDTM)

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%%##

# Import CHM 2020
SBE_2020_CHM <- rast("SBE_2020_CHM.tif")

# Import CHM 2013
SBE_2013_CHM <- rast("SBE_2013_CHM.tif")

# Import DEM
SBE_DEM <- rast("SBE_DEM.tif")

# Load shapefile and covert to sf object
SBE_plots <- vect("SBE_plots.shp")
SBE_plots_sf <- st_as_sf(SBE_plots)

# Crop rasters to the same size
SBE_2013_CHM<-crop(SBE_2013_CHM,SBE_2020_CHM)
SBE_2020_CHM<-crop(SBE_2020_CHM, SBE_2013_CHM)
SBE_DEM <-crop(SBE_DEM, SBE_2020_CHM)

# Reproject SBE_plots to match
SBE_plots <- project(SBE_plots,"epsg:32650")

# Plot
plot(SBE_DEM)
plot(SBE_2020_CHM)
plot(SBE_2013_CHM)
plot(SBE_plots, add=T,border="white")

#----------------------#
### Save the Rasters ###
#----------------------#

writeRaster(SBE_2013_CHM, "SBE_2013_CHM_crop.tif", overwrite=T)
writeRaster(SBE_2020_CHM, "SBE_2020_CHM_crop.tif", overwrite=T)
writeRaster(SBE_DEM, "SBE_DEM_crop.tif", overwrite=T)

##%%%%%%%%%%%%%##
##### DANUM #####
##%%%%%%%%%%%%%##

# Import CHM 2020
Danum_2020_CHM<-rast("Danum_2020_CHM.tif")

# Import CHM 2013
Danum_2013_CHM<-rast("Danum_2013_CHM.tif")

# Import DEM
Danum_DEM <- rast("Danum_DEM.tif")

# Load shapefile and covert to sf object
Danum_plots <- vect("Danum_plots.shp")
Danum_plots_sf<- st_as_sf(Danum_plots)

# Crop rasters to the same size
Danum_2020_CHM<-crop(Danum_2020_CHM,Danum_2013_CHM)
Danum_2013_CHM<-crop(Danum_2013_CHM, Danum_2020_CHM)
Danum_DEM <-crop(Danum_DEM, Danum_2020_CHM)

# Reproject SBE_plots to match
Danum_plots <- project(Danum_plots,"epsg:32650")

# Add danum shapefile to plot
plot(Danum_DEM)
plot(Danum_2020_CHM)
plot(Danum_2013_CHM)
plot(Danum_plots, add=T,border="white")

#----------------------#
### Save the Rasters ###
#----------------------#

writeRaster(Danum_2013_CHM, "Danum_2013_CHM_crop.tif", overwrite=T)
writeRaster(Danum_2020_CHM, "Danum_2020_CHM_crop.tif", overwrite=T)
writeRaster(Danum_DEM, "Danum_DEM_crop.tif", overwrite=T)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### CREATE TOPOGRAPHIC RASTERS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%###

# Re-sample at 5m resolution for speed
SBE_DEM_5<-aggregate(SBE_DEM, fact=5, fun=mean)

#-------------------------------#
### Calculate terrain metrics ###
#-------------------------------#

# Calculate terrain metrics
SBE_Slope<- terrain(SBE_DEM_5, v="slope", unit='degrees')
SBE_TPI <- TPI(SBE_DEM_5, shape = "rectangle")
SBE_TRI <- terrain(SBE_DEM_5, v="TRI")
SBE_flow_dir <- terrain(SBE_DEM_5, v = "flowdir")
SBE_Aspect<- terrain(SBE_DEM_5, v="aspect", unit='degrees')
SBE_DEM_hillshade <- shade(SBE_Slope, SBE_Aspect, angle = 40, direction = 270)
SBE_Roughness<- terrain(SBE_DEM_5, v = "roughness")

# Plot the rasters
plot(SBE_Slope)
plot(SBE_TPI)
plot(SBE_flow_dir)
plot(SBE_Aspect)
plot(SBE_DEM_hillshade)
plot(SBE_Roughness)

#----------------------#
### Save the Rasters ###
#----------------------#

writeRaster(SBE_Slope, "SBE_Slope.tif", overwrite=T)
writeRaster(SBE_TPI, "SBE_TPI.tif", overwrite=T)
writeRaster(SBE_TRI, "SBE_TRI.tif", overwrite=T)
writeRaster(SBE_flow_dir, "SBE_flow_dir.tif", overwrite=T)
writeRaster(SBE_Aspect, "SBE_Aspect.tif", overwrite=T)
writeRaster(SBE_Roughness, "SBE_Roughness.tif",overwrite=T)

##%%%%%%%%%%%%%##
##### Danum #####
##%%%%%%%%%%%%%##

#-------------------------------#
### Calculate terrain metrics ###
#-------------------------------#

# Resample at 5m resolution for speed
Danum_DEM_5_danum<-aggregate(Danum_DEM, fact=5, fun=mean)

# Calculate terrain metrics
Danum_Slope <- terrain(Danum_DEM_5_danum, v="slope", unit='degrees')
Danum_TPI <- TPI(Danum_DEM_5_danum, shape = "rectangle")
Danum_TRI <- terrain(Danum_DEM_5_danum, v="TRI")
Danum_flow_dir <- terrain(Danum_DEM_5_danum, v = "flowdir")
Danum_Aspect<- terrain(Danum_DEM_5_danum, v="aspect", unit='degrees')
Danum_DEM_hillshade <- shade(Danum_Slope, Danum_Aspect, angle = 40, direction = 270)
Danum_Roughness<- terrain(Danum_DEM_5_danum, v = "roughness")

#----------------------#
### Save the Rasters ###
#----------------------#

writeRaster(Danum_Slope, "Danum_Slope.tif", overwrite=T)
writeRaster(Danum_TPI, "Danum_TPI.tif", overwrite=T)
writeRaster(Danum_TRI, "Danum_TRI.tif", overwrite=T)
writeRaster(Danum_Aspect, "Danum_Aspect.tif", overwrite=T)
writeRaster(Danum_flow_dir, "Danum_flow_dir.tif", overwrite=T)
writeRaster(Danum_Roughness, "Danum_Roughness.tif", overwrite=T)
