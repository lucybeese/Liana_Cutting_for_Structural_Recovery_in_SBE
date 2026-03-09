##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SBE/Danum LiDAR Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load libraries
library(terra)
library(sf)
library(MultiscaleDTM)
library(dplyr)

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%%##

# Import CHM 2020
SBE_2020_CHM <- rast("SBE_2020_CHM_crop.tif")

# Import CHM 2013
SBE_2013_CHM <- rast("SBE_2013_CHM_crop.tif")

# Import DEM
SBE_DEM <- rast("SBE_DEM_crop.tif")

# Load shapefile and covert to sf object
SBE_plots <- vect("SBE_plots.shp")
SBE_plots_sf <- st_as_sf(SBE_plots)

# Reproject SBE_plots to match
SBE_plots <- project(SBE_plots,"epsg:32650")

# Load species attributes
species <- st_read('Data_Attributes_SBE.gpkg') 

##%%%%%%%%%%%%%##
##### DANUM #####
##%%%%%%%%%%%%%##

# Import CHM 2020
Danum_2020_CHM<-rast("Danum_2020_CHM_crop.tif")

# Import CHM 2013
Danum_2013_CHM<-rast("Danum_2013_CHM_crop.tif")

# Import DEM
Danum_DEM <- rast("Danum_DEM_crop.tif")

# Load shapefile and covert to sf object
Danum_plots <- vect("Danum_plots.shp")
Danum_plots_sf<- st_as_sf(Danum_plots)

# Reproject SBE_plots to match
Danum_plots <- project(Danum_plots,"epsg:32650")

#------------------------------#
### Load the Rasters (DANUM) ###
#------------------------------#

Danum_Slope <- rast('Danum_Slope.tif')
Danum_TPI <- rast('Danum_TPI.tif')
Danum_TRI <- rast('Danum_TRI.tif')
Danum_Aspect <- rast('Danum_Aspect.tif')
Danum_flow_dir <- rast('Danum_flow_dir.tif')
Danum_Roughness <- rast('Danum_Roughness.tif')

#---------------------------#
### Load the Rasters (SBE)###
#---------------------------#

SBE_Slope <- rast('SBE_Slope.tif')
SBE_TPI<- rast('SBE_TPI.tif')
SBE_TRI<- rast('SBE_TRI.tif')
SBE_flow_dir<- rast('SBE_flow_dir.tif')
SBE_Aspect<- rast('SBE_Aspect.tif')
SBE_Roughness<- rast('SBE_Roughness.tif')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### EXTRACT METRICS FOR EACH PLOT ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%%##

#---------------------------------------------#
### CHM height change between 2013 and 2020 ###
#---------------------------------------------#

SBE_CHM_change<-SBE_2020_CHM-SBE_2013_CHM
plot(SBE_CHM_change)

#-----------------------------------------------------------#
### Extract canopy and topographic data for each SBE plot ###
#-----------------------------------------------------------#

# Create data frame to store data
SBE_plot_data<-data.frame(matrix(dim(SBE_plots)[1],20,data=NA))
names(SBE_plot_data)<-c("Site","Plot_id","Treatment","TCH_2020","TCH_2013","TCH_change","Gap_frac_20m_2020","Gap_frac_20m_2013","Frac_Over_50m","Height_cv_2020",
                        "Height_cv_2013","Elevation","Slope","TRI","Max_Canopy_Height_2020","Max_Canopy_Height_2013","Aspect","TPI", 
                        "flow_dir", 'Roughness')
SBE_plot_data$Site<-"SBE"
SBE_plot_data$Plot_id<-SBE_plots$Plot_id
SBE_plot_data$Treatment<-SBE_plots$Treatment

# Run loop across each plot and extract data
for (i in 1:dim(SBE_plot_data)[1]){
  
  # CHM 2020
  SBE_2020_CHM_crop<-crop(SBE_2020_CHM, SBE_plots[i,])
  # Mean Top of Canopy Height
  SBE_plot_data[i,"TCH_2020"]<-global(SBE_2020_CHM_crop, fun='mean', na.rm=TRUE) 
  # Coefficient of variation
  SBE_plot_data[i,"Height_cv_2020"]<-global(SBE_2020_CHM_crop, fun='sd', na.rm=TRUE)/global(SBE_2020_CHM_crop, fun='mean', na.rm=TRUE) 
  # % pixels below 20m threshold (crude gaps)
  SBE_plot_data[i,"Gap_frac_20m_2020"]<-length(which(as.vector(na.omit(as.vector(SBE_2020_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(SBE_2020_CHM_crop)))) 
  
  # CHM 2013 Gap fraction at 20m and Coefficient of variation
  SBE_2013_CHM_crop<-crop(SBE_2013_CHM,SBE_plots[i,])
  SBE_plot_data[i,"TCH_2013"]<-global(SBE_2013_CHM_crop, fun='mean', na.rm=TRUE)
  SBE_plot_data[i,"Gap_frac_20m_2013"]<-length(which(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(SBE_2013_CHM_crop))))
  SBE_plot_data[i,"Height_cv_2013"]<-global(SBE_2013_CHM_crop, fun='sd', na.rm=TRUE)/global(SBE_2013_CHM_crop, fun='mean', na.rm=TRUE)
  # % pixels above 50m threshold (proxy for logging intensity)
  SBE_plot_data[i,"Frac_Over_50m"]<-length(which(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))>50))/length(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))) 
  
  # Canopy Height change
  SBE_CHM_change_crop<-crop(SBE_CHM_change,SBE_plots[i,])
  SBE_plot_data[i,"TCH_change"]<-global(SBE_CHM_change_crop, fun='mean', na.rm=TRUE)
  
  # Elevation
  SBE_DEM_crop<-crop(SBE_DEM,SBE_plots[i,])
  SBE_plot_data[i,"Elevation"]<-global(SBE_DEM_crop, fun='mean', na.rm=TRUE)
  
  # Aspect
  Aspect_crop<-crop(SBE_Aspect,SBE_plots[i,])
  SBE_plot_data[i,"Aspect"]<-global(Aspect_crop, fun='mean', na.rm=TRUE)
  
  #TPI
  TPI_crop<-crop(SBE_TPI,SBE_plots[i,])
  SBE_plot_data[i,"TPI"]<-global(TPI_crop, fun='mean', na.rm=TRUE)
  
  # Slope
  Slope_crop<-crop(SBE_Slope,SBE_plots[i,])
  SBE_plot_data[i,"Slope"]<-global(Slope_crop, fun='mean', na.rm=TRUE)
  
  # TRI
  TRI_crop<-crop(SBE_TRI,SBE_plots[i,])
  SBE_plot_data[i,"TRI"]<-global(TRI_crop, fun='mean', na.rm=TRUE)
  
  # Flow-dir
  flow_dir_crop<-crop(SBE_flow_dir,SBE_plots[i,])
  SBE_plot_data[i,"flow_dir"]<-global(flow_dir_crop, fun='mean', na.rm=TRUE)
  
  # Roughness
  roughness_crop<-crop(SBE_Roughness,SBE_plots[i,])
  SBE_plot_data[i,"Roughness"]<-global(roughness_crop, fun='mean', na.rm=TRUE)
  
  # Max Canopy Height
  
  # 2020
  MCH_crop_2020<-crop(SBE_2020_CHM,SBE_plots[i,])
  values_MCH_crop_2020 <- values(MCH_crop_2020, na.rm = TRUE)
  Max_Canopy_Height_2020<-stats::quantile(MCH_crop_2020, probs=0.99, na.rm = TRUE)
  SBE_plot_data[i,"Max_Canopy_Height_2020"]<-Max_Canopy_Height_2020
  
  # 2013
  MCH_crop_2013<-crop(SBE_2013_CHM,SBE_plots[i,])
  values_MCH_crop_2013 <- values(MCH_crop_2013, na.rm = TRUE)
  Max_Canopy_Height_2013 <- stats::quantile(values_MCH_crop_2013, probs = 0.99,na.rm = TRUE)
  SBE_plot_data[i,"Max_Canopy_Height_2013"]<-Max_Canopy_Height_2013
  
  # Progress (% of plots completed)
  print((i/dim(SBE_plot_data)[1])*100)
  
}

# Minimum DEM Value (for Height above lowest point HALP)
min(SBE_DEM) #120.74
min_DEM_SBE<- 120.74
SBE_plot_data$HALP<-SBE_plot_data$Elevation-min_DEM_SBE

# Aspect (cos transformed)
SBE_plot_data$Aspect_Cos<-cos(SBE_plot_data$Aspect)

# Changes in Max canopy height
SBE_plot_data$Change_in_MCH<-SBE_plot_data$Max_Canopy_Height_2020-SBE_plot_data$Max_Canopy_Height_2013
head(SBE_plot_data)

# Merge columns from our plot info
species_type <- species %>% rename(Plot_id = plot)
SBE_plot_data <- merge(SBE_plot_data, species_type, by = "Plot_id")
colnames(SBE_plot_data)

# Save the dataframe
save(SBE_plot_data, file ='Results_SBE.rda')

##%%%%%%%%%%%%%##
##### Danum #####
##%%%%%%%%%%%%%##

#---------------------------------------------#
### CHM height change between 2013 and 2020 ###
#---------------------------------------------#

Danum_CHM_change<-Danum_2020_CHM-Danum_2013_CHM
plot(Danum_CHM_change)

#-------------------------------------------------------------#
### Extract canopy and topographic data for each Danum plot ###
#-------------------------------------------------------------#

# Create data frame to store data
Danum_plot_data<-data.frame(matrix(dim(Danum_plots)[1],20,data=NA))
names(Danum_plot_data)<-c("Site","Plot_id","Treatment","TCH_2020","TCH_2013","TCH_change","Gap_frac_20m_2020","Gap_frac_20m_2013","Frac_Over_50m","Height_cv_2020","Height_cv_2013","Elevation","Slope","TRI","Max_Canopy_Height_2020","Max_Canopy_Height_2013","Aspect","TPI", "flow_dir", "Roughness")
Danum_plot_data$Plot_id<-Danum_plots$Plot
Danum_plot_data$Site<-"Danum"
Danum_plot_data$Treatment<- "Old Growth"

# Run loop across each plot and extract data
for (i in 1:dim(Danum_plot_data)[1]){
  
  # CHM 2020
  Danum_2020_CHM_crop<-crop(Danum_2020_CHM, Danum_plots[i,])
  Danum_plot_data[i,"TCH_2020"]<-global(Danum_2020_CHM_crop, fun='mean', na.rm=TRUE)
  Danum_plot_data[i,"Height_cv_2020"]<-global(Danum_2020_CHM_crop, fun='sd', na.rm=TRUE)/global(Danum_2020_CHM_crop, fun='mean', na.rm=TRUE)
  Danum_plot_data[i,"Gap_frac_20m_2020"]<-length(which(as.vector(na.omit(as.vector(Danum_2020_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(Danum_2020_CHM_crop))))
  
  # CHM 2013 Gap fraction at 20m and Coefficient of variation
  Danum_2013_CHM_crop<-crop(Danum_2013_CHM,Danum_plots[i,])
  Danum_plot_data[i,"TCH_2013"]<-global(Danum_2013_CHM_crop, fun='mean', na.rm=TRUE)
  Danum_plot_data[i,"Gap_frac_20m_2013"]<-length(which(as.vector(na.omit(as.vector(Danum_2013_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(Danum_2013_CHM_crop))))
  Danum_plot_data[i,"Height_cv_2013"]<-global(Danum_2013_CHM_crop, fun='sd', na.rm=TRUE)/global(Danum_2013_CHM_crop, fun='mean', na.rm=TRUE)
  # % pixels above 50m threshold (proxy for logging intensity)
  Danum_plot_data[i,"Frac_Over_50m"]<-length(which(as.vector(na.omit(as.vector(Danum_2013_CHM_crop)))>50))/length(as.vector(na.omit(as.vector(Danum_2013_CHM_crop)))) 
  
  # Total Canopy Height change
  Danum_CHM_change_crop<-crop(Danum_CHM_change,Danum_plots[i,])
  Danum_plot_data[i,"TCH_change"]<-global(Danum_CHM_change_crop, fun='mean', na.rm=TRUE)
  
  # Elevation
  Danum_DEM_crop<-crop(Danum_DEM,Danum_plots[i,])
  Danum_plot_data[i,"Elevation"]<-global(Danum_DEM_crop, fun='mean', na.rm=TRUE)
  
  # Aspect
  Danum_Aspect_crop<-crop(Danum_Aspect,Danum_plots[i,])
  Danum_plot_data[i,"Aspect"]<-global(Danum_Aspect_crop, fun='mean', na.rm=TRUE)
  
  # TPI
  Danum_TPI_crop<-crop(Danum_TPI,Danum_plots[i,])
  Danum_plot_data[i,"TPI"]<-global(Danum_TPI_crop, fun='mean', na.rm=TRUE)
  
  # Slope
  Danum_Slope_crop<-crop(Danum_Slope,Danum_plots[i,])
  Danum_plot_data[i,"Slope"]<-global(Danum_Slope_crop, fun='mean', na.rm=TRUE)
  
  # TRI
  Danum_TRI_crop<-crop(Danum_TRI,Danum_plots[i,])
  Danum_plot_data[i,"TRI"]<-global(Danum_TRI_crop, fun='mean', na.rm=TRUE)
  
  # Roughness
  Danum_Roughness_crop<-crop(Danum_Roughness,Danum_plots[i,])
  Danum_plot_data[i,"Roughness"]<-global(Danum_Roughness_crop, fun='mean', na.rm=TRUE)
  
  # flow_dir
  Danum_flow_dir_crop<-crop(Danum_flow_dir,Danum_plots[i,])
  Danum_plot_data[i,"flow_dir"]<-global(Danum_flow_dir_crop, fun='mean', na.rm=TRUE)
  
  # Max Canopy Height 2020
  MCH_danum_crop_2020<-crop(Danum_2020_CHM,Danum_plots[i,])
  values_MCH_danum_crop_2020 <- values(MCH_danum_crop_2020, na.rm = TRUE)
  Max_Canopy_Height_2020_danum<-terra::quantile(values_MCH_danum_crop_2020, probs=0.99, na.rm = TRUE)
  Danum_plot_data[i,"Max_Canopy_Height_2020"]<-Max_Canopy_Height_2020_danum
  
  # Max Canopy Height 2013
  MCH_danum_crop_2013<-crop(Danum_2013_CHM,Danum_plots[i,])
  values_MCH_danum_crop_2013 <- values(MCH_danum_crop_2013, na.rm = TRUE)
  Max_Canopy_Height_2013_danum<-quantile(values_MCH_danum_crop_2013, probs=0.99, na.rm = TRUE)
  Danum_plot_data[i,"Max_Canopy_Height_2013"]<-Max_Canopy_Height_2013_danum
  
  # Progress (% of plots completed)
  print((i/dim(Danum_plot_data)[1])*100)
  
}

# Minimum DEM Value (for Height above lowest point HALP)
min(Danum_DEM) #198.10
min_DEM_Danum<- 198.10

# HALP
Danum_plot_data$HALP<-Danum_plot_data$Elevation-min_DEM_Danum

# Aspect (cos transformed) -check which transformation to use
Danum_plot_data$Aspect_Cos<-cos(Danum_plot_data$Aspect)

# Changes in Max canopy height
Danum_plot_data$Change_in_MCH<-Danum_plot_data$Max_Canopy_Height_2020-Danum_plot_data$Max_Canopy_Height_2013
colnames(Danum_plot_data)

# Save results
save(Danum_plot_data, file ='Results_Danum.rda')
