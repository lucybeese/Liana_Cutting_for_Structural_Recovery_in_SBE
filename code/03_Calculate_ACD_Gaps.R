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

# Load results
load('Results_SBE.rda')
load('Results_Danum.rda')

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%%##

# Import CHM 2020
SBE_2020_CHM <- rast("SBE_2020_CHM_crop.tif")

# Import CHM 2013
SBE_2013_CHM <- rast("SBE_2013_CHM_crop.tif")

# Load shapefile and covert to sf object
SBE_plots <- vect("SBE_plots.shp")
# Reproject SBE_plots to match
SBE_plots <- project(SBE_plots,"epsg:32650")
SBE_plots_sf <- st_as_sf(SBE_plots)

##%%%%%%%%%%%%%##
##### DANUM #####
##%%%%%%%%%%%%%##

# Import CHM 2020
Danum_2020_CHM<-rast("Danum_2020_CHM_crop.tif")

# Import CHM 2013
Danum_2013_CHM<-rast("Danum_2013_CHM_crop.tif")

# Load shapefile and covert to sf object
Danum_plots <- vect("Danum_plots.shp")
Danum_plots_sf<- st_as_sf(Danum_plots)
# Reproject SBE_plots to match
Danum_plots <- project(Danum_plots,"epsg:32650")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 2: LOAD FUNCTIONS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Create and load updated ForestGapR Functions to work with terra

#-------------------#
### GetForestGaps ###
#-------------------#

getForestGaps <- function(chm_layer, threshold = 10, size = c(1, 10^4)) {
  chm_layer[chm_layer > threshold] <- NA
  chm_layer[chm_layer <= threshold] <- 1
  gaps <- terra::patches(chm_layer, directions = 8, allowGaps = FALSE)
  rcl <- terra::freq(gaps)
  rcl$layer<-NULL
  rcl[, 2] <- rcl[, 2] * terra::res(chm_layer)[1]^2
  z <- terra::classify(gaps, rcl = rcl, right = FALSE)
  z[is.na(gaps)] <- NA
  gaps[z > size[2]] <- NA
  gaps[z < size[1]] <- NA
  gaps <- terra::patches(gaps, directions = 8, allowGaps = FALSE)
  names(gaps) <- "gaps"
  return(gaps)
}

#------------------#
### GapChangeDec ###
#------------------#

GapChangeDec <- function(gap_layer1, gap_layer2) {
  gap_layer1[!is.na(gap_layer1)] <- 1
  gap_layer2[!is.na(gap_layer2)] <- 2
  gap_layer1[is.na(gap_layer1)] <- 0
  gap_layer2[is.na(gap_layer2)] <- 0
  gap_diff <- gap_layer2 - gap_layer1
  gap_diff[gap_diff[] != 2] <- NA
  return(gap_diff)
}

CD_summarize_gap_dynamics_short = function(a_chm, a_gaps, b_chm, b_gaps){
  
  # This function summarizes the change in each gap and disturbance based on the rasters only
  # It relies on the order of cells in the raster, and dplyr summarize. This is much faster than polygonizing
  a_gaps[a_gaps>0]=1
  b_gaps[b_gaps>0]=1
  
  tile_stack = terra::rast(list(a_chm,  b_chm,  a_gaps,  b_gaps))
  names(tile_stack) = c("a_chm","b_chm","a_gaps","b_gaps")
  
  # A stack of rasters is created, combining different layers for further analysis
  tile_tb = tile_stack %>% as.data.frame(xy=TRUE)
  
  # Define overlaps
  # This is a hierarchy. The higher numbers overwrite the others
  # These values mark pixels, not gaps. A single gap can contain all three values.
  # I only use numbers so I can save them as rasters & visualize the results
  tile_tb$overlap = "Intact canopy"
  tile_tb$overlap[tile_tb$a_gaps==1]= "Recovered gap"
  tile_tb$overlap[tile_tb$b_gaps==1]= "New gap"
  tile_tb$overlap[tile_tb$a_gaps==1 & tile_tb$b_gaps==1]="Persistent gap"
  
  # Summarize gaps
  # For each gap in year a, these are the characteristics and changes
  # Important: I group by unbuffered gaps / disturbances to preserve the ID + so that they don't overlap and double count.
  # Therefore, I need the disturbances to be buffered for the gaps, and the gaps to be buffered for the disturbances, if that makes sense.
  
  df = tile_tb %>% dplyr::ungroup() %>% dplyr::group_by(overlap) %>%
    dplyr::summarise(
      area_a_gaps = ifelse(is.na(sum(a_gaps)), 0, sum(a_gaps)),
      area_b_gaps = ifelse(is.na(sum(b_gaps)), 0, sum(b_gaps)),
      area_change = area_b_gaps - area_a_gaps,
      CHM1_mean = mean(a_chm, na.rm = TRUE),
      CHM2_mean = mean(b_chm, na.rm = TRUE),
      gap_volume_a_gaps = CHM1_mean * area_a_gaps,
      gap_volume_b_gaps = CHM2_mean * area_b_gaps,
      gap_volume_change = ifelse(area_change > 0, gap_volume_b_gaps, gap_volume_b_gaps - gap_volume_a_gaps),
      height_growth = CHM2_mean - CHM1_mean
    )
  
  return(df)
}

##%%%%%%%%%%%%%%%%%%%%%##

#### CALCULATING ACD ####

##%%%%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%##
##### SBE ######
#%%%%%%%%%%%%%##

# Canopy cover calculations at 20m for 2020 and 2013
SBE_plot_data$Cover_20_2020 <- 1 - SBE_plot_data$Gap_frac_20m_2020
SBE_plot_data$Cover_20_2013 <- 1 - SBE_plot_data$Gap_frac_20m_2013

# Cover Resid (as described in Jucker et al. 2018)
SBE_plot_data$Cover_resid_2013 <- SBE_plot_data$Cover_20_2013 - (1/(1+exp(12.431)*(SBE_plot_data$TCH_2013^-4.061)))
SBE_plot_data$Cover_resid_2020 <- SBE_plot_data$Cover_20_2020 - (1/(1+exp(12.431)*(SBE_plot_data$TCH_2020^-4.061)))

# ACD (as described in Jucker et al. 2018)
SBE_plot_data$ACD_2020<- (0.62369 * SBE_plot_data$TCH_2020 ^1.63899) * (1 + 1.983 * SBE_plot_data$Cover_resid_2020)^1.081
SBE_plot_data$ACD_2013<- (0.62369 * SBE_plot_data$TCH_2013 ^1.63899) * (1 + 1.983 * SBE_plot_data$Cover_resid_2013)^1.081

# ACD change per year
SBE_plot_data$ACD_change<-(SBE_plot_data$ACD_2020-SBE_plot_data$ACD_2013)/7

# Canopy cover change per year (%)
SBE_plot_data$CC_change<- (SBE_plot_data$Cover_20_2020-SBE_plot_data$Cover_20_2013)/7

colnames(SBE_plot_data)[colnames(SBE_plot_data) == "Climbers"] <- "Climbers_cut"

colnames(SBE_plot_data)
save(SBE_plot_data, file ='Results_SBE_update.rda')

##%%%%%%%%%%%%%##
##### Danum #####
##%%%%%%%%%%%%%##

# Canopy cover calculations at 20m for 2020 and 2013
Danum_plot_data$Cover_20_2020 <- 1 - Danum_plot_data$Gap_frac_20m_2020
Danum_plot_data$Cover_20_2013 <- 1 - Danum_plot_data$Gap_frac_20m_2013

# Cover Resid (Jucker et al. 2018)
Danum_plot_data$Cover_resid_2020 <- Danum_plot_data$Cover_20_2013 - (1/(1+exp(12.431)*(Danum_plot_data$TCH_2020^-4.061)))
Danum_plot_data$Cover_resid_2013 <- Danum_plot_data$Cover_20_2020 - (1/(1+exp(12.431)*(Danum_plot_data$TCH_2013^-4.061)))

# ACD
Danum_plot_data$ACD_2020<- (0.62369 * Danum_plot_data$TCH_2020 ^1.63899) * (1 + 1.983 * Danum_plot_data$Cover_resid_2020)^1.081
Danum_plot_data$ACD_2013<- (0.62369 * Danum_plot_data$TCH_2013 ^1.63899) * (1 + 1.983 * Danum_plot_data$Cover_resid_2013)^1.081

# ACD change per year
Danum_plot_data$ACD_change<-(Danum_plot_data$ACD_2020-Danum_plot_data$ACD_2013)/7

# CC change per year
Danum_plot_data$CC_change<- (Danum_plot_data$Cover_20_2020-Danum_plot_data$Cover_20_2013)/7

# Remove NA entries
Danum_plot_data <- Danum_plot_data[complete.cases(Danum_plot_data$Site), ]

colnames(Danum_plot_data)
save(Danum_plot_data, file ='Results_danum_update.rda')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### CALCULATING CANOPY GAPS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%##
##### SBE #####
##%%%%%%%%%%%##

# Looking at gaps @10m
threshold<- 10
size <- c(25, 10000)

# The following two lines take some time (5 mins)
SBE_gaps_2013<-getForestGaps(chm_layer=SBE_2013_CHM, threshold=threshold, size=size) 
SBE_gaps_2020<-getForestGaps(chm_layer=SBE_2020_CHM, threshold=threshold, size=size)

# Crop to equal sizes
SBE_gaps_2013<- crop(SBE_gaps_2013, SBE_gaps_2020)
SBE_gaps_2020<- crop(SBE_gaps_2020, SBE_gaps_2013)

# Write gaps
writeRaster(SBE_gaps_2013, "SBE_gaps_2013.tif", overwrite=T)
writeRaster(SBE_gaps_2020, "SBE_gaps_2020.tif", overwrite=T)

# Gap change detection
SBE_gap_change <- GapChangeDec(SBE_gaps_2013, SBE_gaps_2020)
plot(SBE_gap_change)

# Plotting ALS-derived CHM and forest gaps
par(mfrow = c(1, 3))

plot(SBE_2013_CHM, main="Forest Canopy Gap - 2013")
plot(SBE_gaps_2013, add=T, col="red", legend=FALSE)

plot(SBE_2020_CHM,  main="Forest Canopy Gap - 2020")
plot(SBE_gaps_2020, add=T,col="red", legend=FALSE)

plot(SBE_2020_CHM,main="Forest Gaps Changes Detection")
plot(SBE_gap_change, add=T, col="yellow", legend=FALSE)

#------------------------------------#
### Get Gap Dynamics for each plot ###
#------------------------------------#

# Create an empty list to store the results
SBE_results_list <- list()

# Loop through each row in SBE_plots_sf
for (i in 1:nrow(SBE_plots_sf)) {
  
  # Crop CHM to the extent of the plot for both years
  # SBE contains 124 plots so this will find the gap dynamics for each plot
  SBE_2013_CHM_crop <- crop(SBE_2013_CHM, (SBE_plots[i,]))
  SBE_2020_CHM_crop <- crop(SBE_2020_CHM, (SBE_plots[i,]))
  
  # Crop gaps to the extent of the plot for both years
  SBE_gaps_2013_crop <- crop(SBE_gaps_2013, (SBE_plots[i,]))
  SBE_gaps_2020_crop <- crop(SBE_gaps_2020, (SBE_plots[i,]))
  
  # Calculate gap dynamics for the current plot
  summary_df <- CD_summarize_gap_dynamics_short(SBE_2013_CHM_crop, SBE_gaps_2013_crop, SBE_2020_CHM_crop, SBE_gaps_2020_crop)
  
  # Add 'plot' column from SBE_plots_sf to the resulting dataframe
  summary_df$Plot_id <- SBE_plots_sf$Plot_id[i]
  
  # Store the result in the list
  SBE_results_list[[i]] <- summary_df
}

# Combine all data frames into one and bind them row-wise
SBE_gap_results <- bind_rows(SBE_results_list)

# Because SBE is 4ha while Danum is 1ha plots, we need to divide volumes and areas by 4.
SBE_gap_results$area_a_gaps<-(SBE_gap_results$area_a_gaps)/4
SBE_gap_results$area_b_gaps<-(SBE_gap_results$area_b_gaps)/4
SBE_gap_results$area_change<-(SBE_gap_results$area_b_gaps)-(SBE_gap_results$area_a_gaps)
SBE_gap_results$gap_volume_a_gaps<-(SBE_gap_results$gap_volume_a_gaps)/4
SBE_gap_results$gap_volume_b_gaps<-(SBE_gap_results$gap_volume_b_gaps)/4
SBE_gap_results$gap_volume_change<-(SBE_gap_results$gap_volume_b_gaps)-(SBE_gap_results$gap_volume_a_gaps)

# If a gap spans two plots/ extends beyond the plot boundary, it's split into a gap for each plot

# Save data
colnames(SBE_gap_results)
save(SBE_gap_results, file ='SBE_gap_results.rda')

# Export data frame
write.csv(SBE_gap_results,"SBE_gap_results.csv",row.names = F)

##%%%%%%%%%%%%%##
##### DANUM #####
##%%%%%%%%%%%%%##

# Retrieve gaps
Danum_gaps_2013<-getForestGaps(chm_layer=Danum_2013_CHM, threshold=threshold, size=size)
Danum_gaps_2020<-getForestGaps(chm_layer=Danum_2020_CHM, threshold=threshold, size=size)

# Crop to equal sizes
Danum_gaps_2013<- crop(Danum_gaps_2013, Danum_gaps_2020)
Danum_gaps_2020<- crop(Danum_gaps_2020, Danum_gaps_2013)

# Write gaps
writeRaster(Danum_gaps_2013, "Danum_gaps_2013.tif", overwrite=T)
writeRaster(Danum_gaps_2020, "Danum_gaps_2020.tif", overwrite=T)

Danum_gap_change<- GapChangeDec(Danum_gaps_2013, Danum_gaps_2020)
par(mfrow= c(1,1))
plot(Danum_gap_change)

# Plotting ALS-derived CHM and forest gaps
par(mfrow = c(1, 3))

plot(Danum_2013_CHM, main="Forest Canopy Gap - 2013")
plot(Danum_gaps_2013, add=TRUE, col="red", legend=FALSE)

plot(Danum_2020_CHM,  main="Forest Canopy Gap - 2020")
plot(Danum_gaps_2020, add=TRUE,col="red", legend=FALSE)

plot(Danum_2020_CHM,main="Forest Gaps Changes Detection")
plot(Danum_gap_change, add=TRUE, col="yellow", legend=FALSE)

#------------------------------------#
### Get Gap Dynamics for each plot ###
#------------------------------------#

# Create an empty list to store the results
Danum_results_list <- list()

# Loop through each row in SBE_plots_sf
for (i in 1:nrow(Danum_plots)) {
  
  ## Crop CHM to the extent of the plot for both years
  Danum_2013_CHM_crop <- crop(Danum_2013_CHM , (Danum_plots[i,]))
  Danum_2020_CHM_crop <- crop(Danum_2020_CHM , (Danum_plots[i,]))
  
  ## Crop gaps to the extent of the plot for both years
  Danum_gaps_2013_crop <- crop(Danum_gaps_2013 , (Danum_plots[i,]))
  Danum_gaps_2020_crop <- crop(Danum_gaps_2020 , (Danum_plots[i,]))
  
  # Calculate gap dynamics for the current plot
  summary_df <- CD_summarize_gap_dynamics_short(Danum_2013_CHM_crop, Danum_gaps_2013_crop, Danum_2020_CHM_crop, Danum_gaps_2020_crop)
  
  # Add 'plot' column from SBE_plots_sf to the resulting dataframe
  summary_df$Plot_id <- Danum_plot_data$Plot_id[i]
  
  # Store the result in the list
  Danum_results_list[[i]] <- summary_df
}

# Combine all data frames into one and bind them row-wise
Danum_gap_results <- bind_rows(Danum_results_list)

# Save data
colnames(Danum_gap_results)
save(Danum_gap_results, file ='Danum_gap_results.rda')

sapply(Danum_gap_results, function(x) sum(is.na(x)))

# Export data frame
write.csv(Danum_gap_results,"Danum_gap_results.csv",row.names = F)
