##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SBE/Danum LiDAR Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load libraries
library(terra)
library(sf)
library(MultiscaleDTM)
library(dplyr)
library(tidyr)

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

# Load results
SBE_gap_results<-read.csv("SBE_gap_results.csv")
Danum_gap_results<-read.csv("Danum_gap_results.csv")

load('Results_SBE_update.rda')
load('Results_Danum_update.rda')

# Load shapefile and covert to sf object
Danum_plots <- vect("Danum_plots.shp")
Danum_plots_sf<- st_as_sf(Danum_plots)
# Reproject SBE_plots to match
Danum_plots <- project(Danum_plots,"epsg:32650")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### COMBINE THE DATAFRAMES ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%##
##### GAP DATA #####
##%%%%%%%%%%%%%%%%##

# Add a column specifying the site
SBE_gap_results$Site<- "SBE"
Danum_gap_results$Site<- "Danum"

# Check each dataset has the same variables
colnames(SBE_gap_results)
colnames(Danum_gap_results)

sapply(Danum_gap_results, function(x) sum(is.na(x)))
sapply(SBE_gap_results, function(x) sum(is.na(x)))

SBE_gap_results <- SBE_gap_results %>%
  pivot_wider(
    id_cols = c(Plot_id, Site),
    names_from = overlap,
    values_from = c(area_a_gaps, area_b_gaps, area_change, CHM1_mean, CHM2_mean, gap_volume_a_gaps, gap_volume_b_gaps, gap_volume_change, height_growth),
    names_glue = "{overlap}_{.value}"
  )

Danum_gap_results <- Danum_gap_results %>%
  pivot_wider(
    id_cols = c(Plot_id, Site),
    names_from = overlap,
    values_from = c(area_a_gaps, area_b_gaps, area_change, CHM1_mean, CHM2_mean, gap_volume_a_gaps, gap_volume_b_gaps, gap_volume_change, height_growth),
    names_glue = "{overlap}_{.value}"
  )

# Some persistant gaps have NAs in Danum 
sapply(Danum_gap_results, function(x) sum(is.na(x)))
sapply(SBE_gap_results, function(x) sum(is.na(x)))

# NAs are being introduced here when there was no persistant gaps
# Some persistant gaps have NAs that need to be changed to 0s
Danum_gap_results <- Danum_gap_results %>%
  mutate_all(~replace_na(., 0))

# Verify that there are no NAs left
sapply(Danum_gap_results, function(x) sum(is.na(x)))

# Combine
Gap_Data<-rbind(SBE_gap_results, Danum_gap_results)
colnames(Gap_Data)
head(as.data.frame(Gap_Data))

# Because recovered gap area in a is a number and b is zero, onceptually we need to switch these around 
# New gaps didnt exist in 2013 and recovered gaps didnt exist in 2020

# Clean data
Gap_Data$`Intact canopy_area_a_gaps`<-NULL
Gap_Data$`New gap_area_a_gaps`<-NULL
Gap_Data$`Intact canopy_area_b_gaps`<-NULL
Gap_Data$`Recovered gap_area_b_gaps`<-NULL
Gap_Data$`Intact canopy_area_change`<-NULL
Gap_Data$`New gap_area_change`<-NULL #same as new_gap_area_b_gaps
Gap_Data$`Persistent gap_area_change`<-NULL
Gap_Data$`Recovered gap_area_change`<-NULL # same as recovered_gap_area_a_gaps, but we dont want it inverse
Gap_Data$`Intact canopy_gap_volume_a_gaps`<-NULL
Gap_Data$`New gap_gap_volume_a_gaps`<-NULL
Gap_Data$`Recovered gap_gap_volume_b_gaps`<-NULL
Gap_Data$`Recovered gap_gap_volume_change`<-NULL
Gap_Data$`Intact canopy_gap_volume_b_gaps`<-NULL
Gap_Data$`Intact canopy_gap_volume_change`<-NULL
Gap_Data$`New gap_gap_volume_change`<-NULL

names(Gap_Data)[names(Gap_Data) == "Persistent gap_area_a_gaps"] <- "Persistant_gap_area_2013"
names(Gap_Data)[names(Gap_Data) == "Recovered gap_area_a_gaps"] <- "Recovered_gap_area_2013"
names(Gap_Data)[names(Gap_Data) == "New gap_area_b_gaps"] <- "New_gap_area_2020"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_area_b_gaps"] <- "Persistant_gap_area_2020"
names(Gap_Data)[names(Gap_Data) == "Intact canopy_CHM1_mean"] <- "Intact_canopy_TCH_2013"
names(Gap_Data)[names(Gap_Data) == "Intact canopy_CHM2_mean"] <- "Intact_canopy_TCH_2020"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_CHM1_mean"] <- "Persistant_gap_TCH_2013"
names(Gap_Data)[names(Gap_Data) == "Recovered gap_CHM1_mean"] <- "Recovered_gap_TCH_2013"
names(Gap_Data)[names(Gap_Data) == "New gap_CHM1_mean"] <- "New_gap_TCH_2013"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_CHM2_mean"] <- "Persistant_gap_TCH_2020"
names(Gap_Data)[names(Gap_Data) == "Recovered gap_CHM2_mean"] <- "Recovered_gap_TCH_2020"
names(Gap_Data)[names(Gap_Data) == "New gap_CHM2_mean"] <- "New_gap_TCH_2020"
names(Gap_Data)[names(Gap_Data) == "Recovered gap_gap_volume_a_gaps"] <- "Recovered_gap_volume_2013"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_gap_volume_a_gaps"] <- "Persistant_gap_volume_2013"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_gap_volume_b_gaps"] <- "Persistant_gap_volume_2020"
names(Gap_Data)[names(Gap_Data) == "New gap_gap_volume_b_gaps"] <- "New_gap_volume_2020"
names(Gap_Data)[names(Gap_Data) == "Persistent gap_gap_volume_change"] <- "Persistant_gap_volume_change"

# Create 'Existing gaps, a combined gap of recovered and persistent gaps
Gap_Data$Existing_gap_area_2013<-Gap_Data$Persistant_gap_area_2013+Gap_Data$Recovered_gap_area_2013
Gap_Data$Existing_gap_area_2020<-Gap_Data$Persistant_gap_area_2020

Gap_Data$Existing_gap_volume_2013<-Gap_Data$Persistant_gap_volume_2013+Gap_Data$Recovered_gap_volume_2013
Gap_Data$Existing_gap_volume_2020<-Gap_Data$Persistant_gap_volume_2020

Gap_Data$Existing_proportion_closed <- (Gap_Data$Existing_gap_volume_2013 - Gap_Data$Existing_gap_volume_2020) / Gap_Data$Existing_gap_volume_2013
Gap_Data$Existing_proportion_closed_perc <- (Gap_Data$Existing_proportion_closed)*100

# Will not be able to get a TCH measurement for this as we cannot combine the areas and get that from that after they have been calculated. 

# Save
save(Gap_Data, file ='Gap_Data.rda')

# Export data frame
write.csv(Gap_Data,"Gap_Data.csv",row.names = F)

##%%%%%%%%%%%%%%%%%##
#####  ACD DATA #####
##%%%%%%%%%%%%%%%%%##

# Check each dataset has the same variables
colnames(SBE_plot_data)
colnames(Danum_plot_data)

# Get the names of variables that are in df1 but not in df2
variables_only_in_df1 <- setdiff(names(SBE_plot_data), names(Danum_plot_data))

# Get the names of variables that are in df2 but not in df1
variables_only_in_df2 <- setdiff(names(Danum_plot_data), names(SBE_plot_data))

# Get the names of variables that are common in both data frames
common_variables <- intersect(names(SBE_plot_data), names(Danum_plot_data))

# Print the results
print("Variables only in df1:")
print(variables_only_in_df1)

print("Variables only in df2:")
print(variables_only_in_df2)

print("Common variables:")
print(common_variables)

# Add necessary columns
Danum_plot_data$climber_cu<- "No"
Danum_plot_data$geom<-Danum_plots_sf$geometry
Danum_plot_data$block<-'NA'
Danum_plot_data$species_ri<-'NA'
Danum_plot_data$species_co<-'NA'
Danum_plot_data$generic_di<-'NA'
Danum_plot_data$canopy_typ<-'NA'
Danum_plot_data$canopy_str<-'NA'
Danum_plot_data$random<-'NA'
Danum_plot_data$planted<-'Old Growth'
Danum_plot_data$treatment<-'Old Growth'

# Rename SBE 
names(SBE_plot_data)[names(SBE_plot_data) == "Site.x"] <- "site"
names(SBE_plot_data)[names(SBE_plot_data) == "block.x"] <- "block"
names(SBE_plot_data)[names(SBE_plot_data) == "species_ri.x"] <- "species_ri"
names(SBE_plot_data)[names(SBE_plot_data) == "species_co.x"] <- "species_co"
names(SBE_plot_data)[names(SBE_plot_data) == "climber_cu.x"] <- "climber_cu"
names(SBE_plot_data)[names(SBE_plot_data) == "generic_di.x"] <- "generic_di"
names(SBE_plot_data)[names(SBE_plot_data) == "canopy_typ.x"] <- "canopy_typ"
names(SBE_plot_data)[names(SBE_plot_data) == "canopy_str.x"] <- "canopy_str"
names(SBE_plot_data)[names(SBE_plot_data) == "treatment.x"] <- "treatment"
names(SBE_plot_data)[names(SBE_plot_data) == "random.x"] <- "random"
names(SBE_plot_data)[names(SBE_plot_data) == "geom.x"] <- "geom"
names(SBE_plot_data)[names(SBE_plot_data) == "planted.x"] <- "planted"

SBE_plot_data$block.y <- NULL
SBE_plot_data$species_ri.y<- NULL
SBE_plot_data$species_co.y<- NULL
SBE_plot_data$generic_di.y<- NULL
SBE_plot_data$canopy_typ.y<- NULL
SBE_plot_data$canopy_str.y<- NULL
SBE_plot_data$planted.y<- NULL
SBE_plot_data$treatment.y<- NULL
SBE_plot_data$random.y<- NULL
SBE_plot_data$geom.y<- NULL
SBE_plot_data$climber_cu.y <- NULL

# Combine
All_Data<-rbind(SBE_plot_data,Danum_plot_data)

# Save
colnames(All_Data)
save(All_Data, file ='All_Data.rda')

# Export data frame
write.csv(All_Data,"All_Data.csv",row.names = F)

##%%%%%%%%%%%%%%%%%%##
##### FINAL DATA #####
##%%%%%%%%%%%%%%%%%%##

colnames(Gap_Data)
colnames(All_Data)

# Only want gap volumes now, not area 
Gap_Data$Persistant_gap_area_2013<-NULL
Gap_Data$Recovered_gap_area_2020<-NULL
Gap_Data$Persistant_gap_area_2020<-NULL
Gap_Data$Existing_gap_area_2013<-NULL
Gap_Data$Existing_gap_area_2020<-NULL
Gap_Data$New_gap_area_2020<-NULL
Gap_Data$Recovered_gap_area_2013<-NULL

# Using 'existing gaps' now, not persitent and recovered or intact
Gap_Data$Persistant_gap_TCH_2013<-NULL
Gap_Data$Persistant_gap_TCH_2020<-NULL
Gap_Data$Persistant_gap_volume_2013<-NULL
Gap_Data$Persistant_gap_volume_2020<-NULL
Gap_Data$Persistant_gap_volume_change<-NULL
Gap_Data$`Persistent gap_height_growth`<-NULL

Gap_Data$Recovered_gap_TCH_2013<-NULL
Gap_Data$Recovered_gap_TCH_2020<-NULL
Gap_Data$Recovered_gap_volume_2020<-NULL
Gap_Data$`Recovered gap_height_growth`<-NULL

Gap_Data$`Intact canopy_height_growth`<-NULL
Gap_Data$Intact_canopy_TCH_2013<-NULL
Gap_Data$Intact_canopy_TCH_2020<-NULL

# cant get tch for existing as it is worked out after means from chm are calculated for persist and recovered gaps
Gap_Data$New_gap_TCH_2013<-NULL
Gap_Data$New_gap_TCH_2020<-NULL
Gap_Data$`New gap_height_growth`<-NULL

#Just keep the perc version:
Gap_Data$Existing_proportion_closed<-NULL

# /7 so it is per year
Gap_Data$Existing_proportion_closed_perc_yr<- (Gap_Data$Existing_proportion_closed_perc)/7
Gap_Data$New_gap_volume_yr<-(Gap_Data$New_gap_volume_2020)/7

# Merge ACD and Gap Data
Final_Data <- merge(Gap_Data, All_Data, by = 'Plot_id')

# Clean up data frame
colnames(Final_Data)
names(Final_Data)[names(Final_Data) == "Site.x"] <- "Site"
names(Final_Data)[names(Final_Data) == "flow_dir"] <- "Flow_dir"
names(Final_Data)[names(Final_Data) == "block"] <- "Block"
names(Final_Data)[names(Final_Data) == "species_ri"] <- "Species_richness"
names(Final_Data)[names(Final_Data) == "species_co"] <- "Species_composition"
names(Final_Data)[names(Final_Data) == "climber_cu"] <- "Climbers_cut"
names(Final_Data)[names(Final_Data) == "generic_di"] <- "Generic_diversity"
names(Final_Data)[names(Final_Data) == "canopy_typ"] <- "Canopy_type"
names(Final_Data)[names(Final_Data) == "geom"] <- "Geometry"
names(Final_Data)[names(Final_Data) == "CC_change"] <- "Cover20_change"
names(Final_Data)[names(Final_Data) == "planted"] <- "Planted"
names(Final_Data)[names(Final_Data) == "random"] <- "Random"

Final_Data$canopy_str <- NULL
Final_Data$Site.y <- NULL 
Final_Data$treatment <- NULL

Final_Data$Cover_20_2013_perc<-(Final_Data$Cover_20_2013)*100 
Final_Data$Cover_20_2020_perc<-(Final_Data$Cover_20_2020)*100 
Final_Data$Frac_Over_50m_perc<-(Final_Data$Frac_Over_50m)*100 

# To avoid pseudoreplication ensure the ID is only repeated once per dataset tested
head(Final_Data)
colnames(Final_Data)

# Identify rows with no missing values in the relevant columns
colSums(is.na(Final_Data))

# Save
save(Final_Data, file ='Final_Data.rda')

# Export data frame
write.csv(Final_Data,"Final_Data.csv",row.names = F)
