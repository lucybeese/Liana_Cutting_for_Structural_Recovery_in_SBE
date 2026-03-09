##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SBE/Danum LiDAR Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load libraries
library(terra)
library(sf)
library(MultiscaleDTM)
library(dplyr)
library(tidyr)
library(TeachingDemos)
library(scales)
library(emmeans)
library(AICcmodavg)
library(corrplot)

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

load('Final_DataV2.rda')

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

##%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SUPPLEMENTARY STATS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%##

SBE<- subset(Final_Data, Site== "SBE")
Danum<- subset(Final_Data, Site=="Danum")

mean(SBE$ACD_2013)
mean(Danum$ACD_2013)

Planted <- subset(Final_Data, combined== "Planted")
Liana <- subset(Final_Data, combined== "Planted with lianas removed")
Control <- subset(Final_Data, combined== "Control")
OldG <- subset(Final_Data, combined== "Old Growth")

mean(Planted$ACD_change)
mean(Liana$ACD_change)
# 18.8 times higher mean ACD change in Liana than Planted

mean(Planted$New_gap_volume_yr)
mean(Liana$New_gap_volume_yr)
# 2.17 times more new gaps in planted than liana removed 

mean(Planted$Existing_proportion_closed_perc_yr)
mean(Liana$Existing_proportion_closed_perc_yr)
# 1.31 times more closure in liana plots

#-----------------------------------------------#
## Gap frac at 10m and cover 20/10 correlation ##
#-----------------------------------------------#

# Create data frame to store data
SBE_10<-data.frame(matrix(dim(SBE_plots)[1],6,data=NA))
names(SBE_10)<-c("Site","Plot_id","Gap_frac_10m_2020","Gap_frac_10m_2013", "Gap_frac_20m_2020","Gap_frac_20m_2013")
SBE_10$Site<-"SBE"
SBE_10$Plot_id<-SBE_plots$Plot_id

# Run loop across each plot and extract data
for (i in 1:dim(SBE_10)[1]){
  SBE_2020_CHM_crop<-crop(SBE_2020_CHM, SBE_plots[i,])
  SBE_2013_CHM_crop<-crop(SBE_2013_CHM, SBE_plots[i,])
  SBE_10[i,"Gap_frac_10m_2013"]<-length(which(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))<10))/length(as.vector(na.omit(as.vector(SBE_2013_CHM_crop))))
  SBE_10[i,"Gap_frac_10m_2020"]<-length(which(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))<10))/length(as.vector(na.omit(as.vector(SBE_2013_CHM_crop))))
  SBE_10[i,"Gap_frac_20m_2020"]<-length(which(as.vector(na.omit(as.vector(SBE_2020_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(SBE_2020_CHM_crop)))) 
  SBE_10[i,"Gap_frac_20m_2013"]<-length(which(as.vector(na.omit(as.vector(SBE_2013_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(SBE_2013_CHM_crop))))
}

# Create data frame to store data
Danum_10<-data.frame(matrix(dim(Danum_plots)[1],6,data=NA))
names(Danum_10)<-c("Site","Plot_id","Gap_frac_10m_2020","Gap_frac_10m_2013", "Gap_frac_20m_2020","Gap_frac_20m_2013")
Danum_10$Plot_id<-Danum_plots$Plot
Danum_10$Site<-"Danum"

# Run loop across each plot and extract data
for (i in 1:dim(Danum_10)[1]){
  
  # CHM 2020
  Danum_2020_CHM_crop<-crop(Danum_2020_CHM, Danum_plots[i,])
  Danum_2013_CHM_crop<-crop(Danum_2013_CHM, Danum_plots[i,])
  Danum_10[i,"Gap_frac_20m_2020"]<-length(which(as.vector(na.omit(as.vector(Danum_2020_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(Danum_2020_CHM_crop))))
  Danum_10[i,"Gap_frac_20m_2013"]<-length(which(as.vector(na.omit(as.vector(Danum_2013_CHM_crop)))<20))/length(as.vector(na.omit(as.vector(Danum_2013_CHM_crop))))
  Danum_10[i,"Gap_frac_10m_2020"]<-length(which(as.vector(na.omit(as.vector(Danum_2020_CHM_crop)))<10))/length(as.vector(na.omit(as.vector(Danum_2020_CHM_crop))))
  Danum_10[i,"Gap_frac_10m_2013"]<-length(which(as.vector(na.omit(as.vector(Danum_2013_CHM_crop)))<10))/length(as.vector(na.omit(as.vector(Danum_2013_CHM_crop))))
  
}

data_10<-rbind(SBE_10, Danum_10)
# Canopy cover calculations at 20m for 2020 and 2013
data_10$Cover_20_2020 <- 1 - data_10$Gap_frac_20m_2020
data_10$Cover_20_2013 <- 1 - data_10$Gap_frac_20m_2013
data_10$Cover_10_2020 <- 1 - data_10$Gap_frac_10m_2020
data_10$Cover_10_2013 <- 1 - data_10$Gap_frac_10m_2013


cover_correlation_20 <- cor(data_10$Gap_frac_20m_2020, data_10$Gap_frac_10m_2020, use = "complete.obs")
cover_correlation_20

cover_correlation_13 <- cor(data_10$Gap_frac_20m_2013, data_10$Gap_frac_10m_2013, use = "complete.obs")
cover_correlation_13

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 10: LOGGING INTENSITY/ HALP CORRELATION ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# See if logging intensity is correlated with topography 
# Trees up high are often left unlogged thus they get bigger.

# Check predictors for correlation
corr_halp<-cor(SBE$HALP, SBE$Frac_Over_50m_perc ,use = "complete.obs")
print(corr_halp)
# HALP- Frac over 50 not too highly correlated: 0.25431668

frac_halp <- cor(SBE[, c("ACD_2020", "ACD_change", "Cover_20_2020_perc", "Frac_Over_50m_perc", "HALP", "Existing_proportion_closed_perc_yr", "New_gap_volume_yr")],use = "complete.obs")
print(frac_halp)

#--------------------#
###### Plotting ######
#--------------------#

# set wd to figures folder

cairo_pdf("Correlation.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,1),
    mar=c(6,8,2,2),
    las=1,
    xpd=TRUE,
    cex.axis=1.6,
    ps=14)

# Diverging colour palette
neg <- colorRampPalette(c("#2A7879", "#AEC5C6"))(100)
pos <- colorRampPalette(c("#A793B5", "#481567"))(100)
col <- c(neg, pos)

custom_labels <- c(
  "ACD in 2020",
  "ACD change",
  "Cover20 in 2020",
  "Logging intensity",
  "HALP",
  "Existing gaps closed per year",
  "New gap volume per year"
)

rownames(frac_halp) <- custom_labels
colnames(frac_halp) <- custom_labels

corrplot(frac_halp,
         type="lower",
         tl.col="black",
         tl.cex=1.6,
         tl.srt=70,
         col=col,
         cl.cex=1.3)

dev.off()
