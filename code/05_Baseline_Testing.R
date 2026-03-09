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

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

load('Final_Data.rda')
head(Final_Data)

##%%%%%%%%%%%%%%%%%%%%%%##

#### BASELINE TESTING ####

##%%%%%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

##### Logging intensity (Cover50) #####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Investigating if logging intensity differed between SBE plots

# Create a combined variable of Planting and Climber cutting - do not include old growth.
# Create a new dataframe with rows where 'Site' is 'SBE'
SBE <- subset(Final_Data, Site == "SBE")
Danum <-subset(Final_Data, Site == "Danum")

SBE$combined <- with(SBE, interaction(Climbers_cut, Planted), drop = TRUE)
mylabels <- c("Control", "Planted", "Planted with lianas removed")
SBE$combined  <- factor(SBE$combined, levels = c("No.Control", "No.Planted", "Yes.Planted"), labels = mylabels)
levels(SBE$combined)

#----------------------------------------------------------------------------------#
### Investigating if logging intensity systematically differs between treatments ###
#----------------------------------------------------------------------------------#

# ANOVA
logging_intensity_treatments<-aov(SBE$Frac_Over_50m~SBE$combined) 
summary(logging_intensity_treatments)
# Overall no statistical differences in logging intensity between the treatments at baseline: 
# SBE$combined   2 0.00048 0.0002398   0.842  0.433

# Check test assumptions
# plot(logging_intensity_treatments) # QQ plot not met
# residual vs fitted not met
logging_intensity_treatments_residuals<-residuals(logging_intensity_treatments)
hist(logging_intensity_treatments_residuals) # Left-skewed residuals

# Assumptions of lm don't seem to be met here
# Tried a log transformation and they still weren't met

# Try a Kruskal-Wallis test
KW_logging_intensity_treatments <- kruskal.test(SBE$Frac_Over_50m~SBE$combined)
KW_logging_intensity_treatments
# p-value = 0.3462

#--------------------#
###### Plotting ######
#--------------------#

#///////////////////////#
### Logging_intensity ###
#///////////////////////#

# Set wd to figures folder

# Load violin plot function
source(file="wvioplot.r")

# Plot figure
cairo_pdf("Logging_intensity.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,1),mar=c(6,6,2,2),las=1,xpd=T)

min(SBE$Frac_Over_50m_perc) # 0
max(SBE$Frac_Over_50m_perc) # 10.78119

plot(1,bty="l",xlab="",ylab="Cover50 (%)", type="n",yaxt="n",ylim=c(0,15),xlim=c(0.5,3.3),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(SBE$Frac_Over_50m_perc[SBE$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(SBE$Frac_Over_50m_perc[SBE$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(SBE$Frac_Over_50m_perc[SBE$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)

# these lines get rid of extra lines, make them the same colour as background:
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4) 
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,15,2.5),cex.axis=2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-3)
axis(1, at = seq(1, 3.5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas removed"),
      at = seq(1, 3, 1), padj = 1.5, cex = 1.8, col = "grey20")

dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	Baseline differences in ACD/CC between Danum/SBE #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#$$$$$$$$$$$$$$#
##### 2013 #####
#$$$$$$$$$$$$$$#

#///////////////#
###### ACD ######
#///////////////#

#/////////////////////////////////////////////////#
### Baseline differences in absolute ACD (2013) ###
#/////////////////////////////////////////////////#

# lm 
baseline_ACD_2013<-lm(Final_Data$ACD_2013~Final_Data$Site) 
summary(baseline_ACD_2013)
# Overall statistical differences between the sites at baseline: 
# Final_Data$SiteSBE  -20.131      5.349  -3.763 0.000232 ***

# Danum has 20.13 Mg C/ha more than SBE in 2013
# This is expected because it is likely that the logged site would have 
# less carbon than the Old Growth site at the start of the experiment

# Check test assumptions
plot(baseline_ACD_2013) # QQ plot okay
# residual vs fitted okay
baseline_ACD_2013_residuals<-residuals(baseline_ACD_2013)
hist(baseline_ACD_2013_residuals) # Normally distributed residuals

# Assumptions of lm are met

#///////////////////#
###### Cover20 ######
#///////////////////#

#////////////////////////////////////////////#
### Baseline differences in Cover20 (2013) ###
#////////////////////////////////////////////#

# lm
baseline_cover20<-lm(Final_Data$Cover_20_2013_perc~Final_Data$Site) 
summary(baseline_cover20)
# significant baseline differences :
# Final_Data$SiteSBE   -3.448      1.422  -2.425   0.0164 *
# 3% Greater cov20 in Old growth than SBE as expected

# Check test assumptions
plot(baseline_cover20) # QQ plot okay

# residual vs fitted okay
baseline_cover20_residuals<-residuals(baseline_cover20)
hist(baseline_cover20_residuals) # residuals okay

# Assumptions of lm are met

#//////////////////////////////////////////////////#
### Baseline differences in Canopy height (2013) ###
#//////////////////////////////////////////////////#

# lm controlling for planting and logging intensity
baseline_TCH_2013<-lm(Final_Data$TCH_2013~Final_Data$Site) 
summary(baseline_TCH_2013)
# significant differences in gap volume change:
# Final_Data$SiteSBE  -7.0454     0.5139  -13.71   <2e-16 ***
#7.05m  Greater canopy height in Old growth than SBE as expect due to no logging

# Check test assumptions
plot(baseline_TCH_2013) # QQ okay

# residual vs fitted okay
TCH_2013_residuals<-residuals(baseline_TCH_2013)
hist(TCH_2013_residuals) # residuals okay

#$$$$$$$$$$$$$$$$#
##### Change #####
#$$$$$$$$$$$$$$$$#

#///////////////#
###### ACD ######
#///////////////#

#////////////////////////////////////////#
### Baseline differences in ACD change ###
#////////////////////////////////////////#

# Set Site as a factor
Final_Data$Site<-factor(Final_Data$Site)

# lm 
baseline_ACD_change<-lm(Final_Data$ACD_change~Final_Data$Site) 
summary(baseline_ACD_change)

# Overall no statistical differences between the sites at baseline: 
# Final_Data$SiteSBE -0.02245    0.31006  -0.072   0.9424  

# Check test assumptions of lm
plot(baseline_ACD_change) # QQ okay
# residual vs fitted okay
baseline_ACD_change_residuals<-residuals(baseline_ACD_change)
hist(baseline_ACD_change_residuals) # Normally distributed residuals 

# All assumptions of lm are met

#////////////////#
###### Gaps ######
#////////////////#

#-----------------------------------------------#
### Baseline differences in gap volume change ###
#-----------------------------------------------#

#%%%%%%%%%#
### NEW ###
#%%%%%%%%%#

# lm
baseline_newGV<-lm(Final_Data$New_gap_volume_yr~Final_Data$Site) 
summary(baseline_newGV)
# Overall statistical differences between the sites at baseline: 
#base_New$SiteSBE    82.08      31.54   2.603   0.0101 *  

# SBE has 82.08m3/ha/yr more new gaps form over the experiment than Danum

# Check test assumptions
plot(baseline_newGV) # QQ okay
# residual vs fitted okay
baseline_newGV_residuals<-residuals(baseline_newGV)
hist(baseline_newGV_residuals) # normally distributed residuals (slightly right skewed perhaps but okay)

# Assumptions of lm seem to be met

#%%%%%%%%%%%%%%#
### Existing ###
#%%%%%%%%%%%%%%#

# lm
baseline_ExistingGV<-lm(Final_Data$Existing_proportion_closed_perc_yr~Final_Data$Site) 
summary(baseline_ExistingGV)

# Significant baseline differences:
# Final_Data$SiteSBE  -1.3122     0.3621  -3.624 0.000384 ***
# Gaps at Danum closed more than SBE

# Check test assumptions
plot(baseline_ExistingGV) # QQ plot ok
# residual vs fitted okay
baseline_ExistingGV_residuals<-residuals(baseline_ExistingGV)
hist(baseline_ExistingGV_residuals) # Normal residuals (a bit right skewed but okay)

# Assumptions of lm are met

#/////////////////////////#
###### Canopy Height ######
#/////////////////////////#

#------------------------------------------------#
## Baseline differences in Canopy height change ##
#------------------------------------------------#

# lm
Final_Data$TCH_change_yr<-(Final_Data$TCH_change)/7
baseline_TCH_change_yr<-lm(Final_Data$TCH_change_yr~Final_Data$Site) 
summary(baseline_TCH_change_yr)
# significant baseline differences in gap volume change:
# Final_Data$SiteSBE  0.14732    0.02637   5.586 9.24e-08 ***

# Check test assumptions
plot(baseline_TCH_2013) # QQ plot okay
# residual vs fitted looks okay
TCH_change_yr_residuals<-residuals(baseline_TCH_change_yr)
hist(TCH_2013_residuals) # residuals okay

# Greater canopy height change in SBE than old growth

#$$$$$$$$$$$$$$#
##### 2020 #####
#$$$$$$$$$$$$$$#

#///////////////#
###### ACD ######
#///////////////#

#-------------------------------------------------#
### Baseline differences in absolute ACD (2020) ###
#-------------------------------------------------#

# lm 
baseline_ACD_2020<-lm(Final_Data$ACD_2020~Final_Data$Site) 
summary(baseline_ACD_2020)
# Overall statistical differences between the sites at baseline: 
# Final_Data$SiteSBE  -20.288      5.373  -3.776 0.000221 ***

# Danum has 20.29 Mg C/ha more than SBE in 2020
# This is expected because it is likely that the logged site would have 
# less carbon than the Old Growth site 

# Check test assumptions
plot(baseline_ACD_2020) # QQ plot okay
# residual vs fitted okay
baseline_ACD_2020_residuals<-residuals(baseline_ACD_2020)
hist(baseline_ACD_2020_residuals) # Normally distributed residuals

# Assumptions of lm are met

#///////////////////#
###### Cover20 ######
#///////////////////#

#--------------------------------------------#
### Baseline differences in Cover20 (2020) ###
#--------------------------------------------#

# lm
baseline_cover2020<-lm(Final_Data$Cover_20_2020_perc~Final_Data$Site) 
summary(baseline_cover2020)
# significant baseline differences in gap volume change:
# Final_Data$SiteSBE   -1.704      1.422  -1.198    0.232 
# No difference in cover20

# Check test assumptions
plot(baseline_cover20) # QQ plot okay
# residual vs fitted okay
baseline_cover20_residuals<-residuals(baseline_cover20)
hist(baseline_cover20_residuals) # residuals okay

# Assumptions of lm are met

#/////////////////////////#
###### Canopy Height ######
#/////////////////////////#

#--------------------------------------------------#
### Baseline differences in Canopy height (2020) ###
#--------------------------------------------------#

# lm controlling for planting and logging intensity
baseline_TCH_2020<-lm(Final_Data$TCH_2020~Final_Data$Site) 
summary(baseline_TCH_2020)
# significant differences in gap volume change:
# Final_Data$SiteSBE  -6.0142     0.5270  -11.41   <2e-16 ***
# Greater canopy height in Old growth than SBE as expect due to no logging

# Check test assumptions
plot(baseline_TCH_2013) # QQ okay
# residual vs fitted okay
TCH_2013_residuals<-residuals(baseline_TCH_2013)
hist(TCH_2013_residuals) # residuals okay

#--------------------#
###### Plotting ######
#--------------------#

#///////////////////////#
### Danum vs SBE_2013 ###
#///////////////////////#

# Plot figure
cairo_pdf("Danum vs SBE2013.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(2,2),mar=c(6,8,2,2),las=1,xpd=T)

## TCH 2013 - here just compare in one year (2013)
min(Final_Data$TCH_2013) # 21.04785
max(Final_Data$TCH_2013) # 42.84278

plot(1,bty="l",xlab="",ylab="Canopy height in\n 2013 (m)", type="n",yaxt="n",ylim=c(15,45),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$TCH_2013[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$TCH_2013[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(15,45,5),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## Cover_20 in 2013
min(Final_Data$Cover_20_2013_perc) # 50.61
max(Final_Data$Cover_20_2013_perc) # 91.12746

plot(1,bty="l",xlab="",ylab="Cover20 in 2013 (%)", type="n",yaxt="n",ylim=c(40,100),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$Cover_20_2013_perc[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Cover_20_2013_perc[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(40,100,10),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## ACD 2013
max(Final_Data$ACD_2013) # 237.269
min(Final_Data$ACD_2013) # 66.10165

plot(1,bty="l",xlab="",ylab="ACD in 2013\n (Mg C/ha)", type="n",yaxt="n",ylim=c(50,300),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$ACD_2013[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_2013[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(50,300,50),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex= 1.6, col = "grey20")

dev.off()

#//////////////////////////#
### Danum vs SBE_changes ###
#//////////////////////////#

# Plot figure
cairo_pdf("Danum vs SBE_change.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(2,2),mar=c(6,8,2,2),las=1,xpd=T)

## TCH change
max(Final_Data$TCH_change_yr) # 0.4859788
min(Final_Data$TCH_change_yr) # -0.4442856

plot(1,bty="l",xlab="",ylab="Average canopy height\n change (m/yr)", type="n",yaxt="n",ylim=c(-1,1),xlim=c(0.65,2.1),cex.lab=2,axes=F, las=1)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$TCH_change_yr[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$TCH_change_yr[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(-1,1,0.5),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## New Gap volume  

max(Final_Data$New_gap_volume_yr) # 955.1908
min(Final_Data$New_gap_volume_yr) #18.59573

plot(1,bty="l",xlab="",ylab="Volume of\n New gaps (m3/ha/yr)", type="n",yaxt="n",ylim=c(0,1200),xlim=c(0.65,2.1),cex.lab=2,axes=F, las=1)

clipplot(abline(h=0,col="grey20",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(0,1200,200),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex= 1.6, col = "grey20")


## Existing Gap volume change 
max(Final_Data$Existing_proportion_closed_perc_yr) # 14.28571
min(Final_Data$Existing_proportion_closed_perc_yr) # 3.013316

plot(1,bty="l",xlab="",ylab="Proportion of Existing\n gaps closed (%/yr)", type="n",yaxt="n",ylim=c(0,18),xlim=c(0.65,2.1),cex.lab=2,axes=F, las=1)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(0,18,3),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## ACD change 

max(Final_Data$ACD_change) # 5.388675
min(Final_Data$ACD_change) # -4.778837

plot(1,bty="l",xlab="",ylab="ACD change\n (Mg C/ha/yr)", type="n",yaxt="n",ylim=c(-6,6),xlim=c(0.65,2.1),cex.lab=2,axes=F, las=1)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$ACD_change[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_change[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(-6,6,2),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

dev.off()

#///////////////////////#
### Danum vs SBE_2020 ###
#///////////////////////#

# Plot figure
cairo_pdf("Danum vs SBE_2020.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(2,2),mar=c(6,8,2,2),las=1,xpd=T)

## TCH 2020 - here just compare in one year (2020)
min(Final_Data$TCH_2020) # 22.36554
max(Final_Data$TCH_2020) # 43.20087

plot(1,bty="l",xlab="",ylab="Canopy height in\n 2020 (m)", type="n",yaxt="n",ylim=c(15,45),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$TCH_2020[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$TCH_2020[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(15,45,5),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## Cover_20 in 2020
min(Final_Data$Cover_20_2020_perc) # 42.38
max(Final_Data$Cover_20_2020_perc) # 92.18917

plot(1,bty="l",xlab="",ylab="Cover20 in 2020 (%)", type="n",yaxt="n",ylim=c(40,100),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(40,100,10),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex = 1.6, col = "grey20")

## ACD 2013
max(Final_Data$ACD_2020) # 249.8158
min(Final_Data$ACD_2020) # 88.19072

plot(1,bty="l",xlab="",ylab="ACD in 2020\n (Mg C/ha)", type="n",yaxt="n",ylim=c(50,300),xlim=c(0.65,2.1),cex.lab=2,axes=F)
clipplot(abline(h=0,col="grey20",lty=2),xlim=c(par("usr")[1],par("usr")[2]))
wvioplot(Final_Data$ACD_2020[Final_Data$Site=="SBE"],at=1,add=T,
         col=alpha("#9A3268",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_2020[Final_Data$Site=="Danum"],at=1.75,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(50,300,50),cex.axis=1.6,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at =c(1, 1.75), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("SBE logged\n forest","Primary\n forest"),
      at = c(1, 1.75), padj = 1.5, cex= 1.6, col = "grey20")

dev.off()