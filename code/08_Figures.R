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

##%%%%%%%%%%%%%%%%%##

##### LOAD DATA #####

##%%%%%%%%%%%%%%%%%##

load('Final_DataV2.rda')
Final_Data$TCH_change_yr<-(Final_Data$TCH_change)/7

SBE <- subset(Final_Data, Site == "SBE")

#%%%%%%%%%%%%%%%%%%%%%%%%%#

#### BASELINE TESTING ####

#%%%%%%%%%%%%%%%%%%%%%%%%%#

#//////////////////////////////////////#
###### Logging intensity (Cover50) #####
#//////////////////////////////////////#

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
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)

# these lines get rid of extra lines, make them the same colour as background:
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4) 
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,15,2.5),cex.axis=2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-3)
axis(1, at = seq(1, 3.5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas removed"),
      at = seq(1, 3, 1), padj = 1.5, cex = 1.8, col = "grey20")

dev.off()

#///////////////////////////#
##### Danum vs SBE_2013 #####
#///////////////////////////#

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

#//////////////////////////////#
##### Danum vs SBE_changes #####
#//////////////////////////////#

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

#///////////////////////////#
##### Danum vs SBE_2020 #####
#///////////////////////////#

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### HYPOTHESIS TESTING ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%#

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	HYPOTHESIS 1: ACD ~ Species Richness #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Update species_richness to "Primary forest" where site is "Danum" using dplyr
Spec_Data <- Final_Data %>%
  mutate(Species_richness = if_else(Site == "Danum", "Primary forest", Species_richness))

# Create a combined species richness variable
Spec_Data$combinedsr<-with(Spec_Data, interaction(Climbers_cut, Species_richness), drop = TRUE)

# Set as factors
Spec_Data$combinedsr<-factor(Spec_Data$combinedsr)

# Check the levels before removing
levels(Spec_Data$combinedsr)

# Remove the "Yes.16" level- as we are not looking at liana removal here
Spec_Data <- Spec_Data[Spec_Data$combinedsr != "Yes.16", ]
levels(Spec_Data$combinedsr)

# Define custom labels for the groups
mylabels <- c("Control", "1-species", "4-species", "16-species", "Primary forest")
Spec_Data$combinedsr <- factor(Spec_Data$combinedsr, levels = c("No.0", "No.1", "No.4","No.16", "No.Primary forest"), labels = mylabels)
levels(Spec_Data$combinedsr)

#------------#
# ACD change #
#------------#

# Differences between SBE treatments
cairo_pdf("H1Short.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(2,2),mar=c(6,6,2,2),las=1,xpd=T)

min(Spec_Data$ACD_change) # -4.778837
max(Spec_Data$ACD_change) # 5.388675

plot(1,bty="l",xlab="",ylab="ACD change\n (Mg C/ha/yr)", type="n",yaxt="n",ylim=c(-6,6),xlim=c(0.2,5.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.55,par("usr")[2]))
wvioplot(Spec_Data$ACD_change[Spec_Data$combinedsr=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_change[Spec_Data$combinedsr=="1-species"],at=2,add=T,
         col=alpha("#BAA5C3",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_change[Spec_Data$combinedsr=="4-species"],at=3,add=T,
         col=alpha("#836192",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_change[Spec_Data$combinedsr=="16-species"],at=4,add=T,
         col=alpha("#481F54",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_change[Spec_Data$combinedsr=="Primary forest"],at=5,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4)
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(-6,6,2.4),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control","1\n species","4\n species",
                         "16\n species", "Primary\n forest"),
      at = seq(1, 5, 1), padj = 1.5, cex = 1, col = "grey20")

#---------------------------#
# Change in Gap Volume- NEW #
#---------------------------#

min(Spec_Data$New_gap_volume_yr) # 18.59573
max(Spec_Data$New_gap_volume_yr) # 955.1908

plot(1,bty="l",xlab="",ylab="Volume of\n new gaps (m3/ha/yr)", type="n",yaxt="n",ylim=c(0,1200),xlim=c(0.2,5.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Spec_Data$New_gap_volume_yr[Spec_Data$combinedsr=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$New_gap_volume_yr[Spec_Data$combinedsr=="1-species"],at=2,add=T,
         col=alpha("#BAA5C3",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$New_gap_volume_yr[Spec_Data$combinedsr=="4-species"],at=3,add=T,
         col=alpha("#836192",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$New_gap_volume_yr[Spec_Data$combinedsr=="16-species"],at=4,add=T,
         col=alpha("#481F54",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$New_gap_volume_yr[Spec_Data$combinedsr=="Primary forest"],at=5,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4)
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,1200,200),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control","1\n species","4\n species",
                         "16\n species", "Primary\n forest"),
      at = seq(1, 5, 1), padj = 1.5, cex = 1, col = "grey20")

#--------------------------------#
# Change in Gap Volume- Existing #
#--------------------------------#

min(Spec_Data$Existing_proportion_closed_perc_yr) # 3.013316
max(Spec_Data$Existing_proportion_closed_perc_yr) # 14.28571

plot(1,bty="l",xlab="",ylab="Proportion of existing\n gaps closed (%/yr)", type="n",yaxt="n",ylim=c(0,18),xlim=c(0.2,5.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Spec_Data$Existing_proportion_closed_perc_yr[Spec_Data$combinedsr=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Existing_proportion_closed_perc_yr[Spec_Data$combinedsr=="1-species"],at=2,add=T,
         col=alpha("#BAA5C3",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Existing_proportion_closed_perc_yr[Spec_Data$combinedsr=="4-species"],at=3,add=T,
         col=alpha("#836192",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Existing_proportion_closed_perc_yr[Spec_Data$combinedsr=="16-species"],at=4,add=T,
         col=alpha("#481F54",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Existing_proportion_closed_perc_yr[Spec_Data$combinedsr=="Primary forest"],at=5,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4)
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,18,3),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-2)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control","1\n species","4\n species",
                         "16\n species", "Primary\n forest"),
      at = seq(1, 5, 1), padj = 1.5, cex = 1, col = "grey20")

dev.off()

#////////////////#
### Long-term ###
#////////////////#

# Differences between SBE treatments
cairo_pdf("H1Long.pdf",
          width = 12,
          height = 6,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,2),mar=c(6,6,2,2),las=1,xpd=T)

#------------#
# ACD change #
#------------#

min(Spec_Data$ACD_2020) # 88.19072
max(Spec_Data$ACD_2020) # 249.8158

plot(1,bty="l",xlab="",ylab="ACD 2020 (Mg C/ha)", type="n",yaxt="n",ylim=c(80,260),xlim=c(0.2,5.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(par("usr")[1],4.5))
wvioplot(Spec_Data$ACD_2020[Spec_Data$combinedsr=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_2020[Spec_Data$combinedsr=="1-species"],at=2,add=T,
         col=alpha("#BAA5C3",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_2020[Spec_Data$combinedsr=="4-species"],at=3,add=T,
         col=alpha("#836192",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_2020[Spec_Data$combinedsr=="16-species"],at=4,add=T,
         col=alpha("#481F54",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$ACD_2020[Spec_Data$combinedsr=="Primary forest"],at=5,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(80,260,30),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control","1\n species","4\n species",
                         "16\n species", "Primary\n forest"),
      at = seq(1, 5, 1), padj = 1.5, cex = 1, col = "grey20")

#---------#
# Cover20 #
#---------#

min(Spec_Data$Cover_20_2020_perc) # 42.38
max(Spec_Data$Cover_20_2020_perc) # 92.18917

plot(1,bty="l",xlab="",ylab="Cover20 in 2020 (%)", type="n",yaxt="n",ylim=c(40,100),xlim=c(0.2,5.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(par("usr")[1],4.5))
wvioplot(Spec_Data$Cover_20_2020_perc[Spec_Data$combinedsr=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Cover_20_2020_perc[Spec_Data$combinedsr=="1-species"],at=2,add=T,
         col=alpha("#BAA5C3",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Cover_20_2020_perc[Spec_Data$combinedsr=="4-species"],at=3,add=T,
         col=alpha("#836192",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Cover_20_2020_perc[Spec_Data$combinedsr=="16-species"],at=4,add=T,
         col=alpha("#481F54",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Spec_Data$Cover_20_2020_perc[Spec_Data$combinedsr=="Primary forest"],at=5,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(40,100,10),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control","1\n species","4\n species",
                         "16\n species", "Primary\n forest"),
      at = seq(1, 5, 1), padj = 1.5, cex = 1, col = "grey20")

dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	HYPOTHESIS 2: ACD ~ Planting/Lianas #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#------------------------------------#
# Differences between SBE treatments #
#------------------------------------#

cairo_pdf("H2Short.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(2,2),mar=c(8,6,2,2),las=1,xpd=T)

## ACD change
min(Final_Data$ACD_change) # -4.778837
max(Final_Data$ACD_change) # 5.388675

plot(1,bty="l",xlab="",ylab="ACD change (Mg C/ha/yr)", type="n",yaxt="n",ylim=c(-6,6),xlim=c(0.5,4.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$ACD_change[Final_Data$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_change[Final_Data$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_change[Final_Data$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_change[Final_Data$combined=="Primary forest"],at=4,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
# these lines get rid of extra lines, make them the same colour as background:
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6) 
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(-6,6,2),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 4, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas\n removed", "Primary\n forest"),
      at = seq(1, 4, 1), padj = 1.5, cex = 1, col = "grey20", srt=45)

#---------------------------#
# Change in Gap Volume- New #
#---------------------------#

min(Final_Data$New_gap_volume_yr) # 18.59573
max(Final_Data$New_gap_volume_yr) # 955.1908

plot(1,bty="l",xlab="",ylab="Volume of new\n gaps (m3/ha/yr)", type="n",yaxt="n",ylim=c(0,1200),xlim=c(0.5,4.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$New_gap_volume_yr[Final_Data$combined=="Primary forest"],at=4,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(0,1200,200),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1.5)
axis(1, at = seq(1, 4, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas\n removed", "Primary\n forest"),
      at = seq(1, 4, 1), padj = 1.5, cex = 1, col = "grey20")

#--------------------------------#
# Change in Gap Volume- Existing #
#--------------------------------#

min(Final_Data$Existing_proportion_closed_perc_yr) # 3.013316
max(Final_Data$Existing_proportion_closed_perc_yr) # 14.28571

plot(1,bty="l",xlab="",ylab="Proportion of existing\n gaps closed (%/yr)", type="n",yaxt="n",ylim=c(0,18),xlim=c(0.5,4.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.65,par("usr")[2]))
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Existing_proportion_closed_perc_yr[Final_Data$combined=="Primary forest"],at=4,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6)
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(0,18,3),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1.5)
axis(1, at = seq(1, 4, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas\n removed", "Primary\n forest"),
      at = seq(1, 4, 1), padj = 1.5, cex = 1, col = "grey20")

dev.off()

#///////////////#
### Long-term ###
#///////////////#

#------------------------------------#
# Differences between SBE treatments #
#------------------------------------#

cairo_pdf("H2Long.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,2),mar=c(8,6,2,2),las=1,xpd=T)

## ACD change
min(Final_Data$ACD_2020) # 88.19072
max(Final_Data$ACD_2020) # 249.8158

plot(1,bty="l",xlab="",ylab="ACD in 2020 (Mg C/ha)", type="n",yaxt="n",ylim=c(0,300),xlim=c(0.5,4.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(par("usr")[1],4.5))
wvioplot(Final_Data$ACD_2020[Final_Data$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_2020[Final_Data$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_2020[Final_Data$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$ACD_2020[Final_Data$combined=="Primary forest"],at=4,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6) 
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(0,300,50),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 4, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas\n removed", "Primary\n forest"),
      at = seq(1, 4, 1), padj = 1.5, cex = 1, col = "grey20")

#----------------#
# Cover20 change #
#----------------#

min(Final_Data$Cover_20_2020_perc) # 42.38
max(Final_Data$Cover_20_2020_perc) # 92.18917

plot(1,bty="l",xlab="",ylab="Cover20 in 2020 (%)", type="n",yaxt="n",ylim=c(40,100),xlim=c(0.5,4.5),cex.lab=1.6,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(par("usr")[1],4.5))
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$combined=="Control"],at=1,add=T,
         col=alpha("#8EB2C4",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$combined=="Planted"],at=2,add=T,
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$combined=="Planted with lianas removed"],at=3,add=T,
         col=alpha("#DC8B7F",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(Final_Data$Cover_20_2020_perc[Final_Data$combined=="Primary forest"],at=4,add=T,
         col=alpha("#D59AAE",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=6);abline(h=par("usr")[3],col="white",lwd=6) 
abline(v=par("usr")[1],col="white",lwd=6);abline(v=par("usr")[2],col="white",lwd=6)
axis(2,at=seq(40,100,10),cex.axis=1.2,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 4, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Unplanted\n control", "\nPlanted", "Planted with\n lianas\n removed", "Primary\n forest"),
      at = seq(1, 4, 1), padj = 1.5, cex = 1, col = "grey20")

dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##### HYPOTHESIS 3: Gap Volume~ Planting + HALP ##### 

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#----------------------------#
### Predicting from model  ###
#----------------------------#

# Predict changes
new.data <- data.frame(
  combined = rep(factor(c("Control", "Planted", "Planted with lianas removed")), each = 100),
  HALP = rep(seq(min(SBE$HALP, na.rm = TRUE), max(SBE$HALP, na.rm = TRUE), len = 100), times = 3)
)

levels(new.data$combined)

# Models
ACD_HALP<-lm(ACD_change~combined+HALP, data=SBE)
summary(ACD_HALP) # HALP important for ACD in short term

GV_HALP_short_NEW<-lm(New_gap_volume_yr~combined+HALP, data=SBE)
summary(GV_HALP_short_NEW) # HALP not important for GV in short term

GV_HALP_short_Existing<-lm(Existing_proportion_closed_perc_yr~combined+HALP, data=SBE)
summary(GV_HALP_short_Existing) # HALP  important for GV in short term

# Log transform the data for new gap volumes
SBE$log_New_gap_volume_yr <- log(SBE$New_gap_volume_yr)  

# Fit the model using log-transformed data
GV_HALP_short_NEW_log <- lm(log_New_gap_volume_yr ~ combined + HALP, data = SBE)

# Predict
new.data$ACD_Change<- predict(ACD_HALP, new.data)
new.data$New_Gap_volume_change_hect<- predict(GV_HALP_short_NEW, new.data)
new.data$Existing_proportion_closed_perc_yr<- predict(GV_HALP_short_Existing, new.data)

# Predict on the log-transformed data
new.data$log_New_Gap_volume_change_hect <- predict(GV_HALP_short_NEW_log, new.data)

# Back transform the predicted values to the original scale
new.data$New_Gap_volume_change_hect <- exp(new.data$log_New_Gap_volume_change_hect)

# Get mean and max values for axis lims
min(SBE$Existing_proportion_closed_perc_yr) #4
max(SBE$Existing_proportion_closed_perc_yr) #12

min(SBE$ACD_change) #-4.778837
max(SBE$ACD_change) #5.388675

min(SBE$New_gap_volume_yr) # 40.07824
max(SBE$New_gap_volume_yr) # 835.5748 

# Function to calculate confidence interval
calculate_ci <- function(model, data, level = 0.95) {
  preds <- predict(model, newdata = data, interval = "confidence", level = level)
  data.frame(HALP = data$HALP, ymin = preds[, "lwr"], ymax = preds[, "upr"])
}

# Confidence interval level
ci_level <- 0.95

#////////////////#
### Short-term ###
#////////////////#

cairo_pdf("H3.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")

par(family = "Cambria")
par(mfrow = c(2,2), mar = c(6,8,2,2), las = 1, xpd = TRUE)

# Treatment colours
colors <- c(
  "Control" = "#8EB2C4",
  "Planted" = "#652D76",
  "Planted with lianas removed" = "#DC8B7F"
)

# Safe colour mapping for observed data
point_groups <- as.character(SBE$combined)
point_cols <- colors[point_groups]
point_cols[is.na(point_cols)] <- "grey50"


#================#
### ACD CHANGE ###
#================#

plot(ACD_change ~ HALP,
     SBE,
     pch = 16,
     col = point_cols,
     cex = 1.2,
     xlim = c(0,250),
     ylim = c(-6,6),
     xlab = "",
     ylab = "",
     axes = FALSE)

mtext("Height above nearest drainage (m)",1,line=4,cex=1.2,col="grey20")
mtext("ACD change (Mg C/ha/yr)",2,line=4,cex=1.2,las=3,col="grey20")

axis(1,at=seq(0,250,50),cex.axis=1.2,col="grey20",col.axis="grey20")
axis(2,at=seq(-6,6,2),cex.axis=1.2,col="grey20",col.axis="grey20")

clipplot(abline(h=0,col="grey20",lty=2),
         xlim=c(par("usr")[1],par("usr")[2]))

for (grp in levels(new.data$combined)) {
  
  data_subset <- new.data[new.data$combined == grp,]
  model_preds <- calculate_ci(ACD_HALP, data_subset, level = ci_level)
  
  line_color <- colors[grp]
  
  polygon(
    c(model_preds$HALP, rev(model_preds$HALP)),
    c(model_preds$ymin, rev(model_preds$ymax)),
    col = alpha(line_color,0.3),
    border = NA
  )
  
  lines(ACD_Change ~ HALP,
        data_subset,
        lwd = 3,
        col = line_color)
}


#===========================#
### NEW GAP FORMATION ###
#===========================#

plot(New_gap_volume_yr ~ HALP,
     SBE,
     pch = 16,
     col = point_cols,
     cex = 1.2,
     xlim = c(0,250),
     ylim = c(0,900),
     xlab = "",
     ylab = "",
     axes = FALSE)

mtext("Height above nearest drainage (m)",1,line=4,cex=1.2,col="grey20")
mtext("Volume of new gaps (m3/yr)",2,line=4,cex=1.2,las=3,col="grey20")

axis(1,at=seq(0,250,50),cex.axis=1.2,col="grey20",col.axis="grey20")
axis(2,at=seq(0,900,150),cex.axis=1.2,col="grey20",col.axis="grey20")

for (grp in levels(new.data$combined)) {
  
  data_subset <- new.data[new.data$combined == grp,]
  model_preds <- calculate_ci(GV_HALP_short_NEW_log, data_subset, level = ci_level)
  
  model_preds$ymin <- exp(model_preds$ymin)
  model_preds$ymax <- exp(model_preds$ymax)
  
  line_color <- colors[grp]
  
  polygon(
    c(model_preds$HALP, rev(model_preds$HALP)),
    c(model_preds$ymin, rev(model_preds$ymax)),
    col = alpha(line_color,0.3),
    border = NA
  )
  
  lines(New_Gap_volume_change_hect ~ HALP,
        data_subset,
        lwd = 3,
        col = line_color)
}


#==============================#
### EXISTING GAP CLOSURE ###
#==============================#

plot(Existing_proportion_closed_perc_yr ~ HALP,
     SBE,
     pch = 16,
     col = point_cols,
     cex = 1.2,
     xlim = c(0,250),
     ylim = c(0,18),
     xlab = "",
     ylab = "",
     axes = FALSE)

mtext("Height above nearest drainage (m)",1,line=4,cex=1.2,col="grey20")
mtext("Proportion of existing\n gaps closed (%/yr)",2,line=4.5,cex=1.2,las=3,col="grey20")

axis(1,at=seq(0,250,50),cex.axis=1.2,col="grey20",col.axis="grey20")
axis(2,at=seq(0,18,3),cex.axis=1.2,col="grey20",col.axis="grey20")

clipplot(abline(h=0,col="grey20",lty=2),
         xlim=c(par("usr")[1],par("usr")[2]))

for (grp in levels(new.data$combined)) {
  
  data_subset <- new.data[new.data$combined == grp,]
  model_preds <- calculate_ci(GV_HALP_short_Existing, data_subset, level = ci_level)
  
  line_color <- colors[grp]
  
  polygon(
    c(model_preds$HALP, rev(model_preds$HALP)),
    c(model_preds$ymin, rev(model_preds$ymax)),
    col = alpha(line_color,0.3),
    border = NA
  )
  
  lines(Existing_proportion_closed_perc_yr ~ HALP,
        data_subset,
        lwd = 3,
        col = line_color)
}

dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%##

#### SUPPLEMENTARY STATS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%##

# Check predictors for correlation
corr_halp<-cor(SBE$HALP, SBE$Frac_Over_50m_perc ,use = "complete.obs")
print(corr_halp)
# HALP- Frac over 50 not too highly correlated: 0.25431668

frac_halp <- cor(SBE[, c("ACD_2020", "ACD_change", "Cover_20_2020_perc", "Frac_Over_50m_perc", "HALP", "Existing_proportion_closed_perc_yr", "New_gap_volume_yr")],use = "complete.obs")
print(frac_halp)

#-----------------------#
###### Correlation ######
#-----------------------#

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


