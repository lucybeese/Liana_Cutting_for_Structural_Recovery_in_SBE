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

##%%%%%%%%%%%%%%%##

#### LOAD DATA ####

##%%%%%%%%%%%%%%%##

load('Final_Data.rda')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	HYPOTHESIS 1: ACD ~ Species Richness #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Response ~ species richness (lm where species richness is a factor with 4 levels: controls, 1-sp mix, 4-sp mix, 16-sp mix)

#----------------------#
###### Short term ######
#----------------------#

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

ACD_H3<-lm(Spec_Data$ACD_change~Spec_Data$combinedsr)
summary(ACD_H3) 
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)                      -0.1187     0.4656  -0.255    0.799
# Spec_Data$combinedsr1-species    0.1489     0.5459   0.273    0.785
# Spec_Data$combinedsr4-species    0.4844     0.5459   0.887    0.376
# Spec_Data$combinedsr16-species   0.1855     0.5459   0.340    0.735
# Spec_Data$combinedsrOld Growth   0.6234     0.5240   1.190    0.236

# Check assumptions of lm
plot(ACD_H3) # QQ okay
# Check for homoscedasticity using residuals vs fitted values plot, looks okay
ACD_H3_residuals<-residuals(ACD_H3)
hist(ACD_H3_residuals) # Normally distributed
# Assumptions met

# Obtain estimated marginal means
emmACD_H3 <- emmeans(ACD_H3, ~ combinedsr)

# Perform pairwise comparisons
contrast_emmACD_H3 <- pairs(emmACD_H3)

# Print the results
print(contrast_emmACD_H3)

# contrast                   estimate    SE  df t.ratio p.value
# Control - (1-species)       -0.1489 0.546 148  -0.273  0.9988
# Control - (4-species)       -0.4844 0.546 148  -0.887  0.9012
# Control - (16-species)      -0.1855 0.546 148  -0.340  0.9971
# Control - Old Growth        -0.6234 0.524 148  -1.190  0.7572
# (1-species) - (4-species)   -0.3356 0.403 148  -0.832  0.9202
# (1-species) - (16-species)  -0.0366 0.403 148  -0.091  1.0000
# (1-species) - Old Growth    -0.4746 0.373 148  -1.272  0.7085
# (4-species) - (16-species)   0.2990 0.403 148   0.741  0.9464
# (4-species) - Old Growth    -0.1390 0.373 148  -0.373  0.9959
# (16-species) - Old Growth   -0.4380 0.373 148  -1.174  0.7660

# No significant PAIRWISE differences are observed
# In the short term, there seems to be no  effect of Species richness on ACD. 

#---------------------#
###### Long term ######
#---------------------#

# ANOVA
ACD_H3_long<-lm(Spec_Data$ACD_2020~Spec_Data$combinedsr) 
summary(ACD_H3_long)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         135.962      8.929  15.227  < 2e-16 ***
#   Spec_Data$combinedsr1-species         7.771     10.470   0.742  0.45913    
# Spec_Data$combinedsr4-species         7.095     10.470   0.678  0.49906    
# Spec_Data$combinedsr16-species       13.774     10.470   1.316  0.19037    
# Spec_Data$combinedsrPrimary forest   31.494     10.049   3.134  0.00208 ** 
#   ---

plot(ACD_H3_long) # residuals off
ACD_H3_long_residuals<-residuals(ACD_H3_long)
hist(ACD_H3_long_residuals) # Okay 

# Obtain estimated marginal means
emmACD_H3_long <- emmeans(ACD_H3_long, ~ combinedsr)

# Perform pairwise comparisons
contrast_emmACD_H3_long <- pairs(emmACD_H3_long)

# Print the results
print(contrast_emmACD_H3_long)

# contrast                   estimate    SE  df t.ratio p.value
# Control - (1-species)        -7.771 10.47 148  -0.742  0.9462
# Control - (4-species)        -7.095 10.47 148  -0.678  0.9610
# Control - (16-species)      -13.774 10.47 148  -1.316  0.6820
# Control - Old Growth        -31.494 10.05 148  -3.134  0.0175
# (1-species) - (4-species)     0.676  7.73 148   0.087  1.0000
# (1-species) - (16-species)   -6.003  7.73 148  -0.776  0.9371
# (1-species) - Old Growth    -23.723  7.15 148  -3.317  0.0099
# (4-species) - (16-species)   -6.679  7.73 148  -0.864  0.9097
# (4-species) - Old Growth    -24.399  7.15 148  -3.411  0.0073
# (16-species) - Old Growth   -17.721  7.15 148  -2.478  0.1013

# Significant differences between:
#(1-species) - Old Growth    -23.723  7.15 148  -3.317  0.0099
#(4-species) - Old Growth    -24.399  7.15 148  -3.411  0.0073
# Control - Old Growth        -31.494 10.05 148  -3.134  0.0175

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	HYPOTHESIS 1: Gap Volume ~ Species Richness #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#----------------------#
###### Short term ######
#----------------------#

#/////////#
### NEW ###
#/////////#

# ANOVA
GV_H3_New_lm<-lm(Spec_Data$New_gap_volume_yr~Spec_Data$combinedsr) 
summary(GV_H3_New_lm)

# Spec_Data$combinedsr1-species    -22.38      58.34  -0.384   0.7019    
# Spec_Data$combinedsr4-species    -37.81      58.34  -0.648   0.5180    
# Spec_Data$combinedsr16-species   -37.47      58.34  -0.642   0.5217    
# Spec_Data$combinedsrOld Growth  -141.20      56.00  -2.522   0.0127 *  

plot(GV_H3_New_lm) # QQ okay
# residual vs fitted slightly off
GV_H3_New_residuals<-residuals(GV_H3_New_lm)
hist(GV_H3_New_residuals) # Normal

# Assumptions met

# Obtain estimated marginal means
emmGV_H3_New <- emmeans(GV_H3_New_lm, ~ combinedsr)

# Perform pairwise comparisons
contrast_emmGV_H3_New <- pairs(emmGV_H3_New)

# Print the results
print(contrast_emmGV_H3_New)

# contrast                   estimate   SE  df t.ratio p.value
# Control - (1-species)        22.375 58.3 148   0.384  0.9954
# Control - (4-species)        37.809 58.3 148   0.648  0.9668
# Control - (16-species)       37.470 58.3 148   0.642  0.9678
# Control - Old Growth        141.201 56.0 148   2.522  0.0914
# (1-species) - (4-species)    15.434 43.1 148   0.358  0.9964
# (1-species) - (16-species)   15.095 43.1 148   0.350  0.9967
# (1-species) - Old Growth    118.826 39.9 148   2.981  0.0273
# (4-species) - (16-species)   -0.339 43.1 148  -0.008  1.0000
# (4-species) - Old Growth    103.392 39.9 148   2.594  0.0766
# (16-species) - Old Growth   103.730 39.9 148   2.603  0.0751

# pairwise sig diffs in New gap volume with species richness:
# (1-species) - Old Growth    118.826 39.9 148   2.981  0.0273

#//////////////#
### Existing ###
#//////////////#

# ANOVA
GV_H3_Existing_lm<-lm(Spec_Data$Existing_proportion_closed_perc_yr~Spec_Data$combinedsr)
summary(GV_H3_Existing_lm)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       6.5381     0.5784  11.304  < 2e-16 ***
#   Spec_Data$combinedsr1-species    0.6780     0.6782   1.000  0.31915    
# Spec_Data$combinedsr4-species    0.7993     0.6782   1.178  0.24051    
# Spec_Data$combinedsr16-species   0.4448     0.6782   0.656  0.51293    
# Spec_Data$combinedsrOld Growth   2.1811     0.6510   3.350  0.00102 **   

plot(GV_H3_Existing_lm) # QQ ok
GV_H3_Existing_residuals<-residuals(GV_H3_Existing_lm)
hist(GV_H3_Existing_residuals) # okay

# Assumptions met

# Obtain estimated marginal means
emmGV_H3_Existing <- emmeans(GV_H3_Existing_lm, ~ combinedsr)

# Perform pairwise comparisons
contrast_emmGV_H3_Existing <- pairs(emmGV_H3_Existing)

# Print the results
print(contrast_emmGV_H3_Existing)

# contrast                   estimate    SE  df t.ratio p.value
# Control - (1-species)        -0.678 0.678 148  -1.000  0.8552
# Control - (4-species)        -0.799 0.678 148  -1.178  0.7637
# Control - (16-species)       -0.445 0.678 148  -0.656  0.9653
# Control - Old Growth         -2.181 0.651 148  -3.350  0.0089
# (1-species) - (4-species)    -0.121 0.501 148  -0.242  0.9992
# (1-species) - (16-species)    0.233 0.501 148   0.465  0.9903
# (1-species) - Old Growth     -1.503 0.463 148  -3.244  0.0125
# (4-species) - (16-species)    0.354 0.501 148   0.708  0.9545
# (4-species) - Old Growth     -1.382 0.463 148  -2.982  0.0273
# (16-species) - Old Growth    -1.736 0.463 148  -3.747  0.0023

# Control - Old Growth         -2.181 0.651 148  -3.350  0.0089
# (1-species) - Old Growth     -1.503 0.463 148  -3.244  0.0125
# (4-species) - Old Growth     -1.382 0.463 148  -2.982  0.0273
# (16-species) - Old Growth    -1.736 0.463 148  -3.747  0.0023

#---------------------#
###### Long term ######
#---------------------#

# lm
C20_H3<-lm(Spec_Data$Cover_20_2020_perc~Spec_Data$combinedsr) 
summary(C20_H3)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       70.674      2.338  30.232   <2e-16 ***
#   Spec_Data$combinedsr1-species     2.216      2.741   0.809    0.420    
# Spec_Data$combinedsr4-species     2.339      2.741   0.853    0.395    
# Spec_Data$combinedsr16-species    4.467      2.741   1.629    0.105    
# Spec_Data$combinedsrOld Growth    5.200      2.631   1.976    0.050 *  

# Check assumptions
plot(C20_H3) # QQ okay
C20_H3_residuals<-residuals(C20_H3)
hist(C20_H3_residuals) # right skewed residuals
# Assumptions not met

# Obtain estimated marginal means
emmC20_H3 <- emmeans(C20_H3, ~ combinedsr)

# Perform pairwise comparisons
contrast_C20_H3 <- pairs(emmC20_H3)

# Print the results
print(contrast_C20_H3)

# contrast                   estimate   SE  df t.ratio p.value
# Control - (1-species)        -2.216 2.74 148  -0.809  0.9277
# Control - (4-species)        -2.339 2.74 148  -0.853  0.9132
# Control - (16-species)       -4.467 2.74 148  -1.629  0.4811
# Control - Old Growth         -5.200 2.63 148  -1.976  0.2827
# (1-species) - (4-species)    -0.123 2.02 148  -0.061  1.0000
# (1-species) - (16-species)   -2.250 2.02 148  -1.111  0.8003
# (1-species) - Old Growth     -2.983 1.87 148  -1.593  0.5041
# (4-species) - (16-species)   -2.127 2.02 148  -1.051  0.8311
# (4-species) - Old Growth     -2.861 1.87 148  -1.528  0.5463
# (16-species) - Old Growth    -0.733 1.87 148  -0.392  0.9950

# Non-parametric alternative (not reported in manuscript)
KWtest3 <- kruskal.test(Spec_Data$Cover_20_2020_perc~Spec_Data$combinedsr)
KWtest3  # p-value = 0.06749

# Perform pairwise Wilcoxon rank-sum tests with p-value adjustments
pairwise.wilcox.test(Spec_Data$Cover_20_2020_perc, Spec_Data$combinedsr, p.adjust.method = "bonferroni")

# Control 1-species 4-species 16-species
# 1-species  1.00    -         -         -         
#   4-species  1.00    1.00      -         -         
#   16-species 0.55    1.00      1.00      -         
#   Old Growth 0.20    0.76      0.42      1.00  

#--------------------#
###### Plotting ######
#--------------------#

#////////////////#
### Short-term ###
#////////////////#

plot.new()

# set wd to figures folder

## Load violin plot function
source(file="wvioplot.r")

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
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
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
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
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
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
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
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
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
         col=alpha("#652D76",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
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

#-----------------------#
###### Short term #######
#-----------------------#

# Create combined variable
Final_Data$combined<-with(Final_Data, interaction(Climbers_cut, Planted), drop = TRUE)

# Set as factor
Final_Data$combined<-factor(Final_Data$combined)
levels(Final_Data$combined)

# Define custom labels for the groups
mylabels <- c("Control","Primary forest", "Planted", "Planted with lianas removed")
Final_Data$combined <- factor(Final_Data$combined, levels = c("No.Control", "No.Old Growth", "No.Planted", "Yes.Planted"), labels = mylabels)

# Perform linear model 
ACD_h12<-lm(Final_Data$ACD_change~Final_Data$combined) 
summary(ACD_h12)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     -0.1187     0.4615  -0.257    0.797    
# Final_Data$combinedOld Growth                    0.6234     0.5195   1.200    0.232    
# Final_Data$combinedPlanted                       0.2729     0.4895   0.558    0.578    
# Final_Data$combinedPlanted with lianas removed   3.0201     0.6106   4.946 1.85e-06 ***

# Check assumptions of lm
plot(ACD_h12) # QQ okay
ACD_h12_residuals<-residuals(ACD_h12) # plotted vs residuals okay
hist(ACD_h12_residuals) # Normal residuals

# Assumptions met

# Obtain estimated marginal means
emm_ACD_h12 <- emmeans(ACD_h12, ~ combined)

# Perform pairwise comparisons
contrast_ACD_h12 <- pairs(emm_ACD_h12)

# Print the results
print(contrast_ACD_h12)
# 
# contrast                                 estimate    SE  df t.ratio p.value
# Control - Old Growth                       -0.623 0.519 165  -1.200  0.6276
# Control - Planted                          -0.273 0.490 165  -0.558  0.9444
# Control - Planted with lianas removed      -3.020 0.611 165  -4.946  <.0001
# Old Growth - Planted                        0.351 0.289 165   1.213  0.6192
# Old Growth - Planted with lianas removed   -2.397 0.465 165  -5.150  <.0001
# Planted - Planted with lianas removed      -2.747 0.432 165  -6.363  <.0001

# Significant differences are observed between:
# Control - Planted with lianas removed <.0001
# Planted - Planted with lianas removed <.0001
# Planted with lianas removed - Old Growth <.0001

#----------------------#
###### Long term #######
#----------------------#

# lm
ACD_h12_long<-lm(Final_Data$ACD_2020~Final_Data$combined) 
summary(ACD_h12_long) 

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     135.962      8.762  15.517  < 2e-16 ***
#   Final_Data$combinedOld Growth                    31.494      9.861   3.194  0.00168 ** 
#   Final_Data$combinedPlanted                        9.547      9.293   1.027  0.30581    
# Final_Data$combinedPlanted with lianas removed   29.571     11.591   2.551  0.01164 *  

# Check assumptions
plot(ACD_h12_long) # QQ okay
ACD_h12_long_residuals<-residuals(ACD_h12_long)
hist(ACD_h12_long_residuals) # Normal residuals
# Check for homoscedasticity using residuals vs fitted values plot - looks okay

# Obtain estimated marginal means
emm_h12_long <- emmeans(ACD_h12_long, ~ combined)

# Perform pairwise comparisons
contrast_h12_long <- pairs(emm_h12_long)

# Print the results
print(contrast_h12_long)

# contrast                                 estimate    SE  df t.ratio p.value
# Control - Old Growth                       -31.49  9.86 165  -3.194  0.0090
# Control - Planted                           -9.55  9.29 165  -1.027  0.7337
# Control - Planted with lianas removed      -29.57 11.59 165  -2.551  0.0560
# Old Growth - Planted                        21.95  5.48 165   4.003  0.0005
# Old Growth - Planted with lianas removed     1.92  8.83 165   0.218  0.9963
# Planted - Planted with lianas removed      -20.02  8.20 165  -2.443  0.0731

# Significant differences are observed between:
# Control - Old Growth                       -31.49  9.86 165  -3.194  0.0090
# Old Growth - Planted                        21.95  5.48 165   4.003  0.0005          

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#####	HYPOTHESIS 2: Gap Volume ~ Planting/Lianas #####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Testing the 2 different gap types: New, and Existing

#----------------------#
###### Short term ######
#----------------------#

length(Final_Data$New_gap_volume_yr)
length(Final_Data$Existing_proportion_closed_perc_yr)

#/////////#
### NEW ###
#/////////#

# lm
GV_h12_AOV_New<-lm(Final_Data$New_gap_volume_yr~Final_Data$combined) # significant, <2e-16
summary(GV_h12_AOV_New)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      460.53      48.77   9.443  < 2e-16 ***
#   GV_H12_New$combinedOld Growth                   -141.20      54.89  -2.573    0.011 *  
#   GV_H12_New$combinedPlanted                       -32.55      51.73  -0.629    0.530    
# GV_H12_New$combinedPlanted with lianas removed  -262.89      64.52  -4.075 7.14e-05 ***

# Check Assumptions
plot(GV_h12_AOV_New) # QQ okay
GV_h12_AOV_New_residuals<-residuals(GV_h12_AOV_New)
hist(GV_h12_AOV_New_residuals) # Normal
# Check for homoscedasticity using residuals vs fitted values plot - looks okay

# Obtain estimated marginal means
emm_GV_h12_AOV_New <- emmeans(GV_h12_AOV_New, ~ combined)

# Perform pairwise comparisons
contrast_GV_h12_AOV_New <- pairs(emm_GV_h12_AOV_New)

# Print the results
print(contrast_GV_h12_AOV_New)

# contrast                                 estimate   SE  df t.ratio p.value
# Control - Old Growth                        141.2 54.9 165   2.573  0.0530
# Control - Planted                            32.6 51.7 165   0.629  0.9225
# Control - Planted with lianas removed       262.9 64.5 165   4.075  0.0004
# Old Growth - Planted                       -108.6 30.5 165  -3.560  0.0027
# Old Growth - Planted with lianas removed    121.7 49.2 165   2.475  0.0677
# Planted - Planted with lianas removed       230.3 45.6 165   5.049  <.0001

# Significant differences are observed between:
# Control - Planted with lianas removed       262.9 64.5 165   4.075  0.0004
# Old Growth - Planted                       -108.6 30.5 165  -3.560  0.0027
# Planted - Planted with lianas removed       230.3 45.6 165   5.049  <.0001

#//////////////#
### Existing ###
#//////////////#

# ANOVA
GV_h12_AOV_Existing<- lm(Final_Data$Existing_proportion_closed_perc_yr~Final_Data$combined)
summary(GV_h12_AOV_Existing)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      6.5381     0.5700  11.469  < 2e-16 ***
#   Final_Data$combinedOld Growth                    2.1811     0.6416   3.400 0.000846 ***
#   Final_Data$combinedPlanted                       0.6407     0.6046   1.060 0.290854    
# Final_Data$combinedPlanted with lianas removed   2.8898     0.7541   3.832 0.000180 ***

plot(GV_h12_AOV_Existing) # QQ off
GV_h12_AOV_Existing_residuals<-residuals(GV_h12_AOV_Existing)
hist(GV_h12_AOV_Existing_residuals) # normal
# Assumptions met

# Obtain estimated marginal means
emm_GV_h12_AOV_Existing <- emmeans(GV_h12_AOV_Existing, ~ combined)

# Perform pairwise comparisons
contrast_GV_h12_AOV_Existing <- pairs(emm_GV_h12_AOV_Existing)

# Print the results
print(contrast_GV_h12_AOV_Existing)

# contrast                                 estimate    SE  df t.ratio p.value
# Control - Old Growth                       -2.181 0.642 165  -3.400  0.0047
# Control - Planted                          -0.641 0.605 165  -1.060  0.7145
# Control - Planted with lianas removed      -2.890 0.754 165  -3.832  0.0010
# Old Growth - Planted                        1.540 0.357 165   4.318  0.0002
# Old Growth - Planted with lianas removed   -0.709 0.575 165  -1.233  0.6068
# Planted - Planted with lianas removed      -2.249 0.533 165  -4.218  0.0002

# significant PAIRWISE differences are observed:
# Control - Old Growth                       -2.181 0.642 165  -3.400  0.0047
# Old Growth - Planted                        1.540 0.357 165   4.318  0.0002
# Control - Planted with lianas removed      -2.890 0.754 165  -3.832  0.0010
# Planted - Planted with lianas removed      -2.249 0.533 165  -4.218  0.0002

#----------------------#
###### Long term #######
#----------------------#

# ANOVA
C20_h12_AOV_long<-lm(Final_Data$Cover_20_2020_perc~Final_Data$combined) 
summary(C20_h12_AOV_long) 

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      70.674      2.301  30.718  < 2e-16 ***
#   Final_Data$combinedOld Growth                     5.200      2.589   2.008  0.04626 *  
#   Final_Data$combinedPlanted                        3.008      2.440   1.232  0.21954    
# Final_Data$combinedPlanted with lianas removed    9.050      3.044   2.974  0.00338 ** 

# Check assumptions
plot(C20_h12_AOV_long) # QQ a little off at ends but okay
C20_h12_AOV_long_residuals<-residuals(C20_h12_AOV_long)
hist(C20_h12_AOV_long_residuals) # Slightly right skewed residuals

# Obtain estimated marginal means
emm_C20_h12_AOV_long <- emmeans(C20_h12_AOV_long, ~ combined)

# Perform pairwise comparisons
contrast_C20_h12_AOV_long <- pairs(emm_C20_h12_AOV_long)

# Print the results
print(contrast_C20_h12_AOV_long)

# contrast                                 estimate   SE  df t.ratio p.value
# Control - Old Growth                        -5.20 2.59 165  -2.008  0.1891
# Control - Planted                           -3.01 2.44 165  -1.232  0.6072
# Control - Planted with lianas removed       -9.05 3.04 165  -2.974  0.0176
# Old Growth - Planted                         2.19 1.44 165   1.523  0.4263
# Old Growth - Planted with lianas removed    -3.85 2.32 165  -1.660  0.3484
# Planted - Planted with lianas removed       -6.04 2.15 165  -2.808  0.0283

# Significant differences are observed between:
# Control - Planted with lianas removed       -9.05 3.04 165  -2.974  0.0176
# Planted - Planted with lianas removed       -6.04 2.15 165  -2.808  0.0283

# Non-parametric alternative (not reported in manuscript)
KWtest3 <- kruskal.test(Final_Data$Cover_20_2020_perc~Final_Data$combined)
KWtest3  # p-value = 0.002693

# Perform pairwise Wilcoxon rank-sum tests with p-value adjustments
pairwise.wilcox.test(Final_Data$Cover_20_2020_perc, Final_Data$combined, p.adjust.method = "bonferroni")

# Significant differences are observed between the same as in parametric test 

#--------------------#
###### Plotting ######
#--------------------#

#////////////////#
### Short-term ###
#////////////////#

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
par(mfrow=c(2,2),mar=c(8,6,2,2),las=1,xpd=T)

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

#-----------------------#
# Save the updated data #
#-----------------------#

save(Final_Data, file ='Final_DataV2.rda')
write.csv(Final_Data,"Final_DataV2.csv",row.names = F)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##### HYPOTHESIS 3: ACD ~ Planting + HALP ##### 

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Create new data by selecting specific columns and excluding entries where Site is "Danum"
SBE<-subset(Final_Data, Site == "SBE")

#---------------------------#
###### Model Selection ######
#---------------------------#
# ACD
model1 <- lm(ACD_change ~ combined + HALP, data = SBE)
model2 <- lm(ACD_change ~ combined * HALP, data = SBE)
model3 <- lm(ACD_change ~ HALP, data = SBE)

models <- list(model1, model2, model3)
aictab(cand.set = models, modnames = c("Model1", "Model2", "Model3"))

# New
model1<-lm(New_gap_volume_yr~combined+HALP, data=SBE)
model2<-lm(New_gap_volume_yr~combined*HALP, data=SBE)
model3<-lm(New_gap_volume_yr~HALP, data=SBE)

models <- list(model1, model2, model3)
aictab(cand.set = models, modnames = c("Model1", "Model2", "Model3"))

# Existing
model1<-lm(Existing_proportion_closed_perc_yr~combined+HALP, data=SBE)
model2<-lm(Existing_proportion_closed_perc_yr~combined*HALP, data=SBE)
model3<-lm(Existing_proportion_closed_perc_yr~HALP, data=SBE)

models <- list(model1, model2, model3)
aictab(cand.set = models, modnames = c("Model1", "Model2", "Model3"))

# Model1 is the best-fitting model based on the AICc criteria

#----------------------#
###### Short term ######
#----------------------#

SBE$combined<-factor(SBE$combined)
levels(SBE$combined)

# ACD
ACD_HALP<-lm(ACD_change~combined+HALP, data=SBE)
summary(ACD_HALP) # HALP important for ACD in short term

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          1.298094   0.494814   2.623  0.00984 ** 
#   combinedPlanted                      0.111133   0.418942   0.265  0.79126    
# combinedPlanted with lianas removed  2.715243   0.524745   5.174 9.29e-07 ***
#   HALP                                -0.013058   0.002763  -4.726 6.28e-06 ***

# Check assumptions
plot(ACD_HALP)

coef(ACD_HALP)[[4]]*(max(SBE$HALP,na.rm=T)-min(SBE$HALP,na.rm=T))   # provides a measure of the effect size of HALP in the context of the given data.
# ACD change on ridge is 2.5 C/ha/yr lower than near rivers

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##### HYPOTHESIS 3: Gap Volume~ Planting + HALP ##### 

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#----------------------#
###### Short term ######
#----------------------#

#///////#
## NEW ##
#///////#

GV_HALP_short_NEW<-lm(New_gap_volume_yr~combined+HALP, data=SBE)
summary(GV_HALP_short_NEW) # HALP not important for GV in short term

# Estimate Std. Error t value Pr(>|t|)    
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          498.3181    55.2869   9.013 3.81e-15 ***
#   combinedPlanted                      -36.8666    46.8095  -0.788    0.432    
# combinedPlanted with lianas removed -271.0211    58.6312  -4.622 9.63e-06 ***
#   HALP                                  -0.3482     0.3087  -1.128    0.262                                    -1.393      1.235  -1.128    0.262    

coef(GV_HALP_short_NEW)[['HALP']]*(max(Final_Data$HALP,na.rm=T)-min(Final_Data$HALP,na.rm=T))   
# new gap vol on ridge -67 m3 less than near rivers (fewer new gaps formed near ridges than at valleys)

# Check Assumptions
#plot(GV_HALP_short_NEW) # qq okay
GV_HALP_short_NEW_resid <- residuals(GV_HALP_short_NEW)
hist(GV_HALP_short_NEW_resid) # normal

# Assumptions met

# Log transform the data for new gap volumes
SBE$log_New_gap_volume_yr <- log(SBE$New_gap_volume_yr)  

# Fit the model using log-transformed data
GV_HALP_short_NEW_log <- lm(log_New_gap_volume_yr ~ combined + HALP, data = SBE)

#////////////#
## Existing ##
#////////////#

GV_HALP_short_Existing<-lm(Existing_proportion_closed_perc_yr~combined+HALP, data=SBE)
summary(GV_HALP_short_Existing) # HALP  important for GV in short term

# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         5.845846   0.494556  11.820  < 2e-16 ***
#   combinedPlanted                     0.719747   0.418724   1.719   0.0882 .  
# combinedPlanted with lianas removed 3.038754   0.524471   5.794 5.65e-08 ***
#   HALP                                0.006380   0.002762   2.310   0.0226 *   

coef(GV_HALP_short_Existing)[[4]]*(max(SBE$HALP,na.rm=T)-min(SBE$HALP,na.rm=T))   
# existing gaps on ridges close 1.223786%/yr more than near rivers

# Obtain estimated marginal means
emmGV_HALP_short_Existing <- emmeans(GV_HALP_short_Existing, ~ combined)

# Perform pairwise comparisons
contrast_GV_HALP_short_Existing <- pairs(emmGV_HALP_short_Existing)

# Print the results
print(contrast_GV_HALP_short_Existing)

# contrast                              estimate    SE  df t.ratio p.value
# Control - Planted                        -0.72 0.419 120  -1.719  0.2023
# Control - Planted with lianas removed    -3.04 0.524 120  -5.794  <.0001
# Planted - Planted with lianas removed    -2.32 0.369 120  -6.280  <.0001

#----------------------------#
### Predicting from model  ###
#----------------------------#

# Predict changes
new.data <- data.frame(
  combined = rep(factor(c("Control", "Planted", "Planted with lianas removed")), each = 100),
  HALP = rep(seq(min(SBE$HALP, na.rm = TRUE), max(SBE$HALP, na.rm = TRUE), len = 100), times = 3)
)

levels(new.data$combined)

new.data$ACD_Change<- predict(ACD_HALP, new.data)
new.data$New_Gap_volume_change_hect<- predict(GV_HALP_short_NEW, new.data)
new.data$Existing_proportion_closed_perc_yr<- predict(GV_HALP_short_Existing, new.data)

# Predict on the log-transformed data
new.data$log_New_Gap_volume_change_hect <- predict(GV_HALP_short_NEW_log, new.data)

# Back transform the predicted values to the original scale
new.data$New_Gap_volume_change_hect <- exp(new.data$log_New_Gap_volume_change_hect)

#--------------------#
###### Plotting ######
#--------------------#

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