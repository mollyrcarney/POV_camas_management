# Camas Management on the North American Columbia Plateau
This repository contains all data and code for the forthcoming paper "Harvesting strategies as evidence for 4,000 years of camas (*Camassia quamash*) management in the North American Columbia Plateau." 

Files included:

-Dataset S1. Radiocarbon dates used in this analysis.

-Dataset S2. Master spreadsheet with archaeological bulb data and provenience.

-Dataset S3. Modern, comparative *Camassia quamash* data 3.

-Dataset S4. δ18O data from Cleland lake, British Columbia, used in Figure 4 (Steinman et al. 2016).

-Dataset S5. Northeastern Pacific sea surface temperature (°C) reconstruction used in Figure 4 (Kienast and McKay 2001).

### References

Carney, M. et al. Harvesting strategies as evidence for 4,000 years of camas (Camassia quamash) management in the North American Columbia Plateau. *SocArXiv* November 23. doi:10.31235/osf.io/c752m (2020).

Kienast, S. S. & McKay, J. L. Sea surface temperatures in the subarctic northeast Pacific reflect millennial‐scale climate oscillations during the last 16 kyrs. *Geophysical Research Letters* 28, 1563-1566 (2001).

Steinman, B. A. et al. Oxygen isotope records of Holocene climate variability in the Pacific Northwest. *Quaternary Science Reviews* 142, 40-60, doi:10.1016/j.quascirev.2016.04.012 (2016).


## R code
### Packages
```r
library(dplyr)
library(resphape2)
library(DescTools)
library(rcompanion)
```

### Load data
```r
archbulbs.analysis <- read.csv(file.choose(), fileEncoding="UTF-8-BOM") #SI1_cvap_master
control.camas <- read.csv(file.choose(), fileEncoding="UTF-8-BOM") #SI2_control_camas
cleland <- read.csv(file.choose(), fileEncoding="UTF-8-BOM") #SI4_cleland_lake_d18o
sst <- read.csv(file.choose(), fileEncoding="UTF-8-BOM") #SI5_sea_surface_temp
```

### Statistical analyses for bulb maturity and comparison with archaelogical materials
```r
#Test if immature vs mature are from same sample populations
wilcox.test(NumLeaves ~ DateRange, data = control.camas)

#Test if charred ratios differ between control sample
wilcox.test(CharredRatio ~ DateRange, data = control.camas)

#Prep data for main graphs and analyses
subset_analysis <- archbulbs.analysis %>%
        select(NumLeaves, DateRange, Group)

subset_control <- control.camas %>%
        select(NumLeaves, DateRange, Group)

boxplotdf <- rbind(subset_analysis, subset_control)
boxplotdf <- na.omit(boxplotdf)

#Pairwise Mann-Whitney tests for comparing number of leaves across time/controls
pairwise.wilcox.test(boxplotdf$NumLeaves, boxplotdf$DateRange, p.adjust.method="BH")

#test for change in bulb size over time
kruskal.test(Ratio ~ DateRange, data = archbulbs.analysis) 
epsilonSquared(x = archbulbs.analysis$Ratio, g = archbulbs.analysis$DateRange)
kruskal.test(Weight ~ DateRange, data = archbulbs.analysis)
epsilonSquared(x = archbulbs.analysis$Weight, g = archbulbs.analysis$DateRange)
```

### Figures 3 and 4
```r
#Vector for visualizing leaf scale data by harvesting strategy
colorbygroup <- c("#0072b2", "#0072b2", "#cc79a7", "#0072b2", 
                  "#cc79a7", "#0072b2", "#0072b2", "#cc79a7", "#D3D3D3", "#D3D3D3")

#Figure 3 Boxplot of leaf scale numbers through time
plot.new()
boxplot(NumLeaves ~ DateRange, data = boxplotdf,
        legend("bottomleft", legend = c("Selective Harvesting", "Stripping", "Control"),
                col = c("#0072b2", "#cc79a7", "#D3D3D3"),
                bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE),
        las = 2,
        col = colorbygroup,
        ylab = "Number of Leaves", xlab = "",)
       
#Figure 4 Climate and bulb data graph

#Plot Figure 4 graph
dev.off()
par(mar = c(4, 4, 1, 1))
par(mfrow=c(4,1))
plot(cleland$age_calBP, cleland$d18OcarbVPDB, 
     type= "l", col= "#0072b2", pch=19, lwd = 1.5,
     xlim=c(0,4500), ylim=c(-1, -16),
     xlab = "", ylab="??18O")
ticks = c(500, 1000, 1500, 2500, 3000, 3500)
plot(sst$calyrBP, sst$sst.uk37, 
     type= "l", col= "#71C4EA", pch=19, lwd = 2,
     ylab="SST (°C)", 
     xlab = "", xlim=c(0, 4500), ylim= c(8, 11))
plot(archbulbs.analysis$MedianDate, archbulbs.analysis$Ratio, 
     type= "p", col="#7187EA", pch=19,   xlim=c(0,4500),
     xlab="", ylab = "Size Ratio")
plot(archbulbs.analysis$MedianDate, archbulbs.analysis$NumLeaves, 
     type = "p", col="#EA7187", xlab="Age BP", 
     ylab= "Num Leaves (n)", pch=19, 
     ylim = c(2.8, 4.2), xlim=c(0,4500))
```
