# Factor analysis for redundancy checks
# Max Korbmacher, May 2025
library(psych)
library(dplyr)
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
# filter data for continous vars
data = data %>% select("relapses_12mnths_before_baseline", 
                "edss", "PASAT", "BL_BMI", 
                "BAG_c", "baselineC", "baselineV",
                "Vit_A_0","Vit_D_0", "Vit_E_0")
#                "PF", "BP", "GH", "VT", "SF","MH")
# Perform factor analysis 
factor_analysis <- factanal(na.omit(data),factors = 3)
print(factor_analysis)

# check only for the QoL data
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
data = data %>% select("PF", "BP", "GH", "VT", "SF","MH")
factor_analysis <- factanal(na.omit(data),factors = 3)
# problem with singularity!
# this problem comes down to the coding used when grouping the items.

# we also check whether the original items can be 