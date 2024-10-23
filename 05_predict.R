#!/usr/bin/env Rscript
# This script is for creating both corrected and raw brain age predictions
# Max Korbmacher, 23rd October 2024
noquote("This script is supposed to be run within/with access to the files in the repository at https://github.com/MaxKorbmacher/OFAMS_Brain_Age.")
noquote("A minimum input is the path to your csv file assembled from FreeSurfer recon-all tabular outputs, but you can add output file or repo as well.")
noquote("The script returns the original data frame / csv-file with new columns for predicted age and corrected predicted age.")
noquote("##########")
noquote("Arguments:")
noquote("##########")
noquote("1. Path to data (csv) predictions shall be run on: '/my/path/to/data.csv'")
noquote("2. Output file and location: '/my/output/file.csv'")
noquote("   If not defined, it will be saved at the current directory.")
noquote("3. If not run inside the repository, specify the repository name '/my/path/to/repo/'.")
noquote("##########")
#
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least a data file must be supplied for predictions.", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "predictions.csv"
  args[3] = getwd()
} else if (length(args)==2) {
  args[3] = getwd()
}
noquote("Load dependencies, if not already installed.")
# load dependencies
if (!require("pacman")) install.packages("pacman")
pacman::p_load(mgcv, Metrics, MASS, relaimpo, lmtest, lme4, effectsize, 
               rempsyc, dplyr, ggeffects, ggpubr, marginaleffects,MuMIn,rlist,
               lmerTest,update = F)
noquote("Load data and model.")
df=read.csv(args[1])
m=readRDS(paste0(args[3],"/model.rda"))
params=read.csv(paste0(args[3],"/corr_params.csv"))
noquote("Predict.")
df$brainage=predict(object=m,df)
df$corrected_brainage= df$brainage + (df$age - (params$slope*df$brainage+params$intercept))
noquote("Write data file.")
write.csv(x = df, file = args[2], row.names = F)