# Analyse multiverse outputs
#
# Max Korbmacher, 14 Nov 2024
#
# load stuff
library(dplyr)
library(sjstats)
#devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
library(ggpubr)
library(ggforce)
#
# Then load specific for the two predictions
#
#
#
# model info
model_info = read.csv("/Users/max/Documents/Local/MS/results/EDSS_model_info_optimised.csv")
# parameter info
coef_table = read.csv("/Users/max/Documents/Local/MS/results/EDSS_stratification_optimised.csv")
#
paste("Recap. The number of models run was ", nrow(model_info), sep="")
paste("Median McFadden pseudo R2 = ",round(median(model_info$McFaddenR2),2),"±",round(mad(model_info$McFaddenR2),2),sep="")
#
l = coef_table %>% group_by(Names)%>%summarize(Beta_Md=median((Beta)), Beta_MAD = mad((Beta)), P_Md=median(P)) # estimate medians and MADs
a = props(coef_table%>%group_by(Names),(exp(Beta))>1) # estimate proportions of effects towards a single direction.
b = props(coef_table%>%group_by(Names),(exp(Beta))<1)
l = merge(l,a,by="Names") # merge data keeping the correct order
l = merge(l,b,by="Names")
l$Names = gsub("Age_BL_OFAMS", "Age", l$Names)
l$Names = gsub("baselineV", "Lesion volume", l$Names)
l$Names = gsub("baselineC", "Lesion count", l$Names)
l$Names = gsub("TotalVol", "Brain volume", l$Names)
l$Names = gsub("TIV", "Intracranial volume", l$Names)
l$Names = gsub("geno", "HLA-DRB1 carrier", l$Names)
l$Names = gsub("CH3L.1..mg.ml..mean", "CHI3L1 mg/ml", l$Names)
l$Names = gsub("Mean_BMI_OFAMS", "Mean BMI", l$Names)
l$Names = gsub("edss_baseline", "EDSS", l$Names)
l$Names = gsub("Vit_A_0", "Vitamin A mcmol/L", l$Names)
l$Names = gsub("Vit_D_0", "Vitamin D mcmol/L", l$Names)
l$Names = gsub("Vit_E_0", "Vitamin E mcmol/L", l$Names)
l$Names = gsub("PF", "Physical functioning", l$Names)
l$Names = gsub("RF", "Role-physical", l$Names)
l$Names = gsub("BP", "Bodily pain", l$Names)
l$Names = gsub("GH", "General health", l$Names)
l$Names = gsub("VT", "Vitality", l$Names)
l$Names = gsub("SF", "Social functioning", l$Names)
l$Names = gsub("RE", "Role-emotional", l$Names)
l$Names = gsub("MH", "Mental health", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("NfL..pg.ml.", "NfL pg/ml", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("relapses_12mnths_before_baseline", "Relapses 12 months prior baseline", l$Names)
l$Names = gsub("smoking_OFAMS", "Smoking status", l$Names)
l$Names = gsub("sex", "Sex", l$Names)
l$Names = gsub("Treatment_OFAMS", "Treatment", l$Names)
l$Names = gsub("age","Age",l$Names)
l$Names = gsub("edss","EDSS",l$Names)
l$Names = gsub("Omega3_suppl","Omega 3 supplement received",l$Names)
l=l[order(l$Names),] # order by name
l = l%>%filter(Names != "(Intercept)")
l$group = c("General", "Patient-reported outcome measures","Brain markers",
  "Omics", "Clinical markers", "Patient-reported outcome measures",
  "Omics", "Brain markers", "Brain markers",
  "General", "Patient-reported outcome measures", "Omics",
  "Intervention","Clinical markers","Patient-reported outcome measures",
  "Clinical markers","Patient-reported outcome measures","Patient-reported outcome measures",
  "General","General","Patient-reported outcome measures",
  "Intervention", "Patient-reported outcome measures", "Omics",
  "Omics", "Omics")
plot=ggforestplot::forestplot(
  df = l,
  name = Names,
  estimate = Beta_Md,
  se = Beta_MAD,
  xlab="Median odds ratio ± MAD",
  logodds = TRUE,
  pvalue = P_Md,
  psignif = .05
)+
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
plot=annotate_figure(plot,top = "Predictors of functional loss group membership")
ggsave(plot = plot, filename = "/Users/max/Documents/Local/MS/results/FDG_multiverse.pdf", width = 8, height = 8)
# 
# 
# # PASAT
# df=read.csv("/Users/max/Documents/Local/MS/results/PASAT_effects.csv")
# # get estimates
# l=df %>% group_by(Names)%>%summarize(Beta_Md=median(Beta), Beta_MAD = mad(Beta), P_Md=median(P)) # estimate medians and MADs
# a = props(df%>%group_by(Names),Beta>0) # estimate proportions of effects towards a single direction.
# b = props(df%>%group_by(Names),Beta<0)
# l = merge(l,a,by="Names") # merge data keeping the correct order
# l = merge(l,b,by="Names")
# l$Names = gsub("Age_BL_OFAMS", "Age", l$Names)
# l$Names = gsub("baselineV", "Lesion volume", l$Names)
# l$Names = gsub("baselineC", "Lesion count", l$Names)
# l$Names = gsub("TotalVol", "Brain volume", l$Names)
# l$Names = gsub("TIV", "Intracranial volume", l$Names)
# l$Names = gsub("geno", "HLA-DRB1 carrier", l$Names)
# l$Names = gsub("CH3L.1..mg.ml..mean", "CH3L1", l$Names)
# l$Names = gsub("Mean_BMI_OFAMS", "Mean BMI", l$Names)
# l$Names = gsub("edss_baseline", "EDSS", l$Names)
# l$Names = gsub("Vit_A_0", "Vitamin A", l$Names)
# l$Names = gsub("Vit_D_0", "Vitamin D", l$Names)
# l$Names = gsub("Vit_E_0", "Vitamin E", l$Names)
# l$Names = gsub("PF", "PROM: Physical functioning", l$Names)
# l$Names = gsub("RF", "PROM: Role-physical", l$Names)
# l$Names = gsub("BP", "PROM: Bodily pain", l$Names)
# l$Names = gsub("GH", "PROM: General health", l$Names)
# l$Names = gsub("VT", "PROM: Vitality", l$Names)
# l$Names = gsub("SF", "PROM: Social functioning", l$Names)
# l$Names = gsub("RE", "PROM: Role-emotional", l$Names)
# l$Names = gsub("MH", "PROM: Mental health", l$Names)
# l$Names = gsub("BAG_c", "Brain age gap", l$Names)
# l$Names = gsub("NfL..pg.ml.", "NfL pg/ml", l$Names)
# l$Names = gsub("BAG_c", "Brain age gap", l$Names)
# l$Names = gsub("relapses_12mnths_before_baseline", "Relapses 12 months prior baseline", l$Names)
# l$Names = gsub("smoking_OFAMS", "Smoking status", l$Names)
# l$Names = gsub("sex", "Sex", l$Names)
# l$Names = gsub("Treatment_OFAMS", "Received Omega 3 supplement", l$Names)
# l=l[order(l$Names),] # order by name
# plot=ggforestplot::forestplot(
#   df = l,
#   name = Names,
#   estimate = Beta_Md,
#   se = Beta_MAD,
#   xlab="Median beta ± MAD"
# )
# plot=annotate_figure(plot,top = "Predictors of PASAT")
# ggsave(plot = plot, filename = "/Users/max/Documents/Local/MS/results/PASAT_multiverse.pdf", width = 5, height = 6)
# #
# #
# #
# #
# # EDSS
# df2=read.csv("/Users/max/Documents/Local/MS/results/EDSS_effects.csv")
# # get estimates
# edss=df2 %>% group_by(Names)%>%summarize(Beta_Md=median(Beta), Beta_MAD = mad(Beta), P_Md=median(P)) # estimate medians and MADs
# a = props(df2%>%group_by(Names),Beta>0) # estimate proportions of effects towards a single direction.
# b = props(df2%>%group_by(Names),Beta<0)
# edss = merge(edss,a,by="Names") # merge data keeping the correct order
# edss = merge(edss,b,by="Names")
# edss$Names = gsub("Age_BL_OFAMS", "Age", edss$Names)
# edss$Names = gsub("BL_PASATcorrect", "PASAT", edss$Names)
# edss$Names = gsub("baselineV", "Lesion volume", edss$Names)
# edss$Names = gsub("baselineC", "Lesion count", edss$Names)
# edss$Names = gsub("TotalVol", "Brain volume", edss$Names)
# edss$Names = gsub("TIV", "Intracranial volume", edss$Names)
# edss$Names = gsub("geno", "HLA-DRB1 carrier", edss$Names)
# edss$Names = gsub("CH3L.1..mg.ml..mean", "CH3L1", edss$Names)
# edss$Names = gsub("Mean_BMI_OFAMS", "Mean BMI", edss$Names)
# edss$Names = gsub("edss_baseline", "EDSS", edss$Names)
# edss$Names = gsub("Vit_A_0", "Vitamin A", edss$Names)
# edss$Names = gsub("Vit_D_0", "Vitamin D", edss$Names)
# edss$Names = gsub("Vit_E_0", "Vitamin E", edss$Names)
# edss$Names = gsub("PF", "PROM: Physical functioning", edss$Names)
# edss$Names = gsub("RF", "PROM: Role-physical", edss$Names)
# edss$Names = gsub("BP", "PROM: Bodily pain", edss$Names)
# edss$Names = gsub("GH", "PROM: General health", edss$Names)
# edss$Names = gsub("VT", "PROM: Vitality", edss$Names)
# edss$Names = gsub("SF", "PROM: Social functioning", edss$Names)
# edss$Names = gsub("RE", "PROM: Role-emotional", edss$Names)
# edss$Names = gsub("MH", "PROM: Mental health", edss$Names)
# edss$Names = gsub("BAG_c", "Brain age gap", edss$Names)
# edss$Names = gsub("NfL..pg.ml.", "NfL pg/ml", edss$Names)
# edss$Names = gsub("BAG_c", "Brain age gap", edss$Names)
# edss$Names = gsub("relapses_12mnths_before_baseline", "Relapses 12 months prior baseline", edss$Names)
# edss$Names = gsub("smoking_OFAMS", "Smoking status", edss$Names)
# edss$Names = gsub("sex", "Sex", edss$Names)
# edss$Names = gsub("Treatment_OFAMS", "Received Omega 3 supplement", edss$Names)
# edss=edss[order(edss$Names),] # order by name
# plot2=ggforestplot::forestplot(
#   df = edss,
#   name = Names,
#   estimate = Beta_Md,
#   se = Beta_MAD
# )
# plot2=annotate_figure(plot2,top = "Predictors of EDSS")
# ggsave(plot = plot2, filename = "/Users/max/Documents/Local/MS/results/EDSS_multiverse.pdf", width = 5, height = 6)
