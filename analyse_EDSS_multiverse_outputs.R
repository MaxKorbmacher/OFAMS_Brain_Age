# Analyse functionally stable or improving group (FSIG) vs functional loss group (FLG) multiverse outputs
#
# Max Korbmacher, May 2025
#
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
#
# load stuff
library(dplyr)
library(sjstats)
#devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
library(ggpubr)
library(ggforce)
#
# Load data
# model info
model_info = read.csv("/Users/max/Documents/Local/MS/results/EDSS_model_info_optimised.csv")
model_info = model_info %>% dplyr::filter(Formula != "~") %>% dplyr::filter(Formula != "FLG")
# parameter info
coef_table = read.csv("/Users/max/Documents/Local/MS/results/EDSS_stratification_optimised.csv")
#coef_table = read.csv("/Users/max/Documents/Local/MS/results/model_coefficients.csv")
#
#
# First, some general information from the model_info object
paste("Recap. The number of models run was ", nrow(na.omit(model_info)), " with ", round(sum(ifelse(model_info$McFaddenR2>.2,1,0))/nrow(model_info)*100,2), "% showing an excellent model fit.", sep="")
paste("Good model fit (R2>10%): ", round(sum(ifelse(model_info$McFaddenR2>.1,1,0))/nrow(model_info)*100,2), "%.", sep="")
paste("Median McFadden pseudo R2 = ",round(median(model_info$McFaddenR2),4),"±",round(mad(model_info$McFaddenR2),2),sep="")
paste("AUC =",round(median(model_info$AUC),2),"±",round(mad(model_info$AUC),2))
paste("Brier =",round(median(model_info$Brier),2),"±",round(mad(model_info$Brier),2))
# Second, we can remove mis-specified models for the same display
paste("Now, we look only at well-powered models.")
pwr = read.csv("/Users/max/Documents/Local/MS/results/power.csv")
#pwr$Formula = gsub("FLG ~ ","",pwr$Formula)
model_info = merge(model_info,pwr,by="Formula")
paste("Convergig models for power calculations: ", length(na.omit(model_info$pwr)), sep="")
paste("Median power: ", 100*median(na.omit(model_info$pwr)),"±",100*mad(na.omit(model_info$pwr)), "%.", sep="")
#
#
#
write.csv(x = model_info, "/Users/max/Documents/Local/MS/results/model_info.csv")
#
# Standardize formula column
coef_table = coef_table %>% filter(Names != "(Intercept)")
coef_table$Formula = gsub('[^[:alnum:] ]', '', coef_table$Formula)
coef_table$Formula = gsub('FLGc ', '',coef_table$Formula)
coef_table$Formula = gsub('  ', ' ',coef_table$Formula)
coef_table$Formula = gsub(' ', ' + ',coef_table$Formula)
model_info$Formula = gsub('[^[:alnum:] ]', '', model_info$Formula)
model_info$Formula = gsub('  ', ' + ',model_info$Formula)
coef_table = full_join(coef_table,model_info,by="Formula")
# # Duplicate rows and extract predictor names
# model_info$Formula = as.character(c(model_info$Formula))
# # Ensure Formula column is character type and then split by space
# model_info <- model_info %>%
#   mutate(Formula = as.character(Formula)) %>%  # Ensure Formula is treated as character
#   tidyr::separate_rows(Formula, sep = " ") %>%        # Split Formula by space
#   rename(Names = Formula)                  # Rename Formula to Predictors
# # View the result
# model_info = model_info %>% filter(Names != "+")
# 
# 
# coef_table = cbind(model_info, coef_table %>% dplyr::select(-Names))
rm(model_info)
#
#
#
# These are the stats looking at all model outputs TOGETHER
# CONVERGANCE PROBLEMS ARE IGNORED BY THIS APPROACH!!
#
l = coef_table %>% group_by(Names)%>%dplyr::summarize(Beta_Md=median((Beta)), Beta_MAD = mad((Beta)), P_Md=median(P)) # estimate medians and MADs
coef_table$bigger = ifelse(exp(coef_table$Beta)>1,1,0)
bigger = coef_table%>%group_by(Names)%>%summarize("OR>1" = sum(bigger)/length(bigger))
l = merge(l,bigger,by="Names") # merge data keeping the correct order
l$Names = gsub("Age_BL_OFAMS", "Age", l$Names)
l$Names = gsub("baselineV", "Lesion volume", l$Names)
l$Names = gsub("baselineC", "Lesion count", l$Names)
l$Names = gsub("TotalVol", "Brain volume", l$Names)
l$Names = gsub("TIV", "Intracranial volume", l$Names)
l$Names = gsub("geno", "HLA-DRB1 carrier", l$Names)
#l$Names = gsub("CH3L.1..mg.ml..mean", "Chitinase-3 like-protein-1 mg/ml", l$Names)
l$Names = gsub("BL_BMI", "Body Mass Index", l$Names)
l$Names = gsub("edss_baseline", "EDSS", l$Names)
l$Names = gsub("Vit_A_0", "Vitamin A umol/L", l$Names)
l$Names = gsub("Vit_D_0", "Vitamin D nmol/L", l$Names)
l$Names = gsub("Vit_E_0", "Vitamin E umol/L", l$Names)
l$Names = gsub("PF", "Physical functioning", l$Names)
l$Names = gsub("RF", "Role-physical", l$Names)
l$Names = gsub("BP", "Bodily pain", l$Names)
l$Names = gsub("GH", "General health", l$Names)
l$Names = gsub("VT", "Vitality", l$Names)
l$Names = gsub("SF", "Social functioning", l$Names)
l$Names = gsub("Current_DMT", "Disease modifying treatment", l$Names)
#l$Names = gsub("RE", "Role-emotional", l$Names)
l$Names = gsub("MH", "Mental health", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("NfL..pg.ml.", "Neurofilament light chain pg/ml", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("relapses_12mnths_before_baseline", "Relapses at baseline", l$Names) # 12 months prior
l$Names = gsub("smoking_OFAMS", "Smoking status", l$Names)
l$Names = gsub("sex", "Sex (female)", l$Names)
#l$Names = gsub("Treatment_OFAMS", "Treatment", l$Names)
l$Names = gsub("age","Age",l$Names)
l$Names = gsub("edss","Expanded Disability Status Scale (EDSS)",l$Names)
l$Names = gsub("PASAT","Paced Auditory Serial Addition Test (PASAT)",l$Names)
l$Names = gsub("Omega3_suppl","Omega 3 supplement received",l$Names)
l=l[order(l$Names),] # order by name
l = l%>%filter(Names != "(Intercept)")
l$group = c("General", "Patient-reported outcome measures","General","Brain markers",
  "Intervention", "Clinical markers",
  "Omics", "Brain markers", "Brain markers",
  "Patient-reported outcome measures", "Omics",
  "Intervention","Clinical markers","Patient-reported outcome measures",
  "Clinical markers",
  "General","General",
  "Patient-reported outcome measures", "Omics",
  "Omics", "Omics")
plot=ggforestplot::forestplot(
  df = l,
  name = Names,
  estimate = Beta_Md,
  se = Beta_MAD,
  xlab="Median odds ratio ± median absolute deviation",
  logodds = TRUE,
  pvalue = P_Md,
  psignif = .05
)+
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
#plot=annotate_figure(plot,top = "Predictors of functional loss group membership")
ggsave(plot = plot, filename = "/Users/max/Documents/Local/MS/results/EDSS_multiverse.pdf", width = 8, height = 8)
#
# make Odds ratio table
mktable = function(coef_table){
l = coef_table %>% group_by(Names)%>%dplyr::summarize(Beta_Md=round(median(exp(Beta)),4), Beta_MAD = round(mad(exp(Beta)),4), P_Md=round(median(P),4)) # estimate medians and MADs
coef_table$bigger = ifelse(exp(coef_table$Beta)>1,1,0)
bigger = coef_table%>%group_by(Names)%>%summarize("OR>1" = round(sum(bigger)/length(bigger),4))
l = merge(l,bigger,by="Names") # merge data keeping the correct order
l$Names = gsub("Age_BL_OFAMS", "Age", l$Names)
l$Names = gsub("baselineV", "Lesion volume", l$Names)
l$Names = gsub("baselineC", "Lesion count", l$Names)
l$Names = gsub("TotalVol", "Brain volume", l$Names)
l$Names = gsub("TIV", "Intracranial volume", l$Names)
l$Names = gsub("geno", "HLA-DRB1 carrier", l$Names)
l$Names = gsub("Current_DMT", "Disease modifying treatment", l$Names)
#l$Names = gsub("CH3L.1..mg.ml..mean", "Chitinase-3 like-protein-1 mg/ml", l$Names)
l$Names = gsub("BL_BMI", "Body Mass Index", l$Names)
l$Names = gsub("edss_baseline", "EDSS", l$Names)
l$Names = gsub("Vit_A_0", "Vitamin A umol/L", l$Names)
l$Names = gsub("Vit_D_0", "Vitamin D nmol/L", l$Names)
l$Names = gsub("Vit_E_0", "Vitamin E umol/L", l$Names)
l$Names = gsub("PF", "Physical functioning", l$Names)
l$Names = gsub("RF", "Role-physical", l$Names)
l$Names = gsub("BP", "Bodily pain", l$Names)
l$Names = gsub("GH", "General health", l$Names)
l$Names = gsub("VT", "Vitality", l$Names)
l$Names = gsub("SF", "Social functioning", l$Names)
l$Names = gsub("RE", "Role-emotional", l$Names)
l$Names = gsub("MH", "Mental health", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("NfL..pg.ml.", "Neurofilament light chain pg/ml", l$Names)
l$Names = gsub("BAG_c", "Brain age gap", l$Names)
l$Names = gsub("relapses_12mnths_before_baseline", "Relapses at baseline", l$Names) # 12 months prior
l$Names = gsub("smoking_OFAMS", "Smoking status", l$Names)
l$Names = gsub("sex", "Sex (female)", l$Names)
#l$Names = gsub("Treatment_OFAMS", "Treatment", l$Names)
l$Names = gsub("age","Age",l$Names)
l$Names = gsub("edss","Expanded Disability Status Scale",l$Names)
l$Names = gsub("PASAT","Paced Auditory Serial Addition Test",l$Names)
l$Names = gsub("Omega3_suppl","Omega 3 supplement received",l$Names)
l=l[order(l$Names),] # order by name
l = l%>%dplyr::filter(Names != "(Intercept)")
l$group = c("General", "Patient-reported outcome measures","General","Brain markers",
            "Intervention", "Clinical markers",
            "Omics", "Brain markers", "Brain markers",
            "Patient-reported outcome measures", "Omics",
            "Intervention","Clinical markers","Patient-reported outcome measures",
            "Clinical markers",
            "General","General",
            "Patient-reported outcome measures", "Omics",
            "Omics", "Omics")
return(l)
}
l = mktable(coef_table)
write.csv(x = l,"/Users/max/Documents/Local/MS/results/EDSS_multiverse.csv")
#
# # Consider only well-powered models ####
# coef_table$Power = coef_table$pwr
# cor(coef_table$Power,coef_table$McFaddenR2,use = "pairwise.complete.obs") # just for my own interest: R2 and power correlation
# # note, there is no correspondence when considering low-power outcomes
# c1 = coef_table %>% dplyr::filter(Power >= 0.8) # filter for power ≥ 80%
# cor(c1$Power,c1$McFaddenR2,use = "pairwise.complete.obs") # Now, there is a better correspondence between R2 and power
# l = mktable(c1)
# write.csv(x = l,"/Users/max/Documents/Local/MS/results/FDG_multiverse_well_powered.csv")
# 
# # Consider only high R2 models ####
# c1 = coef_table %>% dplyr::filter(McFaddenR2>.2)
# cor(c1$Power,c1$McFaddenR2,use = "pairwise.complete.obs") # This time a negative correlation between R2 and power?!
# l = mktable(c1)
# write.csv(x = l,"/Users/max/Documents/Local/MS/results/FDG_multiverse_high_R2.csv")
# 
# # Consider only well-powered and high R2 models ####
# c1 = c1 %>% dplyr::filter(Power>.8)
# cor(c1$Power,c1$McFaddenR2,use = "pairwise.complete.obs") # WORSE correspondence between R2 and power
# l = mktable(c1)
# write.csv(x = l,"/Users/max/Documents/Local/MS/results/FDG_multiverse_well_powered_and_high_R2.csv")
