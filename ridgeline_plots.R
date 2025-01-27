# Synthesised ridgeline plot

# Analyse functionally stable or improving group (FSIG) vs functional loss group (FLG) multiverse outputs
#
# Max Korbmacher, 14 Nov 2024
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
library(ggridges)
library(viridis)
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

# Second, we can remove mis-specified models for the same display
paste("Now, we look only at well-powered models.")
pwr = read.csv("/Users/max/Documents/Local/MS/results/power.csv")
pwr$Formula = gsub("FLG ~ ","",pwr$Formula)
model_info = merge(model_info,pwr,by="Formula")
paste("Convergig models for power calculations: ", length(na.omit(pwr$pwr)), sep="")
paste("Median power: ", 100*median(na.omit(pwr$pwr)),"±",100*mad(na.omit(pwr$pwr)), "%.", sep="")

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
l$Names = gsub("CH3L.1..mg.ml..mean", "Chitinase-3 like-protein-1 mg/ml", l$Names)
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
l$Names
l$group = c("General", "Patient-reported outcome measures","General","Brain markers",
            "Omics", "Clinical markers", "Patient-reported outcome measures",
            "Omics", "Brain markers", "Brain markers",
            "Patient-reported outcome measures", "Omics",
            "Intervention","Clinical markers","Patient-reported outcome measures",
            "Clinical markers",
            "General","General","Patient-reported outcome measures",
            "Patient-reported outcome measures", "Omics",
            "Omics", "Omics")











coef_table$Names = gsub("Age_BL_OFAMS", "General: Age", coef_table$Names)
coef_table$Names = gsub("baselineV", "Brain markers: Lesion volume", coef_table$Names)
coef_table$Names = gsub("baselineC", "Brain markers: Lesion count", coef_table$Names)
coef_table$Names = gsub("TotalVol", "Brain markers: Brain volume", coef_table$Names)
coef_table$Names = gsub("TIV", "Brain markers: Intracranial volume", coef_table$Names)
coef_table$Names = gsub("geno", "Omics: HLA-DRB1 carrier", coef_table$Names)
coef_table$Names = gsub("CH3L.1..mg.ml..mean", "Omics: Chitinase-3 like-protein-1 mg/ml", coef_table$Names)
coef_table$Names = gsub("BL_BMI", "General: Body Mass Index", coef_table$Names)
coef_table$Names = gsub("edss_baseline", "Clinical markers: EDSS", coef_table$Names)
coef_table$Names = gsub("Vit_A_0", "Omics: Vitamin A umol/L", coef_table$Names)
coef_table$Names = gsub("Vit_D_0", "Omics: Vitamin D nmol/L", coef_table$Names)
coef_table$Names = gsub("Vit_E_0", "Omics: Vitamin E umol/L", coef_table$Names)
coef_table$Names = gsub("PF", "PROMs: Physical functioning", coef_table$Names)
coef_table$Names = gsub("RF", "PROMs: Role-physical", coef_table$Names)
coef_table$Names = gsub("BP", "PROMs: Bodily pain", coef_table$Names)
coef_table$Names = gsub("GH", "PROMs: General health", coef_table$Names)
coef_table$Names = gsub("VT", "PROMs: Vitality", coef_table$Names)
coef_table$Names = gsub("SF", "PROMs: Social functioning", coef_table$Names)
#coef_table$Names = gsub("RE", "Role-emotional", coef_table$Names)
coef_table$Names = gsub("MH", "PROMs: Mental health", coef_table$Names)
coef_table$Names = gsub("BAG_c", "Brain markers: Brain age gap", coef_table$Names)
coef_table$Names = gsub("NfL..pg.ml.", "Omics: Neurofilament light chain pg/ml", coef_table$Names)
coef_table$Names = gsub("relapses_12mnths_before_baseline", "Clinical markers: Relapses at baseline", coef_table$Names) # 12 months prior
coef_table$Names = gsub("smoking_OFAMS", "General: Smoking status", coef_table$Names)
coef_table$Names = gsub("sex", "General: Sex (female)", coef_table$Names)
coef_table$Names = gsub("age","General: Age",coef_table$Names)
coef_table$Names = gsub("edss","Clinical markers: Expanded Disability Status Scale (EDSS)",coef_table$Names)
coef_table$Names = gsub("PASAT","Clinical markers: Paced Auditory Serial Addition Test (PASAT)",coef_table$Names)
coef_table$Names = gsub("Omega3_suppl","Intervention: Omega 3 supplement received",coef_table$Names)
coef_table=coef_table[order(coef_table$Names),] # order by name
coef_table = coef_table%>%filter(Names != "(Intercept)")


p1 = coef_table %>% group_by(Names)%>% 
  #mutate(Names = fct_reorder(Names)) %>%
  ggplot(aes(x = (Beta), y = Names, fill = ..x..)) +
  geom_density_ridges_gradient(rel_min_height = 0.01) + xlim(-3,3)+
  #scale_fill_viridis(name = "OR") +
  labs(title = 'Functional loss group membership') +
  theme_bw() + coord_cartesian(clip = "off") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) + xlab("Log odds ratio")+ ylab("")# + xlim(rev(levels(factor(coef_table$Names))))

#
#
#
#
#
# That was the first part. Now, we do the second





# Load data
# model info
model_info = read.csv("/Users/max/Documents/Local/MS/results/PASAT_model_info_optimised.csv")
model_info = model_info %>% dplyr::filter(Formula != "~") %>% dplyr::filter(Formula != "FLG")
# parameter info
coef_table = read.csv("/Users/max/Documents/Local/MS/results/PASAT_stratification_optimised.csv")


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
l$Names = gsub("CH3L.1..mg.ml..mean", "Chitinase-3 like-protein-1 mg/ml", l$Names)
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
l$Names
l$group = c("General", "Patient-reported outcome measures","General","Brain markers",
            "Omics", "Clinical markers", "Patient-reported outcome measures",
            "Omics", "Brain markers", "Brain markers",
            "Patient-reported outcome measures", "Omics",
            "Intervention","Clinical markers","Patient-reported outcome measures",
            "Clinical markers",
            "General","General","Patient-reported outcome measures",
            "Patient-reported outcome measures", "Omics",
            "Omics", "Omics")











coef_table$Names = gsub("Age_BL_OFAMS", "General: Age", coef_table$Names)
coef_table$Names = gsub("baselineV", "Brain markers: Lesion volume", coef_table$Names)
coef_table$Names = gsub("baselineC", "Brain markers: Lesion count", coef_table$Names)
coef_table$Names = gsub("TotalVol", "Brain markers: Brain volume", coef_table$Names)
coef_table$Names = gsub("TIV", "Brain markers: Intracranial volume", coef_table$Names)
coef_table$Names = gsub("geno", "Omics: HLA-DRB1 carrier", coef_table$Names)
coef_table$Names = gsub("CH3L.1..mg.ml..mean", "Omics: Chitinase-3 like-protein-1 mg/ml", coef_table$Names)
coef_table$Names = gsub("BL_BMI", "General: Body Mass Index", coef_table$Names)
coef_table$Names = gsub("edss_baseline", "Clinical markers: EDSS", coef_table$Names)
coef_table$Names = gsub("Vit_A_0", "Omics: Vitamin A umol/L", coef_table$Names)
coef_table$Names = gsub("Vit_D_0", "Omics: Vitamin D nmol/L", coef_table$Names)
coef_table$Names = gsub("Vit_E_0", "Omics: Vitamin E umol/L", coef_table$Names)
coef_table$Names = gsub("PF", "PROMs: Physical functioning", coef_table$Names)
coef_table$Names = gsub("RF", "PROMs: Role-physical", coef_table$Names)
coef_table$Names = gsub("BP", "PROMs: Bodily pain", coef_table$Names)
coef_table$Names = gsub("GH", "PROMs: General health", coef_table$Names)
coef_table$Names = gsub("VT", "PROMs: Vitality", coef_table$Names)
coef_table$Names = gsub("SF", "PROMs: Social functioning", coef_table$Names)
#coef_table$Names = gsub("RE", "Role-emotional", coef_table$Names)
coef_table$Names = gsub("MH", "PROMs: Mental health", coef_table$Names)
coef_table$Names = gsub("BAG_c", "Brain markers: Brain age gap", coef_table$Names)
coef_table$Names = gsub("NfL..pg.ml.", "Omics: Neurofilament light chain pg/ml", coef_table$Names)
coef_table$Names = gsub("relapses_12mnths_before_baseline", "Clinical markers: Relapses at baseline", coef_table$Names) # 12 months prior
coef_table$Names = gsub("smoking_OFAMS", "General: Smoking status", coef_table$Names)
coef_table$Names = gsub("sex", "General: Sex (female)", coef_table$Names)
coef_table$Names = gsub("age","General: Age",coef_table$Names)
coef_table$Names = gsub("edss","Clinical markers: Expanded Disability Status Scale (EDSS)",coef_table$Names)
coef_table$Names = gsub("PASAT","Clinical markers: Paced Auditory Serial Addition Test (PASAT)",coef_table$Names)
coef_table$Names = gsub("Omega3_suppl","Intervention: Omega 3 supplement received",coef_table$Names)
coef_table=coef_table[order(coef_table$Names),] # order by name
coef_table = coef_table%>%filter(Names != "(Intercept)")

p2 = coef_table %>% group_by(Names)%>% 
  #mutate(Names = fct_reorder(Names)) %>%
  ggplot(aes(x = (Beta), y = Names, fill = ..x..)) +
  geom_density_ridges_gradient() + xlim(-3,3)+
  #scale_fill_viridis(name = "OR") +
  labs(title = 'Cognitive decline group membership') +
  theme_bw() + coord_cartesian(clip = "off") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) + xlab("Log odds ratio")+ ylab("")# + xlim(rev(levels(factor(coef_table$Names))))
plot = ggarrange(p1,p2)

# plot=ggforestplot::forestplot(
#   df = l,
#   name = Names,
#   estimate = Beta_Md,
#   se = Beta_MAD,
#   xlab="Median odds ratio ± median absolute deviation",
#   logodds = TRUE,
#   #pvalue = P_Md,
#   psignif = .05
# )+
#   ggforce::facet_col(
#     facets = ~group,
#     scales = "free_y",
#     space = "free"
#   )
# #plot=annotate_figure(plot,top = "Predictors of functional loss group membership")
ggsave(plot = plot, filename = "/Users/max/Documents/Local/MS/results/ridgeline.pdf", width = 15, height = 6)
ggsave(plot = plot, filename = "/Users/max/Documents/Local/MS/results/ridgeline.png", width = 15, height = 6)

#

