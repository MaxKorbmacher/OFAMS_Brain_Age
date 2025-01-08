# OFAMS brain age analysis
# 07 Jan 2025
# Max Korbmacher
#
# Note: this is my messy file for data cleaning/wrangling (and some descriptives) 
# The most important is the produced interrim_data.csv which provides the input for the multiverse analyses.
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# set savepath
savepath = "/Users/max/Documents/Local/MS/results/"
# 0.1 Prep ####
# load pkgs
# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(haven,car,dplyr,lme4,lmerTest,ggeffects,sjPlot,
               data.table,factoextra,simputation,ggpubr,psych,reshape2,
               gplots,rstatix,ggseg,marginaleffects,MuMIn,tidyr,#glmnet,ipflasso,parameters,
               olsrr,parallel,gghighlight,cowplot,zoo,readxl)
#
#
# load data
data = read.csv("/Users/max/Documents/Local/MS/GAM_preds.csv") # contains brain age predictions
edss = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/edss.sav") # contains Expanded Disability Status Scale (EDSS) scores
prom = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/sf36.sas7bdat") # patient reported outcome measures (PROM)
fati = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/fatigue.sas7bdat") # fatigue scores
demo = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/Demographics.sav") # Demographics
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
relapse = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/relapsereport.sas7bdat") # relapse info
geno = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/hla_drb1_15_Juni_2011_alle pas.sav") # genotype 
icv = read.csv("/Users/max/Documents/Local/MS/data/icv.csv")
lesion_count = read.csv("/Users/max/Documents/Local/MS/data/lesion_count.csv")
lesion_vol = read.csv("/Users/max/Documents/Local/MS/data/lesion_vol.csv")
blood = read.csv("/Users/max/Documents/Local/MS/data/Copy of Paired_serum-MRI.csv")
background = read_excel("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/Bakgrunnsdata-OFAMS.xls")
#
#
#
# copy original data frame for later analysis of brain structure 
anat = data
# caluculate BAG
data = data %>% dplyr::select(eid, session, age, brainage, corrected_brainage)
data$BAG_c = data$corrected_brainage-data$age
data$BAG = data$brainage-data$age
#
#
######## ADD COVARIATES!
# start with edss
#edss %>% select(Patnr,contains("EDSS")) %>% melt(id.vars="Patnr")
edss_long = melt(edss, id.vars="Patnr")
names(edss_long) = c("eid","session","edss")
print("Note that we only have EDSS scores for baseline and certain follow-up times.")
print("We can note an increase in EDSS scores over time. (Yet, with increasing variablility.)")
edss_long %>% group_by(session) %>% summarize(M=mean(na.omit(edss)), SD=sd(na.omit(edss))) # show mean+sd for each session
levels(edss_long$session) = c(0,6,12,18,24) # standardize session names
data$session = factor(data$session) # treat session as the factor it is also in the brain age data frame
data = left_join(data,edss_long, by=c("eid","session")) # join edss and brain age data frames.
#
# demographics: demo$Gender # 0 = Male
# genotype: geno$HLA_1501_1 # HLADRB1501 positive = 1
# relapse: relapse$relapseno # nb of relapses > var of interest
demo = full_join(demo, geno, by = "Patnr")
relapse$Patnr = relapse$patno
demo = full_join(demo, relapse, by = "Patnr")
print("This completes the cross-sectional / baseline data in the demo frame.")
# patient recorded outcome measures (prom)
prom %>% dplyr::select(patno,contains("spm")) %>% is.na() %>% colSums() # check number of NAs (none)
#print("There are some (very few) NAs. We impute using sequential hot deck imputation.")
#prom=prom %>% dplyr::select(patno,VISIT,contains("spm"))
#prom=impute_shd(prom, .~1, backend="VIM")
visits = levels(factor(prom$VISIT))
prom = prom %>% rename(eid = patno, session = VISIT)
prom$session=factor(prom$session)
levels(prom$session)=c(0,12,24,6)
data=full_join(data,prom,by=c("eid","session"))
#
# Add sex and genotype to demo data.
geno = geno %>% rename(eid=Patnr,geno=HLA_1501_1)
sex = demo %>% dplyr::select(Patnr,Gender) %>% rename(eid=Patnr,sex=Gender)
data = left_join(data,geno, by="eid")
data = left_join(data,sex, by="eid")
#
#
#
#
# 0.2 Functions ####
simple_range_extracter <- function(p, scale) {
  d <- ggplot_build(p)
  d$plot$scales$get_scales(scale$aesthetics)$range$range
}
get_shared_scale <- function(..., scale) {
  plots <- list(...)
  ranges <- purrr::map(plots, ~simple_range_extracter(., scale))
  single_range <- range(unlist(ranges))
  scale$limits <- single_range
  scale
}
# Main function
set_scale_union <- function(..., scale) {
  exprs <- rlang::enexprs(...)
  scale <- get_shared_scale(..., scale = scale)
  var_nms <- purrr::map_chr(exprs, rlang::as_name)
  edit_plots_in_place(var_nms, env = parent.frame(),
                      scale = scale)
  # Invisibly return the scale, in case you need it later
  invisible(scale)
}

# Sub-function
edit_plots_in_place <- function(names, env, scale) {
  vars <- rlang::env_has(env = env, nms = names)
  if (!all(vars))
    stop("Environment does not have variables for ",
         paste(names(vars[!vars]), collapse=", "))
  purrr:::walk(names, function(nm) {
    og_plot <- rlang::env_get(env, nm = nm)
    message("Changing plot `", nm, "`")
    # Muffles messages about already having scales
    withCallingHandlers(
      assign(x = nm, envir = env,
             value = og_plot + scale),
      message = function(err) {
        if (grepl("already present", err$message))
          invokeRestart("muffleMessage")
      })
  })
}
# 1.1 Descriptives ####
# age
mean(na.omit(demo10$Age_BL_OFAMS)) # BL
sd(na.omit(demo10$Age_BL_OFAMS))
mean(na.omit(demo10$Age_OFAMS10)) # 10 year follow-up
sd(na.omit(demo10$Age_OFAMS10))
# sex
table(demo10$sex)

# correlation matrix of proms
x = cor(data%>%dplyr::select(contains("spm")),use = "complete.obs")
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = x, col = col, symm = TRUE)
print("The heatmap underscores the findings by Hobart et al. (2001) https://doi.org/10.1136/jnnp.71.3.363 that the subscale sum scores are not that straight forward.")
print("We will anyways go ahead with summing into 8 categories.")
# add proms sum scores (for interpretability divided by nb of item per category)
# first, the proms have to be transformed into percentage scores
trans1 = function(spm){
  a = ifelse(spm == 1,100,spm)
  a = ifelse(a == 2,75,a)
  a = ifelse(a == 3,50,a)
  a = ifelse(a == 4,25,a)
  a = ifelse(a == 5,0,a)
  return(a)
}
translist1 = c("spm1","spm2","spm6","spm8","spm11b","spm11d")
for(i in 1:length(translist1)){
  data[grep(paste("^",translist1[i],"$",sep=""), colnames(data))] = c(as.integer(trans1(unlist(data[grep(paste("^",translist1[1],"$",sep=""), colnames(data))]))))
}
trans2 = function(spm){
  a = ifelse(spm == 1,0,spm)
  a = ifelse(a == 2,50,a)
  a = ifelse(a == 3,100,a)
  return(a)
}
translist2 = data %>% dplyr::select(contains("spm3")) %>% names()
for(i in 1:length(translist2)){
  data[grep(paste("^",translist2[i],"$",sep=""), colnames(data))] = c(as.integer(trans2(unlist(data[grep(paste("^",translist2[1],"$",sep=""), colnames(data))]))))
}
trans3 = function(spm){
  a = ifelse(spm == 1,100,spm)
  a = ifelse(a == 2,80,a)
  a = ifelse(a == 3,60,a)
  a = ifelse(a == 4,40,a)
  a = ifelse(a == 5,20,a)
  a = ifelse(a == 6,0,a)
  return(a)
}
translist3 = c("spm4b", "spm4d","spm9d", "spm9e", "spm9h")
for(i in 1:length(translist3)){
  data[grep(paste("^",translist3[i],"$",sep=""), colnames(data))] = c(as.integer(trans3(unlist(data[grep(paste("^",translist3[1],"$",sep=""), colnames(data))]))))
}
trans4 = function(spm){
  a = ifelse(spm == 1,0,spm)
  a = ifelse(a == 2,20,a)
  a = ifelse(a == 3,40,a)
  a = ifelse(a == 4,60,a)
  a = ifelse(a == 5,80,a)
  a = ifelse(a == 6,100,a)
  return(a)
}
translist3 = c("spm9b", "spm9c", "spm9f", "spm9g", "spm9i")
for(i in 1:length(translist3)){
  data[grep(paste("^",translist3[i],"$",sep=""), colnames(data))] = c(as.integer(trans4(unlist(data[grep(paste("^",translist3[1],"$",sep=""), colnames(data))]))))
}
trans5 = function(spm){
  a = ifelse(spm == 1,0,spm)
  a = ifelse(a == 2,25,a)
  a = ifelse(a == 3,50,a)
  a = ifelse(a == 4,75,a)
  a = ifelse(a == 5,100,a)
  return(a)
}
translist3 = c("spm10","spm11a","spm11c")
for(i in 1:length(translist3)){
  data[grep(paste("^",translist3[i],"$",sep=""), colnames(data))] = c(as.integer(trans5(unlist(data[grep(paste("^",translist3[1],"$",sep=""), colnames(data))]))))
}
trans2 = function(spm){ # role physical and role emotional
  a = ifelse(spm == 0,100,spm)
  a = ifelse(a == 1,0,a)
  return(a)
}
translist2 = data %>% dplyr::select(contains("spm4"), contains("spm5")) %>% names()
for(i in 1:length(translist2)){
  data[grep(paste("^",translist2[i],"$",sep=""), colnames(data))] = c(as.integer(trans2(unlist(data[grep(paste("^",translist2[1],"$",sep=""), colnames(data))]))))
}
# note that the scale levels differ and absolute sum scores have to be used.
data$PF = data%>%dplyr::select(contains("spm3")) %>% rowSums()/10 # Physical functioning (PF)
data$RF = data%>%dplyr::select(contains("spm4")) %>% rowSums()/4 # Role-physical (RF)
data$BP = data%>%dplyr::select(spm7,spm8) %>% rowSums()/2 # Bodily pain (BP)
data$GH = data%>%dplyr::select(spm1,contains("spm11")) %>% rowSums()/4 # General health (GH)
data$VT = data%>%dplyr::select(spm9a,spm9b,spm9c,spm9d) %>% rowSums()/5 # Vitality (VT)
data$SF = data%>%dplyr::select(spm6,spm10) %>% rowSums()/2 # Social functioning (SF)
data$RE = data%>%dplyr::select(contains("spm5")) %>% rowSums()/3 # Role-emotional (RE)
data$MH = data%>%dplyr::select(spm9b, spm9c, spm9d, spm9f, spm3h) %>% rowSums()/5 # Mental health (MH)
#
# correlation matrices of variables of interest
visits = c("0","6","12","24") # at 120 months follow-up, there are no clinical measures taken.
df = unique(data)
for (i in 1:length(visits)){
  x = cor(df%>% filter(session == visits[i])%>%dplyr::select(age,BAG,BAG_c,edss,PF,RF,BP,GH,VT,SF,RE,MH),use = "complete.obs")
  pdf(paste(savepath,"heatmap_",visits[i],"_months",".pdf",sep=""), width=8, height=8)
  heatmap.2(Rowv = F, x = x, col = col, symm = T, cellnote=round(x,2),notecol="black", dendrogram = "none",key=T,density.info="density",trace = "none",notecex = 1.25)
  dev.off()
}
# since there are only very few brain age values at 12 months, we make a cor matrix without it
x = cor(df%>% filter(session == visits[i])%>%dplyr::select(edss,PF,RF,BP,GH,VT,SF,RE,MH),use = "complete.obs")
pdf(paste(savepath,"heatmap_",visits[i],"_months",".pdf",sep=""), width=8, height=8)
heatmap.2(Rowv = F, x = x, col = col, symm = T, cellnote=round(x,2),notecol="black", dendrogram = "none",key=T,density.info="density",trace = "none",notecex = 1.25)
dev.off()
#
print("There are multiple duplicates in the dataset")
print(paste("To be exact: ",nrow(data)-nrow(unique(data))," rows.",sep=""))
print("Another problem are the spread of participants in different session bins:")
data %>% group_by(session) %>% summarize(N=length(na.omit((BAG))))
print("We combat these problems in the representation of initial descriptives by excluding duplicates and low-N sessions.")
df2 = filter(df,session == 0 | session == 6 | session == 24 |session == 120)
print("All time points.")
df$session
ggplot(df, aes(session, BAG_c, group = eid, color = ""))+ #geom_point() + 
  geom_line() + scale_colour_manual(values = c("brown")) +
  theme_classic() + theme(legend.position = "none")
 print("Only time points with sufficient data.")
ggplot(df2, aes(session, BAG_c, group = eid, color = ""))+ #geom_point() + 
   geom_line() + scale_colour_manual(values = c("brown")) +
   theme_classic() + theme(legend.position = "none")
print("All time points.")
ggplot(data= df, aes(x = session, y = BAG_c)) + 
  geom_point(size = 2, alpha= 1, color = "gray85") + #,  aes(color = eid) colour points by group
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="blue",fill = "red")+ #method = "lm"
  ylab("Brain Age Gap") + xlab("Session") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
print("Only time points with sufficient data.")

BAG_descriptive = ggplot(data= df2, aes(x = session, y = BAG_c)) + 
  #scale_shape_manual(values=1:nlevels(factor(df2$eid))) +
  geom_point(size = 1.5, alpha= 1, color = "gray85") + #,  aes(color = eid) colour points by group # shape by eid aes(shape=factor(eid)),
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Brain Age Gap (years)") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
ggsave(paste(savepath,"BAG_descriptive.pdf",sep=""),BAG_descriptive)
df%>%group_by(session)%>%summarize(M=mean(na.omit(edss)))
df3 = filter(df,session == 0 | session == 6 | session == 12 |session == 18 |session == 24)
level_order = c(0,6,12,18,24)
# EDSS_descriptive = ggplot(data= df3, aes(x = factor(session,level=level_order), y = edss)) + 
#   geom_point(size = 1.5, alpha= 1, color = "gray85") +
#   geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
#   stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
#   ylab("EDSS") + xlab("Session (months)") +
#   stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
#   theme_bw()
# #setwd(paste(savepath,sep=""))
df4 = filter(df,session == 0 | session == 6 | session == 12 |session == 24)
PF_descriptive = ggplot(data= df4, aes(x = session, y = PF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Physical functioning") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
RF_descriptive = ggplot(data= df4, aes(x = session, y = RF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Role-physical") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
BP_descriptive = ggplot(data= df4, aes(x = session, y = BP)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Bodily pain") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
GH_descriptive = ggplot(data= df4, aes(x = session, y = GH)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("General health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
Physical_Health = ggarrange(PF_descriptive,RF_descriptive,BP_descriptive,GH_descriptive)
ggsave(paste(savepath,"Physical_Health.pdf",sep=""),Physical_Health)
VT_descriptive = ggplot(data= df4, aes(x = session, y = VT)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("General health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
SF_descriptive = ggplot(data= df4, aes(x = session, y = SF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Social functioning") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
RE_descriptive = ggplot(data= df4, aes(x = session, y = RE)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Role-emotional") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
MH_descriptive = ggplot(data= df4, aes(x = session, y = MH)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Mental health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
Mental_Health = ggarrange(VT_descriptive,SF_descriptive,RE_descriptive,MH_descriptive)
ggsave(paste(savepath,"Mental_Health.pdf",sep=""),Mental_Health)
#
# check atrophy pattern by comparing the first (baseline) to final (10-years) timepoints
# tf0 = anat %>% filter(session == "0") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
#   dplyr::select(-eid, -age, -brainage, -corrected_brainage, -X, -session)
# tf120 = anat %>% filter(session == "120") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
#   dplyr::select(-eid, -age, -brainage, -corrected_brainage, -X, -session)
#
#
#
#
# 1.2 Do any variables change significantly over time? ####
# Variables interesting for final analyses:
# Brain atrophy
# Brain age
# EDSS
# Cumulative relapses
# BL_Age
# BMI, BMI_change
# Omega3_suppl
# Treatment
# BL_BMI, rate of BMI change
# Vitamins
# lesion count & vol
# some data processing:
df = unique(data)
df = df %>% dplyr::select(!starts_with("spm"))
df1 = df %>% filter(session==0)
df2 = df %>% filter(session==120)
df1 = df1[!duplicated((df1)$eid),]
df2 = df2[!duplicated((df2)$eid),]
df1 = df1[df1$eid %in% df2$eid,]
df2 = df2[df2$eid %in% df1$eid,]
df2.1 = rbind(df1%>%dplyr::select(eid,session,BAG_c),df2%>%dplyr::select(eid,session,BAG_c))
df2.1$session = ifelse(df2.1$session == 0,0,120)
df2.1 %>% group_by(session) %>% summarise(M = mean(na.omit(BAG_c)))
df1$BAG_c_diff = (df2$BAG_c-df1$BAG_c)
BAGdf = df1 %>% dplyr::select(eid, BAG_c, BAG_c_diff)
#
demo10$edss_baseline = as.numeric(gsub(",",".",demo10$edss_baseline))
demo10$edss_month_12 = as.numeric(gsub(",",".",demo10$edss_month_12))
demo10$edss_month18 = as.numeric(gsub(",",".",demo10$edss_month18))
demo10$edss_month_24 = as.numeric(gsub(",",".",demo10$edss_month_24))
demo10$EDSSdiff_FUvsM24 = as.numeric(gsub(",",".",demo10$EDSSdiff_FUvsM24))
demo10$edss_month_120 = demo10$EDSSdiff_FUvsM24+demo10$edss_month_24
t.test(demo10$edss_month_120,demo10$edss_baseline,paired=T) # EDSS is higher at 120 months than at baseline
cohens_d(rbind(data.frame(session = "120", edss = demo10$edss_month_120),data.frame(session = "0", edss = demo10$edss_baseline)),edss~session,paired=T,ci=T,ref.group="120")
# showcasing the increase:
mean(na.omit(demo10$edss_month_120))
mean(na.omit(demo10$edss_baseline))
demo10$edss_ARC = (demo10$edss_month_120-demo10$edss_baseline)/(demo10$Age_OFAMS10-demo10$Age_BL_OFAMS)
# 
# Relapse
mean(na.omit(demo10$Cum_relapses))
t.test(na.omit(as.numeric(gsub(",",".",demo10$Relapserate_OFAMS))), mu = 0, alternative = "greater") # relapse rate > 0
t.test(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))), mu = 0, alternative = "greater") # total relapses > 0
min(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
max(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
median(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
mad(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
hist((na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses)))),breaks = 20)
#
# Body mass index (BMI)
t.test(as.numeric(gsub(",",".",demo10$BL_BMI)), as.numeric(gsub(",",".",demo10$BMI_OFAMS10)), paired = T) # no sig differences in BMI
#
#
#
#
#
#
# EDSS plotting #######
demo10$eid=demo10$Patnr
edss_df0 = demo10%>% dplyr::select(eid,Age_OFAMS10,edss_month_120,sex)
edss_df0$session = "144"
names(edss_df0) = c("eid", "age", "edss", "sex", "session")
edss_df = df3 %>% dplyr::select(eid, age, edss, sex, session)
edss_df = rbind(edss_df,edss_df0)
# IMPORTANT: FILTER EMPTY DATA!!
edss_df = edss_df %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
edss_df = edss_df[order(edss_df$eid,edss_df$session),]
#write.csv(edss_df,paste(savepath,"tmp_table.csv",sep=""))
edss_df = read.csv(paste(savepath,"edss_age_table.csv",sep=""))
#
#
#
#
# cross-check edss df
edss.0 = melt(demo10 %>% dplyr::select(eid,edss_baseline,edss_month_12,edss_month18,edss_month_24,edss_month_120),id.vars = "eid")
levels(edss.0$variable)=c(0,12,18,24,144)
names(edss.0) = c("eid","session","edss")
edss_df$session = factor(edss_df$session)
edss.0 = full_join(edss_df,edss.0,by=c("eid","session"))
edss.0$edss = ifelse(is.na(edss.0$edss.x)==T,edss.0$edss.y,edss.0$edss.x)
#
#
idlist = demo10%>%filter(Inclusion==1) %>% dplyr::select(eid)
edss_df = edss.0[edss.0$eid %in% idlist$eid,]
# edss_df=edss_df[order(edss_df$eid,as.numeri(edss_df$ession), edss_df$age),]
# df$age<-na.locf(df$age)
#edss_df %>% group_by(eid,session) %>% arrange(eid,session,age) %>% fill(age)
# 
edss_df$session = factor(edss_df$session)
# to estimate the age and sex adjusted EDSS score, we use an LME
## these are the participant with an edss equal or bigger 4
idlist=edss_df[edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,]%>% dplyr::select(eid) %>% unique  %>% unlist() %>% c()
#211 and 604 improve; 215 and 806 are stable
idlist = idlist[!idlist %in% c(211,215,604,806)]
ids=length(idlist)
# Specify characteristics of the functional decline group
edss_df$FLG = ifelse(edss_df$eid %in% unlist(idlist), 1, 0)
mml = lmer(edss ~ session + age + sex + (1|eid),edss_df[edss_df$eid %in% idlist,])
# get the time point 144 model coefficient
cofff=summary(mml)$coefficients[5]
# plot it
edss_df$session = factor(edss_df$session, levels=c('0', '6', '12', '18','24', '144'))
e1 = edss_df %>% ggplot(aes(x = (session), y = edss)) + 
  geom_line(data = edss_df %>% dplyr::filter(FLG == 0), aes(group = eid),color = "lightgrey") + #gghighlight(use_group_by = T ,use_direct_label = F) +
  geom_point(data = edss_df %>% dplyr::filter(FLG == 0), mapping = aes(x = (session), y = edss), size = 1.5, alpha= 1, color = "lightgrey") + 
  geom_line(data = edss_df %>% dplyr::filter(FLG == 1), aes(group = eid),color = "grey46") + #gghighlight(use_group_by = T ,use_direct_label = F) +
  geom_point(data = edss_df %>% dplyr::filter(FLG == 1), mapping = aes(x = (session), y = edss), size = 1.5, alpha= 1, color = "grey46") + 
  stat_smooth(data = edss_df %>% dplyr::filter(FLG == 1), aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("EDSS") + xlab("Session (months)") +
  stat_summary(data = edss_df %>% dplyr::filter(FLG == 1),aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=3,y=8,label= paste("Adjusted 12-year difference =",round(cofff,2),"N =",ids),cex=4) + 
  ggtitle("Functional loss group")
mml2 = lmer(edss ~ session + age + sex + (1|eid),edss_df[!edss_df$eid %in% idlist,])
cofff=summary(mml2)$coefficients[5]
#ids=edss_df[!edss_df$eid %in% (edss_df %>% filter(session ==144 & edss<4))$eid,]%>% dplyr::select(eid) %>% unique  %>% unlist() %>% length
ids = length(unique(edss_df[!edss_df$eid %in% idlist,]$eid))
#ids=(edss_df[!edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,])$eid%>%unique() %>% unlist() %>% length
level_order = c(0,6,12,18,24,144)
e2 = edss_df %>% ggplot(aes(x = (session), y = edss)) + 
  geom_line(data = edss_df %>% dplyr::filter(FLG == 1), aes(group = eid),color = "lightgrey") + #gghighlight(use_group_by = T ,use_direct_label = F) +
  geom_point(data = edss_df %>% dplyr::filter(FLG == 1), mapping = aes(x = (session), y = edss), size = 1.5, alpha= 1, color = "lightgrey") + 
  geom_line(data = edss_df %>% dplyr::filter(FLG == 0), aes(group = eid),color = "grey46") + #gghighlight(use_group_by = T ,use_direct_label = F) +
  geom_point(data = edss_df %>% dplyr::filter(FLG == 0), mapping = aes(x = (session), y = edss), size = 1.5, alpha= 1, color = "grey46") + 
  stat_smooth(data = edss_df %>% dplyr::filter(FLG == 0), aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("EDSS") + xlab("Session (months)") +
  stat_summary(data = edss_df %>% dplyr::filter(FLG == 0),aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=3,y=8,label= paste("Adjusted 12-year difference =",round(cofff,2),"N =",ids),cex=4) + 
  ggtitle("Functionally stable and improving group")
EDSS_descriptive = ggarrange(e1,e2)
# col = ifelse(edss_df$edss>=4,"red","gray85")
# EDSS_descriptive = ggplot(data= edss_df, aes(x = factor(session,level=level_order), y = edss)) + 
#   geom_point(size = 1.5, alpha= 1, color = "gray85") +
#   geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
#   stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
#   ylab("EDSS") + xlab("Session (months)") +
#   stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
#   theme_bw()
ggsave(paste(savepath,"EDSS_descriptive.pdf",sep=""),EDSS_descriptive, device = cairo_pdf,width = 9.5, height = 4)
#
#
#
# check whether the FLG edss cases can be counted as confirmed disability progression
edss = dcast(edss_df %>% dplyr::select(eid,session, edss), eid ~session)
edss[edss$eid %in% idlist,]
#
# Note: We cannot confirm the disability progression for many of the subjects in the FLG.
# This is partly due to the fact that many of the FLG members reach an EDSS ≥ 4 at the 10-year follow-up.
#
#
#
# Another approach to test the robustness of the edss scores are the number of attacks in close proximity to the measurements
# That way, it might be confirmed whether a peak in edss score might be due to relapse activity
background$eid = as.numeric(background$patno)
background = background %>% dplyr::select(eid,SC_DATEOFVISIT)
#background[as.numeric(background$eid) %in% idlist,]
demo10= merge(demo10,background,by="eid")
# here we can check whether the start date overlaps with 
# date of first visit: SC_DATEOFVISIT
# date of last OFAMS visit: Date_last_OFAMS
# date of follow-up: Date_study_visit
demo10[demo10$eid %in% idlist,] %>% dplyr::select(eid,SC_DATEOFVISIT,Date_last_OFAMS,Date_study_visit,Date_relapse1,Date_relapse2,Date_relapse3,Date_relapse4,Date_relapse5,Date_relapse6,Date_relapse7,Date_relapse8,Date_relapse9)
#
# there are only 3 participants in the functional loss group (FLG) who had attacks during the study period
edss %>% filter(eid %in% c(1503,809,813))
#edss %>% filter(eid %in% idlist)
#
#
#
#
# Paced auditory serial addition test (PASAT) assessment ####
t.test(demo10$PASAT_OFAMS10,demo10$BL_PASATcorrect,paired = T) # PASAT score decreases
# However, the simple t-test does not consider the fluctuations in PASAT scores over time - we need LMMs
pasat1 = demo10%>%dplyr::select(Patnr,sex,BL_PASATcorrect,PASAT_24M, PASAT_OFAMS10)
pasat1 = melt(pasat1,id.vars = c("Patnr","sex"))
names(pasat1) = c("eid","sex","session","PASAT")
pasat1$session=as.numeric(pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 1,0,pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 2,24,pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 3,144,pasat1$session)
pasat1$session=factor(pasat1$session)
pasat1 = pasat1 %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
pasat1 = edss_df %>% right_join(pasat1, by = c("eid","session","sex"))
#
#
pasat1$FLG = ifelse(pasat1$eid %in% unlist(idlist), 1, 0)
# 
#
#
# timed 25-foot walk test (T25FW) at 24 and 120 months
t.test(as.numeric(gsub(",",".",demo10$T_25FW_2)), as.numeric(gsub(",",".",demo10$T_25FW_1)), paired = T) # sig differences in T25FW
cohens_d(rbind(data.frame(session="24",T_25FW=as.numeric(gsub(",",".",demo10$T_25FW_1))),
               data.frame(session="120",T_25FW=as.numeric(gsub(",",".",demo10$T_25FW_2)))),
         T_25FW~session, ref.group="120", paired=T, ci=T)
demo10$T_25FW_ARC = (as.numeric(gsub(",",".",demo10$T_25FW_2))-as.numeric(gsub(",",".",demo10$T_25FW_1)))/14
#
# check proms (show no difference between BL and FU) [BL to 24 months]
df1 = df %>% filter(session==0)
df2 = df %>% filter(session==24)
df1 = df1[!duplicated((df1)$eid),]
df2 = df2[!duplicated((df2)$eid),]
df1 = df1[df1$eid %in% df2$eid,]
df2 = df2[df2$eid %in% df1$eid,]


#
# test total volume degeneration
icv %>% filter(session == 0 | session == 120) %>% t_test(TotalVol~session)
#icv %>% filter(session == 0 | session == 120) %>% cohens_d(TotalVol~session,ref.group="120", paired=T, ci=T)
t.test((icv %>% filter(session == 0))$TotalVol,(icv %>% filter(session == 120))$TotalVol)
#
icv$session = ifelse(icv$session == 120, 144, icv$session)
icv$session = factor(icv$session)
#vol.m = lmer(TotalVol~session+TIV+(1|eid),icv %>% filter(session == 0 | session == 120))
#summary(vol.m) # estimate also the TIV corrected model (also reported in the article)
ndf = full_join(icv,edss_df, by = c("eid","session"))
ndf = full_join(ndf,geno,by="eid")
ndf = ndf %>% dplyr::select(-X)
anat$session=ifelse(anat$session==120,144,anat$session)
anat$session = factor(anat$session)
ndf = left_join(ndf,anat %>% dplyr::select(-age), by = c("eid","session"),relationship = "many-to-many")
#df%>%filter(eid==201 & session == 0) %>% dplyr::select(TotalVol)
blood$session = factor(blood$Visit.nr)
blood$eid = (blood$Sub)
blood = blood %>% dplyr::select(eid,session,NfL..pg.ml.,CH3L.1..mg.ml..mean)

ndf = full_join(ndf,blood, by = c("eid","session"),relationship = "many-to-many")

write.csv(ndf,paste(savepath,"eTIV_TotalVol.csv",sep=""))


# fatigue 
fati$fatigue = rowSums(fati[9:ncol(fati)])/9 # calculate fatigue sum score (higher = more)
fati$eid = fati$patno
fati1 = fati %>% filter(VISIT == "Month 24") %>% dplyr::select(eid, fatigue)
fati2 = fati %>% filter(VISIT == "Baseline") %>% dplyr::select(eid, fatigue)
fati1 = right_join(fati1,fati2,by="eid")
fati1$fatigue = fati1$fatigue.y-fati1$fatigue.x
names(fati1) = c("eid","fatigue_M24", "fatigue_BL", "fatigue_diff24")

prom = right_join(prom,data %>% dplyr::select(eid,session,PF,RF,BP,GH,VT,SF,RE,MH),by = c("eid","session"))
demo10$eid = demo10$Patnr
prom = prom %>% filter(session == 0) # dplyr::select only baseline measure (there is anyways no change)
prom = right_join(demo10,prom, by="eid")
prom$BL_BMI = as.numeric(gsub(",",".",prom$BL_BMI))
prom$BMI_OFAMS10 = as.numeric(gsub(",",".",prom$BMI_OFAMS10))
prom$FatigueSS_score = as.numeric(gsub(",",".",prom$FatigueSS_score))
prom$Vit_A_0 = as.numeric(gsub(",",".",prom$Vit_A_0))
prom$Vit_D_0 = as.numeric(gsub(",",".",prom$Vit_D_0))
prom$Vit_E_0 = as.numeric(gsub(",",".",prom$Vit_E_0))
#
#remove duplicates
#full_join(icv %>% filter(session==0),icv %>% filter(session==120),by=c("eid","session"))
icv2 = icv %>% filter(session == "0" | session == "120") %>% unique()
remove_c = c()
for (i in 1:nrow(icv2)){ # loop provides Boolean-ish for exclusion based on duplicates (seciond is being kept)
  remove_c[i]=(ifelse(icv2$eid[i] == icv2$eid[i-1],F,T)) # T/F = keep/remove
}
remove_c = ifelse(is.na(remove_c) == T, T,remove_c) # one row in beginning must be NA due to loop implementation
icv = icv2[remove_c,]
icv2 = icv%>%filter(session==0)
icv3 = icv%>%filter(session==120)
icv_change = right_join(icv2,icv3,by="eid")
icv_change$TotalVol_change=(icv_change$TotalVol.y-icv_change$TotalVol.x)
icv_change$TIV_change=(icv_change$TIV.y-icv_change$TIV.x)
#
# prom = right_join(vols,prom,by="eid") # proms
# prom$vol_ARC = prom$vol_diff/(prom$Age_OFAMS10-prom$Age_BL_OFAMS)
prom = right_join(BAGdf,prom,by="eid")
prom=right_join(prom,geno,by="eid") # merge relevant dfs
prom$Mean_BMI_OFAMS = as.numeric(gsub(",",".",prom$Mean_BMI_OFAMS)) # add bmi
prom$age_gap = prom$Age_OFAMS10-prom$Age_BL_OFAMS # add age gap
prom = right_join(icv2,prom,by="eid")
# estimate cohen's d for changes
l = right_join(lesion_count %>% dplyr::select(eid,baseline),lesion_vol %>% dplyr::select(eid,baseline),by="eid")
names(l) = c("eid","baselineC","baselineV")
prom = right_join(l,prom,by="eid")
BL_blood = blood %>% filter(session==0)
prom=right_join(prom,BL_blood,by="eid")
prom=prom%>%dplyr::select(-Sex_OFAMS10,-sex)
jj=full_join(edss_df %>% filter(session==0),pasat1 %>% filter(session==0) %>% dplyr::select(eid,PASAT), by="eid")
prom=full_join(jj,prom,by="eid")
prom = prom %>% dplyr::select(eid,FLG,age,sex,edss,PASAT,geno,relapses_12mnths_before_baseline,
                       CH3L.1..mg.ml..mean,NfL..pg.ml.,smoking_OFAMS,BL_BMI,
                       Treatment_OFAMS,Omega3_suppl,baselineC,baselineV,
                       PF,RF,BP,GH,VT,SF,RE,MH,Vit_A_0,Vit_D_0,Vit_E_0)
prom=full_join((df %>% filter(session == 0) %>% dplyr::select(eid,BAG_c)),prom,by=c("eid"))
prom = (unique(prom))
prom = prom[prom$eid %in% unique(edss_df$eid),]
write.csv(prom,paste(savepath,"interrim_data.csv",sep=""))
#
# Paced auditory serial addition test (PASAT) plotting ####
# check the numbers for PASAT decline below a considerably and meaningfully low threshold
pasat1 = pasat1[pasat1$eid %in% prom$eid,]

levels(pasat1$session)
summary(lm(PASAT~age+sex,data = pasat1%>%filter(session==0)))
summary(lm(PASAT~age+sex,data = pasat1%>%filter(session==24)))
summary(lm(PASAT~age+sex,data = pasat1%>%filter(session==144)))

hist((pasat1%>%filter(session==0))$PASAT,breaks = 100)
hist((pasat1%>%filter(session==24))$PASAT,breaks = 100)
hist((pasat1%>%filter(session==144))$PASAT,breaks = 100)
pasat_ids=(pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40| session ==18 & PASAT<=40| session ==12 & PASAT<=40| session ==6 & PASAT<=40))$eid,])$eid%>%unique() %>% unlist()
pasat_ids%>% length
#pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT>=4 | session ==24 & PASAT>=4| session ==18 & PASAT>=4| session ==12 & PASAT>=4| session ==6 & PASAT>4))$eid,]%>% dplyr::select(eid) %>% unique  %>% unlist() %>% c()
#pasat1 %>% filter(session == 144 & PASAT<=40) %>% summarize(mean(PASAT)) -pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40))$eid,]%>% filter(session ==0) %>% summarize(mean(PASAT))
pasat1$FSL = pasat1$eid %in% pasat_ids
mml = lmer(PASAT ~ session + age + sex + (1|eid),pasat1[pasat1$eid %in% pasat_ids,])
cofff=summary(mml)$coefficients[3] # get corrected coefficient / difference

# get ids 
pasat_ids=(pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=30 | session ==24 & PASAT<=30| session ==18 & PASAT<=30| session ==12 & PASAT<=30| session ==6 & PASAT>4))$eid,])$eid%>%unique() %>% unlist() %>% length
e1 = ggplot(data= pasat1, aes(x = factor(session,level=level_order), y = PASAT)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(min(PASAT)<=30,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("PASAT") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=2,y=65,label= paste("Adjusted 10-year difference=",round(cofff,2),"N =",pasat_ids),cex=4) + 
  ggtitle("Reached PASAT≤40")
mml = lmer(PASAT ~ session + age + sex + (1|eid),pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT>40 | session ==24 & PASAT>40| session ==0 & PASAT>40))$eid,])
cofff=summary(mml)$coefficients[3]

pasat_ids=(pasat1[!pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=30 | session ==24 & PASAT<=30| session ==18 & PASAT<=30| session ==12 & PASAT<=30| session ==6 & PASAT<=30))$eid,])$eid%>%unique() %>% unlist() %>% length

e2 = ggplot(data= pasat1, aes(x = factor(session,level=level_order), y = PASAT)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(min(PASAT)>30,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("PASAT") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=2,y=65,label= paste("Adjusted 10-year difference=",round(cofff,2),"N =",pasat_ids),cex=4) + 
  ggtitle("Reached PASAT≤40")
ggarrange(e1,e2)
# check N in stratified group
pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40))$eid,] %>% group_by(session) %>% summarize(N = length(session))

#
pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40))$eid,] %>% dplyr::select(eid) %>% unique


# Relationship between PASAT and EDSS ####
# Interestingly, PASAT and EDSS appear unrelated in our sample.
summary(lmer(PASAT~edss+age+sex+(1|eid),pasat1))
# This underscores the importance of testing PASAT (cognitive test) separately.


# Relationship between EDSS and lesions ####
#
#
# in addition to baseline and m120, we take the sum of the new lesions within the nearest 6 months
lesion_count$baseline = ifelse(is.na(lesion_count$baseline),0,lesion_count$baseline)
lesion_count$m6s = ifelse(is.na(lesion_count[order(lesion_count$eid),] %>%
  dplyr::select(m1,m2,m3,m4,m5,m6) %>%rowMeans(na.rm = T)),0,lesion_count[order(lesion_count$eid),] %>%
    dplyr::select(m1,m2,m3,m4,m5,m6) %>%rowMeans(na.rm = T))
lesion_count$m12s = ifelse(is.na(lesion_count[order(lesion_count$eid),] %>%
  dplyr::select(m7,m8,m9,m12) %>%rowMeans(na.rm = T)),0,lesion_count[order(lesion_count$eid),] %>%
    dplyr::select(m7,m8,m9,m12) %>%rowMeans(na.rm = T))
lesion_count$m18s = ifelse(is.na(lesion_count[order(lesion_count$eid),] %>%
  dplyr::select(m15,m18) %>%rowMeans(na.rm = T)),0,lesion_count[order(lesion_count$eid),] %>%
    dplyr::select(m15,m18) %>%rowMeans(na.rm = T))
lesion_count$m120 = ifelse(is.na(lesion_count$m120),0,lesion_count$m120)
# merge edss and lesion data frames for time-point by time-point associations
names(edss)=c("eid","edss0","edss6","edss12","edss18","edss24","edss144")
lc = merge(lesion_count,edss,by="eid")
summary(lm(edss0~baseline,lc))
summary(lm(edss6~(m6s-baseline),lc))
summary(lm(edss12~(m12s-m6s),lc))
summary(lm(edss18~m18s,lc))
summary(lm(edss144~m120,lc))
#211 and 604 improve; 215 and 806 are stable
#edss %>% filter(eid %in% idlist) %>%filter(!eid %in% c(211,215,604,806)) %>% nrow()
#
# on can do lmers as well, but this is not relevant here
lesion_count = lesion_count[1:16]
lc = melt(lesion_count,id.vars = "eid")
levels(lc$variable) = c(0,1,2,3,4,5,6,7,8,9,12,24,144,15,18)
e = melt(edss,id.vars = "eid")
levels(e$variable) = c(0,6,12,18,24,144)
lc = full_join(lc,e,by=c("eid","variable"))
lc$FLG = ifelse(lc$eid %in% unlist(idlist), 1, 0) # add flg for group comp
# show the assiciation between edss and lesion count
summary(lmer(value.y~value.x+(1|eid),lc)) # y = edss, x = count