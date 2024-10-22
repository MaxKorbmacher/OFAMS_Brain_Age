# merge data with to train brain age models for predictions in MS
# Max Korbmacher 11 Oct 2024
#
#
#
# prep (load data and pkgs) ####
setwd("/cluster/projects/p33/users/maxk/Lifespan/data")
#
# lib
library(dplyr)
library(MatchIt)
library(readxl)
# load labels
ABCD_labels = read.csv("ABCD_labels.csv")
HBN_labels = read.csv("HBN_labels.csv")
HCP_labels = read.csv("HCP_labels.csv")
UKB_labels = read.csv("UKB_labels.csv")
TOP_labels = read.csv("TOP_labels.csv")
Rockland_labels = read.csv("RocklandSample_labels.csv")
demgen_labels = read.csv("demgen_labels.csv")
AddNeuroMed_labels = read.csv("AddNeuroMed_labels.csv")
MS_labels = read.csv("ms_labels.csv")
ADNI_young = read.csv("ADNI-young_labels.csv")
ADNI_old = read.csv("ADNI-old_labels.csv")
ADNI_orig = read.csv("ADNI_orig.csv")
BBSC_labels = read.csv("BBSC_labels.csv")
#FTHP_labels = read_excel("/cluster/projects/p33/groups/imaging/FTHP/FTHP_Metadata.xlsx")
hubin_labels = read.csv("HUBIN_KASP_clinical_data/HUBIN_long_clean.csv")
hubin_labels2 = read.csv("HUBIN/Covariates.csv")
kasp_labels = read.csv("HUBIN_KASP_clinical_data/KASP_basic.csv")
#
#
# load data
abcd_thick = read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/imaging/mri_y_smr_thk_dsk.csv")
abcd_area = read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/imaging/mri_y_smr_area_dsk.csv")
abcd_vol = read.csv("/cluster/projects/p33/groups/imaging/abcd/pheno/ABCDStudyNDA_5.1/core/imaging/mri_y_smr_vol_dsk.csv")
HBN = read.csv("HBN_data.csv")
HCP = read.csv("HCP.csv")
UKB = read.csv("UKB.csv")
top_GE3T = read.csv("top_GE3T.csv")
top15 = read.csv("top15.csv")
top3 = read.csv("top3.csv")
# top_hubin_surf = read.csv("HUBIN/CorticalMeasuresENIGMA_SurfAvg.csv")
# top_hubin_vol = read.csv("HUBIN/LandRvolumes.csv")
# top_hubin_thick = read.csv("HUBIN/CorticalMeasuresENIGMA_ThickAvg.csv")
top_hubin = read.csv("top_hubin.csv")
top_kasp = read.csv("top_kasp.csv")
Rockland = read.csv("Rockland.csv")
demgen3T = read.csv("demgen3T.csv")
demgenGE750 = read.csv("demgenGE750.csv")
AddNeuroMed = read.csv("AddNeuroMed.csv")
MS = read.csv("ms.csv")
ADNI = read.csv("ADNI.csv")
BBSC1=read.csv("BBSC1.csv")
BBSC2=read.csv("BBSC2.csv")
BBSC3=read.csv("BBSC3.csv")
#FTHP=read.csv("FTHP.csv") # might be later used for validation
#
#
# standardize and prep for merging ####
# abcd
abcd_thick=abcd_thick %>% filter(eventname == "baseline_year_1_arm_1") %>% dplyr::select(-eventname) # filter for first time point (has most time points)
abcd_area=abcd_area %>% filter(eventname == "baseline_year_1_arm_1") %>% dplyr::select(-eventname)
abcd_vol=abcd_vol %>% filter(eventname == "baseline_year_1_arm_1") %>% dplyr::select(-eventname)
#
names(abcd_area)=gsub("smri_area_cdk_","",names(abcd_area)) # remove useless info from column names
names(abcd_thick)=gsub("smri_thick_cdk_","",names(abcd_thick)) # remove useless info from column names
names(abcd_vol)=gsub("smri_vol_cdk_","",names(abcd_vol)) # remove useless info from column names
hemi_prefix = function(df){ # function to return vector indicating the hemisphere which will become the prefix of all columns
  left=ifelse(grepl("lh",names(df)),"lh_","")
  right=ifelse(grepl("rh",names(df)),"rh_","")
  prefix=ifelse(left=="",right,left)
  return(prefix)
}
area_prefix=hemi_prefix(abcd_area) # vector of the prefixes
thick_prefix=hemi_prefix(abcd_thick)
vol_prefix=hemi_prefix(abcd_vol)
#
correctnames = function(data){ # make a function for correct naming of the brain regions
  names(data)=gsub("bankssts","bankssts_",names(data))
  names(data)=gsub("cdacate","caudalanteriorcingulate_",names(data))
  names(data)=gsub("cdmdfr","caudalmiddlefrontal_",names(data))
  names(data)=gsub("cuneus","cuneus_",names(data))
  names(data)=gsub("ehinal","entorhinal_",names(data))
  names(data)=gsub("fusiform","fusiform_",names(data))
  names(data)=gsub("ifpl","inferiorparietal_",names(data))
  names(data)=gsub("iftm","inferiortemporal_",names(data))
  names(data)=gsub("ihcate","isthmuscingulate_",names(data))
  names(data)=gsub("locc","lateraloccipital_",names(data))
  names(data)=gsub("lobfr","lateralorbitofrontal_",names(data))
  names(data)=gsub("lingual","lingual_",names(data))
  names(data)=gsub("mobfr","medialorbitofrontal_",names(data))
  names(data)=gsub("mdtm","middletemporal_",names(data))
  names(data)=gsub("parahpal","parahippocampal_",names(data))
  names(data)=gsub("paracn","paracentral_",names(data))
  names(data)=gsub("parsopc","parsopercularis_",names(data))
  names(data)=gsub("parsobis","parsorbitalis_",names(data))
  names(data)=gsub("parstgris","parstriangularis_",names(data))
  names(data)=gsub("pericc","pericalcarine_",names(data))
  names(data)=gsub("postcn","postcentral_",names(data))
  names(data)=gsub("ptcate","posteriorcingulate_",names(data))
  names(data)=gsub("precn","precentral_",names(data))
  names(data)=gsub("pc","precuneus_",names(data))
  names(data)=gsub("rracate","rostralanteriorcingulate_",names(data))
  names(data)=gsub("rrmdfr","rostralmiddlefrontal_",names(data))
  names(data)=gsub("sufr","superiorfrontal_",names(data))
  names(data)=gsub("supl","superiorparietal_",names(data))
  names(data)=gsub("sutm","superiortemporal_",names(data))
  names(data)=gsub("sm","supramarginal_",names(data))
  names(data)=gsub("frpole","frontalpole_",names(data))
  names(data)=gsub("tmpole","temporalpole_",names(data))
  names(data)=gsub("trvtm","transversetemporal_",names(data))
  names(data)=gsub("insula","insula_",names(data))
  names(data)=gsub("rh","",names(data))
  names(data)=gsub("lh","",names(data))
  names(data)=gsub("entoinal","entorhinal",names(data))
  return(data)
}
names(abcd_area)=paste(area_prefix,correctnames(abcd_area) %>% names(),"area",sep="") # now, rename columns
names(abcd_thick)=paste(thick_prefix,correctnames(abcd_thick) %>% names(),"thickness",sep="")
names(abcd_vol)=paste(vol_prefix,correctnames(abcd_vol) %>% names(),"volume",sep="")
abcd_area = abcd_area %>% dplyr::select(-contains("total")) # remove averages (left,right,whole brain)
abcd_thick = abcd_thick %>% dplyr::select(-contains("total")) 
abcd_vol = abcd_vol %>% dplyr::select(-contains("total")) 
colnames(abcd_area)[colnames(abcd_area) == 'src_subject_idarea'] <- 'eid' # label id variable correctly
colnames(abcd_thick)[colnames(abcd_thick) == 'src_subject_idthickness'] <- 'eid'
colnames(abcd_vol)[colnames(abcd_vol) == 'src_subject_idvolume'] <- 'eid'
ABCD = merge(abcd_thick,abcd_vol,by="eid") # merge volumentrics, surf, thickness estimates
ABCD = merge(ABCD,abcd_area,by="eid")
ABCD$eid = gsub("NDAR_","",ABCD$eid) # modify eid to match demographics df
ABCD_labels$subject = gsub("sub-NDAR","",ABCD_labels$subject)
colnames(ABCD_labels)[colnames(ABCD_labels) == 'subject'] <- 'eid'
#
# hbn
HBN$eid=substr(HBN$eid,1,12)# keep first twelve chars
HBN_labels$subject=gsub("sub-","",HBN_labels$subject)
colnames(HBN_labels)[colnames(HBN_labels) == 'subject'] <- 'eid'
#
# hpc
colnames(HCP_labels)[colnames(HCP_labels) == 'image_id'] <- 'eid'
#
# ukb
colnames(UKB_labels)[colnames(UKB_labels) == 'image_id'] <- 'eid'
#
# top
top15$eid=gsub("/","-15t",top15$eid)
top3$eid=gsub("/","-3t",top3$eid)
top_GE3T$eid=gsub("/","-3tge750",top_GE3T$eid)
top_kasp$eid=gsub("/","",top_kasp$eid)
top_hubin$eid=gsub("/","",top_hubin$eid)
#
#
#
#
#
#
# NOTE: TOP DATA FOR KASP ARE NOT COMPLETE
# Update is running!
#
#
#
#
#
#
TOP = rbind(top15,top3,top_GE3T,top_kasp,top_hubin)
#
TOP = TOP %>% dplyr::select(!contains("aparc")) %>% select(!contains("BrainSeg")) %>% select(!contains("eTIV"))
colnames(TOP_labels)[colnames(TOP_labels) == 'image_id'] <- 'eid'
kasp_labels$eid = paste("KASP_",kasp_labels$SubjID,sep="")
kasp_labels$subject = kasp_labels$eid
kasp_labels$session = 1
kasp_labels$run = 1
kasp_labels$scanner = "KarolinskaGE15SignaHDxt"
kasp_labels$age = kasp_labels$Age
kasp_labels$sex = kasp_labels$Sex
kasp_labels$diagnosis = ifelse(kasp_labels$Psychotic_disorder_current==1, "PSY", "HC")
kasp_labels$handedness = NA
kasp_labels$bmi = NA
kasp_labels$socioeconomic_status = kasp_labels$SubjectSES_Scale
kasp_labels$fluid_intelligence = NA
kasp_labels$neuroticism = NA
kasp_labels$iq = kasp_labels$IQ
# labels should contain:
# [1] "eid"                  "subject"              "session"              "run"                 
# [5] "scanner"              "age"                  "sex"                  "diagnosis"           
# [9] "handedness"           "bmi"                  "socioeconomic_status" "fluid_intelligence"  
# [13] "neuroticism"          "iq" 
hubin_labels$eid = paste("h",hubin_labels$ScanIDU2,sep="")
hubin_labels$subject = paste("h",hubin_labels$ScanIDU2,sep="")
hubin_labels$session = 1
hubin_labels$run = 1
hubin_labels$scanner = "3tge750"
hubin_labels$age = hubin_labels$Age_at_MRI_1
hubin_labels$sex = hubin_labels$Sex
hubin_labels$diagnosis = ifelse(hubin_labels$Diagnosis_patient_all=="", "HC", hubin_labels$Diagnosis_patient_all)
hubin_labels$handedness = hubin_labels$Handedness
hubin_labels$bmi = NA

hubin_labels2$SubjID = gsub("hz","",hubin_labels2$SubjID)
hubin_labels2$SubjID = gsub("h","",hubin_labels2$SubjID)
hubin_labels2$HUBIN_ID = hubin_labels2$SubjID
#
hubin_labels$HUBIN_ID = gsub("hz","",hubin_labels$HUBIN_ID)
hubin_labels$HUBIN_ID = gsub("h","",hubin_labels$HUBIN_ID)
hubin_labels$HUBIN_ID = gsub("HZ","",hubin_labels$HUBIN_ID)
hubin_labels$HUBIN_ID = gsub("HU","",hubin_labels$HUBIN_ID)
hubin_labels = merge(hubin_labels,hubin_labels2,by="HUBIN_ID")
hubin_labels$socioeconomic_status = NA
hubin_labels$fluid_intelligence = NA
hubin_labels$neuroticism = NA
hubin_labels$iq = hubin_labels$IQ
# merge top labels
TOP_labels = rbind(TOP_labels,
      hubin_labels[,names(hubin_labels) %in% names(TOP_labels)],
      kasp_labels[,names(kasp_labels) %in% names(TOP_labels)])
#
#
# Rockland
Rockland = Rockland[grepl("BAS1",Rockland$lh.aparc.thickness),] # select only first time point
Rockland$eid=substr(Rockland$eid,5,13) # extract id from scanner id string
colnames(Rockland_labels)[colnames(Rockland_labels) == 'subject'] <- 'eid'
#
# demgen
demgen=rbind(demgen3T, demgenGE750)
colnames(demgen_labels)[colnames(demgen_labels) == 'subject'] <- 'eid'
demgen$eid=substr(demgen$eid,10,14) # extract id from scanner id string
#
# AddNeuroMed
AddNeuroMed$image_id = substr(AddNeuroMed$eid,13,48) # extract image id from scanner id string
#
# MS
MS$eid=gsub("/","",MS$eid) # remove unwanted chars from id string
colnames(MS_labels)[colnames(MS_labels) == 'id'] <- 'eid' # rename id var in demographics df
MS_labels = MS_labels %>% select(eid,age,sex,diagnosis,scanner)
MS_labels$session = MS_labels$run = MS_labels$handedness = MS_labels$bmi = MS_labels$socioeconomic_status = MS_labels$fluid_intelligence = MS_labels$neuroticism = MS_labels$iq = ""
#
# adni
ADNI$eid=gsub("/","",ADNI$eid) # remove unwanted chars from id string
colnames(ADNI_orig)[colnames(ADNI_orig) == 'id'] <- 'eid' # rename id var in demographics df
#
#
# merge demographics and respective brain stats ####
ABCD=merge(ABCD,ABCD_labels,by="eid")
HBN=merge(HBN,HBN_labels,by="eid")
HCP=merge(HCP,HCP_labels,by="eid")
UKB = merge(UKB %>% select(-c(age,X)),UKB_labels,by="eid")
TOP=merge(TOP,TOP_labels,by="eid")
Rockland=merge(Rockland,Rockland_labels,by="eid")
demgen = merge(demgen,demgen_labels,by="eid")
AddNeuroMed = merge(AddNeuroMed,AddNeuroMed_labels,by="image_id")
MS=merge(MS,MS_labels,by="eid")
rm(MS_labels,AddNeuroMed_labels,demgen_labels,ABCD_labels,UKB_labels,TOP_labels,HBN_labels,HCP_labels,abcd_area,abcd_thick,abcd_vol,top_GE3T,top_hubin,top_kasp,top15,top3)#remove clutter
#
# merge final data frames with each other
ABCD = ABCD %>% select(-contains("mean")) %>% select(-image_id) %>% select(-contains("WhiteSurf")) 
HBN = HBN %>% select(!contains("aparc")) %>% select(-contains("mean"))%>% select(-contains("WhiteSurf"))  %>% select(!image_id)
HCP = HCP %>% select(-subject)%>% select(!contains("aparc")) %>% 
  select(-contains("mean"))%>% select(-contains("WhiteSurf"))
UKB = UKB %>% select(-contains("mean")) %>% select(-contains("Left"))%>% select(-contains("Right"))%>% select(-subject)
TOP=TOP %>% select(!contains("aparc")) %>% select(!contains("eTIV")) %>% 
  select(!contains("BrainSegVol")) %>% select(-contains("mean")) %>% 
  select(-subject) %>% select(-contains("WhiteSurf"))
UKB = UKB %>% select(!contains("Cortex"))%>% select(!contains("WhiteSurf"))%>% select(!contains("CorticalWhiteMatterVol"))
Rockland=Rockland[,names(Rockland) %in% names(UKB)]
#Rockland=Rockland[,names(UKB)]
demgen=demgen[,names(demgen) %in% names(UKB)]
AddNeuroMed=AddNeuroMed[,names(AddNeuroMed) %in% names(UKB)]
MS=MS[,names(MS) %in% names(UKB)]
#
# adni
ADNI = merge(ADNI,ADNI_orig,by='eid')
ADNI$eid = ADNI$subject
colnames(ADNI_old)[colnames(ADNI_old) == 'subject'] <- 'eid' # rename id var in demographics df
ADNI=ADNI %>% dplyr::select(-c(session,run,scanner,age,sex,diagnosis))
ADNI=merge(ADNI, ADNI_old, by = "eid")
ADNI=ADNI[,names(ADNI) %in% names(UKB)]
#
# BBSC
BBSC = rbind(BBSC1,BBSC2,BBSC3)
BBSC = BBSC %>% select(!contains("aparc")) %>% select(!contains("BrainSeg")) %>% select(!contains("eTIV"))
BBSC$eid=gsub("\\..*","",BBSC$eid)
BBSC_labels$image_id = BBSC_labels$subject = BBSC_labels$session = BBSC_labels$eid = gsub("anat_","",BBSC_labels$eid)
BBSC_labels$run = 1
BBSC_labels$session = gsub("_T1w","",BBSC_labels$image_id)
BBSC_labels$session=gsub("ses-*","",BBSC_labels$session)
BBSC_labels$sex="M"
BBSC_labels$handedness="R"
BBSC_labels$scanner = "3tge750"
BBSC_labels$diagnosis = "HC"
BBSC_labels$bmi = BBSC_labels$socioeconomic_status = BBSC_labels$fluid_intelligence = BBSC_labels$neuroticism = BBSC_labels$iq = NA
BBSC=merge(BBSC_labels,BBSC,by="eid")
#
#
#
# label data
ABCD$data = "ABCD"
HBN$data = "HBN"
HCP$data = "HCP"
UKB$data = "UKB"
TOP$data = "TOP"
Rockland$data = "Rockland"
demgen$data = "demgen"
AddNeuroMed$data = "AddNeuroMed"
MS$data = "MS"
# note: the following datasets are not included when establishing age-matches, as they are longitudinal
ADNI$data = "ADNI"
BBSC$data = "BBSC"
BBSC = BBSC[,names(BBSC) %in% names(UKB)]
ADNI = ADNI[,names(ADNI) %in% names(UKB)]
long_validation = rbind(ADNI, BBSC)
#FTHP$data = "FTHP" # FTHP data might be later used for validation
#
# check for differences in case of interest.
#setdiff(names(TOP),names(HBN))
#setdiff(names(HCP),names(HBN))
# otherwise merge&save the file
final=rbind(ABCD,HBN,HCP,UKB,Rockland,AddNeuroMed,TOP,MS)
setwd("/cluster/projects/p33/users/maxk/MS/")
final=final[!is.na(final$age),]
final = final %>% filter(age < 100)
training = final[(final$eid %in% unique(final$eid)),] %>% filter(diagnosis=="HC")
MS_validation = final %>% filter(diagnosis == "MS")
#
# we want to match the MS data
# for that, we need a matching list
match_frame = MS_validation %>% select(eid,age,sex,diagnosis)
match_training = training %>% select(eid,age,sex,diagnosis)
match_training = rbind(match_frame,match_training)
m.out = matchit(as.factor(diagnosis) ~ age + sex, data = match_training,
                 method = "nearest", distance = "glm")
matched_samples = match.data(m.out)
MS_validation=training[training$eid %in% matched_samples$eid,] # contains MS patients matched controls
MS_validation = rbind(MS_validation, final %>% filter(diagnosis == "MS"))
training=training[!training$eid %in% matched_samples$eid,] # training sample without MS validation set
validation = training %>%
  #group_by(age) %>%
  sample_frac(size=.1)
training=training[!training$eid %in% validation$eid,]
write.csv(training,"training_set.csv")
write.csv(validation,"HC_validation_set.csv")
write.csv(MS_validation,"MS_validation_set.csv")
write.csv(long_validation,"long_validation_set.csv")
#
print("All done.")
#
#
#