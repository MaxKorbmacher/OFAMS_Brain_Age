library(haven)
library(rlist)
library(dplyr)
library(tidyr)
# load demographics
test=read_sas("/Users/max/Documents/Local/MS/data/screening.sas7bdat")
demg = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/MRI_A_E_D.sav")
demo = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/Demographics.sav") # Demographics for inclusion age
icv1 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats/aparc_stats_area_lh_20210816145055.csv")
icv2 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220214_update/manual_edits_batch_20210913/aparc_stats_area_lh_20220214103330.csv")
#icv3 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_wmedits_rerun_20211201/aparc_stats_area_lh_20220209170109.csv")
icv4 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_edits_batch_20210913/aparc_stats_area_lh_20220209165901.csv")
icv5 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210127_stats/aparc_stats_area_lh_20210210092415.csv")
icv6 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210304/stats/aparc_stats_area_lh_20210315160653.csv")
#
# load FS stats
filenames = list.files("/Users/max/Documents/Local/MS/tabular/", pattern = ".csv", full.names = T)
ldf = lapply(filenames, read.csv)
ldf = list.rbind(ldf)
#
# save the merged FreeSurfer table
write.csv(ldf,"/Users/max/Documents/Local/MS/data/FSdata.csv",row.names = F)
#
# add age to the demographics
test$age = as.numeric(difftime(test$SC_DATEOFVISIT,test$DATEOFBIRTH)/365)
# demo$patno=demo$Patnr
# test = merge(test,demo,by="patno")
# test %>% filter(age=NA)
# test$age
# nrow(demo)
# make an anonymous data frame which will be used to match the validation set to the study sample (MS patients)
t1 = test %>% select(age,Sex)
names(t1)=c("age","sex")
t1$diagnosis = "MS"
write.csv(t1,"/Users/max/Documents/Local/MS/data/match_frame.csv",row.names = F)
#
#
#
#
########### AGE CALC
#tmp=ldf%>% filter(grepl("baseline",eid)) # select the ids we have available
#
tmp=ldf # create tmp df
session=ifelse(grepl("m24",tmp$eid),"24",tmp$eid)
session=ifelse(grepl("m120",session),"120",session)
session=ifelse(grepl("m12",session),"12",session)
session=ifelse(grepl("m11",session),"11",session)
session=ifelse(grepl("m10",session),"10",session)
session=ifelse(grepl("m9",session),"9",session)
session=ifelse(grepl("m8",session),"8",session)
session=ifelse(grepl("m7",session),"7",session)
session=ifelse(grepl("m6",session),"6",session)
session=ifelse(grepl("m5",session),"5",session)
session=ifelse(grepl("m4",session),"4",session)
session=ifelse(grepl("m3",session),"3",session)
session=ifelse(grepl("m2",session),"2",session)
session=ifelse(grepl("m1",session),"1",session)
tmp$session=as.numeric(ifelse(grepl("baseline",session),"0",session))
tmp$eid = as.numeric(gsub("\\_.*","",gsub("sub-","",tmp$eid))) # standardize ids
test$eid = test$patno
demg$eid=demg$Patnr # equal to test$eid for merging
#
# identify age for each visit using demg$Visitnr and demg$Visit_date
months = levels(factor(demg$Visitnr))
age_list = list()
for (i in 1:length(months)){
  #new = full_join(demg %>% filter(Visitnr==months[i]),test, by = c("eid","Visitnr"))
  new = demg %>% filter(Visitnr==months[i])
  new = merge(new,test,by="eid")
  new$age = as.numeric(difftime(new$Visit_date,new$DATEOFBIRTH)/365)
  age_list[[i]] = new %>% select(eid,Visitnr,age)
  #print(new$age)
}
age_list=list.rbind(age_list)
age_list$session = as.numeric(age_list$Visitnr)
new = full_join(tmp,age_list)
#
# quick and dirty solution to make sure all age fields are filled
## first, remove NA age for baseline (unfortunately, these guys are lost, there is no way of recovering this data)
new = new %>% filter(!(eid %in% (new %>% filter(session==0 & is.na(age)==T) %>% select(eid))$eid & session==0))
#new %>% group_by(session) %>% summarize(NA_sum = sum(is.na(age)), Rows = length(age))
# now, we add the number of months at each visit to the baseline age (in years)
# this is not the perfect solution, but a feasible solution.
for (i in 1:length(new$age)){
  if (is.na(new$age[i])==TRUE){
    if (new %>% filter(eid == eid[i] & session == 0) %>% nrow() > 0){
      a0 = new %>% filter(eid == eid[i] & session == 0) %>% select(age)
      new$age[i] = as.numeric(a0$age)+new[i,]$session/12
      }}}
new = new %>% select(-Visitnr)
new = na.omit(new)
#
# save csv
write.csv(new,"/Users/max/Documents/Local/MS/data/final_dat.csv")
#
# finally make a csv with total volume values
icv_proc = function(icv_df){
  icv = icv_df
  icv$eid = icv[,1]
  icv = icv %>% select(eid,BrainSegVolNotVent,eTIV)
  session=ifelse(grepl("m24",icv$eid),"24",icv$eid)
  session=ifelse(grepl("m120",session),"120",session)
  session=ifelse(grepl("m12",session),"12",session)
  session=ifelse(grepl("m11",session),"11",session)
  session=ifelse(grepl("m10",session),"10",session)
  session=ifelse(grepl("m9",session),"9",session)
  session=ifelse(grepl("m8",session),"8",session)
  session=ifelse(grepl("m7",session),"7",session)
  session=ifelse(grepl("m6",session),"6",session)
  session=ifelse(grepl("m5",session),"5",session)
  session=ifelse(grepl("m4",session),"4",session)
  session=ifelse(grepl("m3",session),"3",session)
  session=ifelse(grepl("m2",session),"2",session)
  session=ifelse(grepl("m1",session),"1",session)
  icv$session=as.numeric(ifelse(grepl("baseline",session),"0",session))
  icv$eid = as.numeric(gsub("\\_.*","",gsub("sub-","",icv$eid))) # standardize ids
  names(icv) = c("eid", "TotalVol", "TIV", "session")
  return(icv)
}
icv = rbind(icv_proc(icv1),icv_proc(icv4),icv_proc(icv5),icv_proc(icv6)) # icv_proc(icv2),icv_proc(icv3) # values not included
write.csv(icv,"/Users/max/Documents/Local/MS/data/icv.csv",row.names = F)