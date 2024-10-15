library(haven)
library(rlist)
# load demographics
test=read_sas("/Users/max/Documents/Local/screening.sas7bdat")
# load FS stats
filenames = list.files("/Users/max/Documents/Local/MS/tabular/", pattern = ".csv", full.names = T)
ldf = lapply(filenames, read.csv)
ldf = list.rbind(ldf)
#
# save the merged FreeSurfer table
write.csv(ldf,"/Users/max/Documents/Local/MS/data/FSdata.csv",row.names = F)
