# Simple script to merge FS output tables 
# Author: Max Korbmacher
# 14 October 2024
#
#######################
# The script can be run the following way:
# python3 merge.py "path/where/recon-all/output/tables/are" "suffix_to_filen_name/session_id'" "path/where/files/shall/be/saved/"
# output file will be saved according to the suffix, to identify the session id
#
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
from functools import reduce
import os
import sys
#
#paths
T1path=sys.argv[1]# example: "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/output/"
proc=sys.argv[2] # example: "20220209165345"
savepath=sys.argv[3] # example: "/Users/max/Documents/Local/MS/tabular/"
#
dict = {}
#
#Make list of MRI variables
MRI_list = [#"thicknessstd",
	"aparc_stats_thickness",
	"aparc_stats_volume",
	"aparc_stats_area"
]
hem_list = ["lh","rh"]
#List of MRI data frames
df_MRI_list = []
i = 0
for MRI in MRI_list:
		for hem in hem_list:
			df_MRI_list.append(pd.read_csv(T1path+'%s_%s_%s.csv' % (MRI,hem,proc), delimiter = ","))
			df_MRI_list[i] = df_MRI_list[i].rename(index=str, columns={"%s.%s" % (hem,MRI): 'eid'})
			df_MRI_list[i] = df_MRI_list[i].replace({'FS_': ''}, regex =True)
			df_MRI_list[i]['eid'] = df_MRI_list[i][df_MRI_list[i].columns[0]] #.values
			df_MRI_list[i].drop(columns=df_MRI_list[i].columns[0], axis=1,  inplace=True)
			df_MRI_list[i]['BrainSegVolNotVent'] = df_MRI_list[i]['eid'] # just to create these two, if not there
			df_MRI_list[i]['eTIV'] = df_MRI_list[i]['eid']
			df_MRI_list[i] = df_MRI_list[i].drop('BrainSegVolNotVent', axis = 1)
			df_MRI_list[i] = df_MRI_list[i].drop('eTIV', axis = 1)
			print (df_MRI_list[i].head(5))
			i += 1

#Merge all the MRI data frames in our df_MRI_list into one list. Merge based on eid
df_MRI_merged = reduce(lambda  left,right: pd.merge(left,right,on='eid'), df_MRI_list)


print ('mri columns:')
print (len(df_MRI_merged.columns))

#df_MRI_merged = df_MRI_merged.drop(['lh_WhiteSurfArea_area','rh_WhiteSurfArea_area','rh_???_volume','lh_???_volume', 'lh_???_area','rh_???_area','lh_???_thickness','rh_???_thickness'], axis=1)

print ('mri columns after drop:')
print (len(df_MRI_merged.columns))
print (df_MRI_merged.head(4))
print (df_MRI_merged.columns)
df=df_MRI_merged
print ('N in final file:')
print (len(df))
print (df.head(5))
print ("Writing output CSV.")
df.to_csv(savepath+proc+'.csv', sep=',',index=None)
