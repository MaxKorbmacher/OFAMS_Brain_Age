# make command list for recon-all
# Author: Max Korbmacher
# Date: 14 October 2024
#
# WE START WITH TIME POINTS m1, m3, m9, and m12 (baseline, m6, m24, and m120 are already processed)
#
# function to list files and desired output
mkcmdlist = function(nii_foldername){
  files = list.files(path = paste("/Users/max/Documents/Local/freesurfer/",nii_foldername,"/raw/",sep=""), pattern = "T1w.nii.gz", full.names = T)
  out = list.files(path = paste("/Users/max/Documents/Local/freesurfer/",nii_foldername,"/raw/",sep=""), pattern = "T1w.nii.gz", full.names = F)
  out=gsub(".nii.gz","",out)
  cmd_list = c()
  for (i in 1:length(files)){
    cmd_list[i] = paste0("recon-all -i ", files[i], " -s /Users/max/Documents/Local/freesurfer/m12/recon-all/", out[i], " -all", sep = "")
  }
  write.table(cmd_list,paste("/Users/max/Documents/Local/MS/scripts/",nii_foldername,".txt",sep=""), row.names = F, quote = F, col.names = F)
}
# we apply the function to m1, m3, m9, and m12
mkcmdlist("m1")
mkcmdlist("m3")
#mkcmdlist("m9")
#mkcmdlist("m12")