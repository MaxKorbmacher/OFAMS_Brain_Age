echo This is a list of commands which can be run using the merge.py script.
echo The folders used in the list might be subject to change.
echo Hence, define them yourself. Double check. Be safe.
echo #
echo Note, there are two waves in the study and the data were processed in several stages.
echo This is why we use these many steps.
echo #
echo #
echo STUDY HALF I
echo load all available data
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/output/" "20220209165345" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo load first update
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_edits_batch_20210913/" "20220209165901" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo load second updated data set
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_wmedits_rerun_20211201/" "20220209170109" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo add remaining data set
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220214_update/manual_edits_batch_20210913/" "20220214103330" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo #
echo #
echo STUDY HALF II
echo load all available data from first FS run
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210127_stats/" "20210210092415" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo then after correction
python3 merge.py "/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210304/stats/" "20210315160653" "/Users/max/Documents/Local/MS/tabular/"
echo #
echo You might need to add more data folder or processing, based on data which will be processed.
echo #
echo And that is it. 
