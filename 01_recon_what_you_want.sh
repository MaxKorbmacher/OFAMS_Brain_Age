echo Reconstruct brains using FreeSurfer in parallel
echo Created by Max Korbmacher, 14 October 2024
echo First, make lists to run parallel on
Rscript make_cmd_list.R
echo #
echo Second, run the commands in the lists in parallel
echo #
cat m*.txt | parallel
echo #
echo The end.