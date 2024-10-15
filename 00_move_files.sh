echo Copy files into the correct folders.
echo This concerns only m1 m3 m9 m12.
echo #
cd /Users/max/Documents/Local/freesurfer/raw
echo Start copying.
for i in m1 m3 m9 m12; do mv /Users/max/Documents/Local/freesurfer/raw/bids_approved/*/${i}/*T1w.nii.gz /Users/max/Documents/Local/freesurfer/${i}/raw/; done
echo All done
