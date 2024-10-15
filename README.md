# OFAMS longitudinal analysis of brain age
Analyses on the longitudinal OFAMS (ω-3 Fatty Acid Treatment in Multiple Sclerosis (MS)) MRI data using Brain Age

## Background
The study is based on a clinical trial where omega 3 fatty acids were provided to MS patients. While the omega 3 treatment was not successful, the patients were frequently followed up, every month during the first year of the study (baseline, m1-m12), and then more sporadically, at m24 and m120 (10-year follow-up). All patients received  interferon beta-1a treatment, which seemed successful based on the MRI data.

Original study MRI protocol (T1w-only):
T1 (TR 500-750 mm, TE 10-20 ms, matrix size 384x512, NEX 2, FOV 250 mm2, slice thickness 4 mm) and T2 weighted fast spin echo (TSE/FSE) (TR 4661 ms, TE 100 ms, matrix size 384x512, NEX 3, FOV 250 mm2, slice thickness 4 mm), in addition to sagittal 3D T1 weighted spoiled gradient echo (FFE/FLASH) (TR 20 ms, TE 4·6 ms, flip angle 25, matrix size 256, NEX 1, FOV 250 mm2, slice thickness 1 mm). All patients were given 15 ml contrast medium containing gadolinium.

Original report of the clinical trial: 
Torkildsen, Ø., Wergeland, S., Bakke, S., Beiske, A. G., Bjerve, K. S., Hovdal, H., ... & Myhr, K. M. (2012). ω-3 fatty acid treatment in multiple sclerosis (OFAMS Study): a randomized, double-blind, placebo-controlled trial. Archives of neurology, 69(8), 1044-1051. https://doi.org/10.1001/archneurol.2012.283

Finally, there is one previously published study looking at brain age in MS patients https://doi.org/10.3389/fneur.2019.00450
Limitations here: Healthy controls had no follow-ups, explainability of the model, relatively small N in training set, no further validations in an external validation set.
Advancements: brain age particularly sensitive to MS when estimated from subcortical regions and cerebellum.

## MRI acquisition 
A Gadolinium tracer was used during MRI (T1-weighted) aquisition. Most of the acquisitions were at relatively low resolution with the exception of fast low angle shot (FLASH) acquisitions at baseline, m6, m24, and m120. In some cases, due to adverse events, patients were scanned ±1 month from the scheduled appointment.

T2-weighted and T1-weighted gadolinium enhanced MRI scans were conducted using standard head coil with a 1·5 Tesla MRI unit. Imaging at the 10-year follow-up visit was performed at
different study sites, on a 3-Tesla (T) MRI scanner if available, alternatively using a 1.5 T MRI scanner, with a standard head coil. More information on the scanner parameters can be found at https://doi.org/10.1007/s00330-021-08405-8

## MRI processing
The data were then reconstructed using FreeSurfer version 7.1.1, and the tabularised data averaging across the parcels of the Desikan-Kelliany atlas used as input data for a trained machine learning brain age model to predict brain age.

The model was previously trained on N = 60,000 scans from healthy controls from 7 independent datasets: ABCD, AddNeuroMed, HBN, HCP, Rockland, TOP, UKB.
Validation steps:
- validation in N = XYZ healthy controls (cross-sectionally). (XYZ datasets).
- validation in N = 748 MS patients (cross-sectionally. (A single dataset from the Oslo University Hospital, Norway).
- validation in N = 27 healthy controls with a total of 1206 scans. (Two datasets ADNI and BBSC).

The cross-sectional healthy control validation sample was age-matched to the current MS sample (using the baseline age).

