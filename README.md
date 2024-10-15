# OFAMS longitudinal analysis of brain age
Analyses on the longitudinal OFAMS (ω-3 Fatty Acid Treatment in Multiple Sclerosis (MS)) MRI data using Brain Age

## Background
The study is based on a clinical trial where omega 3 fatty acids were provided to MS patients. While the omega 3 treatment was not successful, the patients were frequently followed up, every month during the first year of the study (baseline, m1-m12), and then more sporadically, at m24 and m120 (10-year follow-up). All patients received  interferon beta-1a treatment, which seemed successful based on the MRI data.

Original report of the clinical trial: 
Torkildsen, Ø., Wergeland, S., Bakke, S., Beiske, A. G., Bjerve, K. S., Hovdal, H., ... & Myhr, K. M. (2012). ω-3 fatty acid treatment in multiple sclerosis (OFAMS Study): a randomized, double-blind, placebo-controlled trial. Archives of neurology, 69(8), 1044-1051. https://doi.org/10.1001/archneurol.2012.283

## MRI data
A Gadolinium tracer was used during MRI (T1-weighted) aquisition. Most of the acquisitions were at relatively low resolution with the exception of fast low angle shot (FLASH) acquisitions at baseline, m6, m24, and m120. In some cases, due to adverse events, patients were scanned ±1 month from the scheduled appointment.

The data were then reconstructed using FreeSurfer version 7.1.1, and the tabularised data averaging across the parcels of the Desikan-Kelliany atlas used as input data for a trained machine learning brain age model to predict brain age.

The model was previously trained on N = 60,000 scans from healthy controls and validated in N = 10,000 scans of healthy controls. The validation sample was age-matched to the current MS sample (using the baseline age).
