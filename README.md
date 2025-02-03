# Analysis of multimodal MS clinical trial data: the OFAMS longitudinal study
Analyses on the longitudinal OFAMS (ω-3 Fatty Acid Treatment in Multiple Sclerosis (MS)) entails MRI data using Brain Age, multiomics, clinical and behavioural data.

We start with explaining the _cortical_ brainage model and how to apply it to your own data.

Then, more OFAMS data/study-specific details are being described.

## The utilized brain age model: How to apply it to your own data?
1. After the FreeSurfer recon-all call or similar, e.g. using FastSurfer, tabular information are available per subject. We use the averages of the parcels of the Desikan-Killiany atlas as specified in the FreeSurfer output per subject.
2. Whether a brain age shall be predicted for a single subject or several, the obtained tables can be merged to a table per metric (e.g., thickness) using stats2table_bash.sh
3. Now, the resulting tables need to be merged (again!) into a single table. There is a provided merge.py script in the repository that can be used for that purpose. You end up with tables which can be used as input to predict age / to obtain the desired brain age predictions. Simply add "age" as a variable and you will also obtain corrected brain age predictions. The data frame can contain additional variables which will not be used by the model. The output file will contain all input data + the brain age predictions.
4. Download the model (sim_model.rda) from https://doi.org/10.17605/OSF.IO/AE5RX
6. Run 'Rscript 05_predict.R /path/to/your/data.csv' and you will obtain the predictions in the working directory. You can add additional arguments to change paths and the output file name.


## OFAMS study background
The study is based on a clinical trial where omega 3 fatty acids were provided to MS patients. While the omega 3 treatment was not successful, the patients were frequently followed up, every month during the first year of the study (baseline, m1-m12), and then more sporadically, at m24 and m120 (10-year follow-up). All patients received  interferon beta-1a treatment, which seemed successful based on the MRI data.

Original study MRI protocol (T1w-only):
T1 (TR 500-750 mm, TE 10-20 ms, matrix size 384x512, NEX 2, FOV 250 mm2, slice thickness 4 mm) and T2 weighted fast spin echo (TSE/FSE) (TR 4661 ms, TE 100 ms, matrix size 384x512, NEX 3, FOV 250 mm2, slice thickness 4 mm), in addition to sagittal 3D T1 weighted spoiled gradient echo (FFE/FLASH) (TR 20 ms, TE 4·6 ms, flip angle 25, matrix size 256, NEX 1, FOV 250 mm2, slice thickness 1 mm). All patients were given 15 ml contrast medium containing gadolinium.

Original report of the clinical trial: 
Torkildsen, Ø., Wergeland, S., Bakke, S., Beiske, A. G., Bjerve, K. S., Hovdal, H., ... & Myhr, K. M. (2012). ω-3 fatty acid treatment in multiple sclerosis (OFAMS Study): a randomized, double-blind, placebo-controlled trial. Archives of neurology, 69(8), 1044-1051. https://doi.org/10.1001/archneurol.2012.283

Finally, there is one previously published study looking at longitudinal brain age in MS patients https://doi.org/10.3389/fneur.2019.00450
Limitations here: Healthy controls had no follow-ups, explainability of the model, relatively small N in training set, no further validations in an external validation set.
Advancements: brain age particularly sensitive to MS when estimated from subcortical regions and cerebellum.

## MRI acquisition 
A Gadolinium tracer was used during MRI (T1-weighted) aquisition. Most of the acquisitions were at relatively low resolution with the exception of fast low angle shot (FLASH) acquisitions at baseline, m6, m24, and m120. In some cases, due to adverse events, patients were scanned ±1 month from the scheduled appointment.

T2-weighted and T1-weighted gadolinium enhanced MRI scans were conducted using standard head coil with a 1·5 Tesla MRI unit. Imaging at the 10-year follow-up visit was performed at
different study sites, on a 3-Tesla (T) MRI scanner if available, alternatively using a 1.5 T MRI scanner, with a standard head coil. More information on the scanner parameters can be found at https://doi.org/10.1007/s00330-021-08405-8

## MRI processing
The data were then reconstructed using FreeSurfer version 7.1.1, and the tabularised data averaging across the parcels of the Desikan-Kelliany atlas used as input data for a trained machine learning brain age model to predict brain age.

## Model training and validation steps: overview
The model was trained on N = 58,317 scans from healthy controls from 7 independent datasets: ABCD, AddNeuroMed, HBN, HCP, Rockland, TOP, UKB.
Validation steps:
 a) validation in N = 6,608 healthy controls (cross-sectionally). (7 independent datasets).
 b) validation in N = 748 MS patients and N = 751 matched healthy controls (cross-sectionally. (A single patient dataset from the Oslo University Hospital, Norway, and participants from six datasets to match the MS patients: AddNeuroMed, HBN, HCP, Rockland, TOP, UKB).
 c) validation in N = 3 healthy controls with a total of 103 scans. (BBSC: Bergen Breakfast Scanning Club).

## Training
Fit in training data (N=150,517) and validation data (N = 6,608, mean age = 49.49±25.04, range: 5.28-86.7) for corrected brain age (c) and uncorrected brain age (u). Corrected brain age, is the individual level brain age where the training sample bias was regressed out (see: https://doi.org/10.1002/hbm.25837).

All training and validation data were obtained from healthy controls! Fit metrics are comparable to other models trained on similar data and corrected brain age estimates produce naturally better fit indices. However, the advantage of the presented models is their generalizability to other, external, datasets (see below) and their explainablity, since the models have a simple architecture of added splines which allowing polynomials up to the 4th order / k=4 knots.

We also want to highlight that hemisphere-specific models perform similar to models of both hemispheres (also highlighted previously: https://doi.org/10.1038/s41467-024-45282-3).

|    Sample and BA    | Pearson's r	|   R2   |	 MAE  |	 RMSE  |
| :---------: |   :---------: |  :---: |  :---: |  :---: |
|Training u   |   0.9568|0.9154|5.0935|6.5013|
|Training c   |   0.9634|0.9281|4.8933|6.2186|
|Test u |   0.9622|0.9259|5.2827|6.8292|
|Test c |   0.9680|0.9370|5.2034|6.6508|

(R2 = Variance explained, MAE = Mean Absolute Error, RMSE = Root Mean Squared Error)

## Validation
Fit in external validation data, which are all healthy controls (N = 751, mean age=38.83±9.77, range: 18.63-87.5):

| Sample & Correction | Pearson's r	|   R2   |	 MAE  |	 RMSE  |
|  :----------------------: | :---------: |  :---: |  :---: |  :---: |
|   Validation u    | 0.7561|0.5717|5.9837|7.5362|
|    Validation c    | 0.7849|0.6161|5.9662|7.4548|

An example of data fit for a single disease group, here multiple sclerosis patients (N = 748, mean age=38.63±9.46, range: 18.46-70.34):

| Sample & Correction | Pearson's r	|   R2   |	 MAE  |	 RMSE  |
|  :----------------------: | :---------: |  :---: |  :---: |  :---: |
|    MS u    | 0.4626|0.2140|10.4083|12.4970|
|    MS c    | 0.5357|0.2870|9.5808|11.5617|

Finally, we were interested in how our model performs in longitudinal data looking at three healthy individuals scanned in total 103 times over a 1.5 years period (mean age per subject: 30.66, 28.09, 40.66).
|Subject|	Uncorrected BA Pearson's r|	p	| Corrected BA Pearson's r|	p	|
| :---: |:---: |:---: |:---: |:---: |
|sub-1 u	|0.2232	|0.1780|	0.2466|	0.1356	|
|sub-2 u	|0.2086|	0.1965|0.2290|	0.1552	|
|sub-3 u	|0.5236|	0.0072	|0.5310|	0.0063	|


Longitudinal results suggest better individual-level model fit than previous models (e.g., https://doi.org/10.1002/brb3.3219).

Interrim conclusion from the model validation: We recommend using corrected brain age estimates, especially in the case of longitudinal analyses.

Another reasoning for using corrected estimates is boosting effect sizes for group differences. Let's take the mentioned MS sample as an example again:
As an example how corrected brain age estimates boost predictions:
Corrected BAG differences between MS (mean age = 38.63±9.46 years) and HC (mean age = 38.73±9.61) are more clearly expressed when correcting for the training age bias:
See for that first the corrected BAG differences between MS and HC (paired samples t-tests):

| BAG | Cohen's d & 95% CI | t(df) | p|
|:---:|:---:|:---:|:---:|
|BAG c| d=-1.00[-1.11;-0.89] | t(1462)=-19.40| p<2.2*10^-16|
|BAG u| d=-0.89[-1.15;-0.79]| t(1466.7)=-17.29| p<2.2*10^-16|

We see clear differences between the outlined group differences, where corrected brain ages produce larger effect sizes differentiating healthy controls and MS patients.

## Explainability
An important advantage of the utilzed GAM is that the model coefficients can be fairly easily interpreted using ANOVA test, the respective degrees of freedom and by knocking out regions by setting all feature values of a specific region to 0 or varying the values for each region while holding the other regions constant. See a summary of the findings below. The data used for the visualisation are the training data. However, lesion-styled interpretability is also possible in any unseen data by simpy setting the values of a single region to zero.

Degrees of freedom per region and feature weight maps showing the non-linearity of the features and estimated model coefficiens, respectively
![imp_df](https://github.com/user-attachments/assets/01849aeb-f422-4e39-b712-df784aad3ed6)

Knockout importance (setting all features in a region = 0 across the entire sample)
![knock](https://github.com/user-attachments/assets/3984a13e-c650-4cf5-862c-cb0e365d866e)

Individual level variability (each feature can range from the minimum to the maximum observed value in the training sample in subjects with most accurate brain age predictions. The value for the region is correlated with the resulting brain age [Note: the results are constant across age groups]) 
![single_sub](https://github.com/user-attachments/assets/7a353fd5-b83f-4eb1-863d-04b163e23d00)


## Results
The results will soon be presented in a pre-print and following scientific publication.
Multiverse summary statistics can be found here: https://doi.org/10.17605/OSF.IO/AE5RX


