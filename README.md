# OFAMS longitudinal analysis of brain age
Analyses on the longitudinal OFAMS (ω-3 Fatty Acid Treatment in Multiple Sclerosis (MS)) MRI data using Brain Age

We start with explaining the brainage model and how apply it to your own data.

## The utilized brain age model: How to apply it to your own data?
1. After the FreeSurfer recon-all call or similar, e.g. using FastSurfer, tabular information are available per subject. We use the averages of the parcels of the Desikan-Killiany atlas as specified in the FreeSurfer output per subject.
2. Whether a brain age shall be predicted for a single subject or several, the obtained tables can be merged to a table per metric (e.g., thickness) using stats2table_bash.sh
3. Now, the resulting tables need to be merged (again!) into a single table. There is a provided merge.py script in the repository that can be used for that purpose. You end up with tables which can be used as input to predict age / to obtain the desired brain age predictions. Simply add "age" as a variable and you will also obtain corrected brain age predictions. The data frame can contain additional variables which will not be used by the model. The output file will contain all input data + the brain age predictions.
4. Run 'Rscript 05_predict.R /path/to/your/data.csv' and you will obtain the predictions in the working directory. You can add additional arguments to change paths and the output file name.


## OFAMS study background
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

## Model training and validation
The model was trained on N = 58,317 scans from healthy controls from 7 independent datasets: ABCD, AddNeuroMed, HBN, HCP, Rockland, TOP, UKB.
Validation steps:
 a) validation in N = 6,608 healthy controls (cross-sectionally). (7 independent datasets).
 b) validation in N = 748 MS patients and N = 751 matched healthy controls (cross-sectionally. (A single patient dataset from the Oslo University Hospital, Norway, and participants from six datasets to match the MS patients: AddNeuroMed, HBN, HCP, Rockland, TOP, UKB).
 c) validation in N = 3 healthy controls with a total of 103 scans. (BBSC: Bergen Breakfast Scanning Club).
 d) validation in N = 24 adults in their senescence transitioning from HC to MCI to AD with a total of 1103 scans. (ADNI: Alzheimer's Disease Neuroimaging Initiative).

**a) Validation in healthy controls**

A simple generalized additive model was trained with each input feature’s relationship with age being modelled as a spline with k = 4 possible knots. Uncorrected brain age predictions were strongly correlated with age (r = 0.9681, MAE=4.849 years, RMSE=6.279 years), and test sample age bias corrected predictions showed excellent age associations (r = 0.9999, MAE=1.333 years, RMSE=1.552 years; see Supplement 999 for training and validation performance metrics).

**b) Validation in multiple sclerosis (MS) patients matched controls**

Brain age predictions were more accurate in HC (r=0.702, MAE=7.335 years, RMSE=9.259 years), compared to MS (r=0.461, MAE=12.554 years, RMSE=14.649 years). The raw difference between MS and HC BAGs were 0.406 years, d=0.69, 95% CI[0.58,0.79], p=4.278*10-38, indicated by a two-samples t-test. A linear regression, with sex (betafemales=0.157±0.028 years, p=3.713*10-8), age (beta=0.029±0.001 years,p=6.851*10-73), and scanner (showing multiple significant differences at p<.05) as covariates, also indicated a higher BAG in (beta=0.534±0.075 years, p=2.536*10-12). Due to the observed sex-effect, we added the sex-diagnosis interaction to the model, showing a significant effect (betasex-diag=0.155±0.065, p=0.016), which can be interpreted as a larger difference between MS and HC in males (BAGHC= -0.58, 95%PI[-0.74,-0.42], BAGMS=-0.05, 95%PI[-0.11,0.02]) compared to females (BAGHC= -0.42, 95%PI[-0.58,-0.27], BAGMS=0.11, 95%PI[0.06,0.16]), as indicated by marginal effects and respective prediction intervals, holding age and scanner constant (age=38.50 years, scanner=GE750; see also Fig.1a).

**c) Validation in repeated measures of healthy midlife subjects**

N=3 midlife subjects aged 27.75, 30.38, and 40.54 years at baseline were followed over a period of 1.5 years leading to 103 scans. Bias-corrected brain age was near-perfect related with age (betaage=1.06 years) using linear random slopes at the level of the individual (Fig.1b), or within-subject age-correlations (rsub1=0.991, rsub2=0.983, rsub3=0.965). Age-correlations with uncorrected brain age estimates were weaker (rsub1=0.275, rsub2=0.272, rsub3=0.462), yet stronger than in a previous study using another model on the same data set. [1]
 
**d) Validation in repeated measures of elderly subjects transitioning to Alzheimer’s Disease**

In N = 38 elderly adults (baseline: mean=81.63, min=75.70, max=87.40 years of age) with 1872 scans, where subjects were followed over an average time of 8.62 years (min=1.2 years, max=15.2), marginal means by the best fitting model (marginal R2=3.36%, conditional R2=97.41%, see Supplement), indicated the lowest BAG in healthy controls (BAGmales=0.42, BAGmales=0.72), in contrast to mild cognitive impairment (BAGmales=0.61, BAGmales=0.90) and Alzheimer’s Disease (BAGmales=0.73, BAGmales=1.03) at the age of 85.45. For ageing trajectories see Fig.1c-d. 

**Evaluation of the model**

The ageing trajectories indicate that our brain age models can both represent diverse samples’ healthy ageing and simultaneously discriminate individuals with a disorder by an indicated higher brain age, i.e. higher modelling error. The steeper/quicker increase of brain age in mild cognitive impairment (Fig.999) indicates that our model might also be used to evaluate longitudinal brain ageing trajectories.
 
**References**

[1]	Leonardsen, E. H. et al. Deep neural networks learn general and clinically relevant representations of the ageing brain. NeuroImage 256, 119210 (2022).


