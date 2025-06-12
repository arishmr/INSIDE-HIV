# INSIDE-HIV
## Neuroimmunometabolic alterations and severity of depressive symptoms in people with HIV: An exploratory diffusion-weighted MRS study

This code repository contains code used for analysis for the INSIDE-HIV study, carried out at Brighton & Sussex Medical School (BSMS), University of Sussex, Brighton, UK. The objective of this study was to explore the neuroimaging signatures of severe and persistent depression among people with HIV.

## ASSOCIATED PUBLICATION
Mudra Rakshasa-Loots A, Diteko G, Dowell NG, Ronen I, Vera JH. Neuroimmunometabolic alterations and severity of depressive symptoms in people with HIV: An exploratory diffusion-weighted MRS study. Brain and Neuroscience Advances. 2025 Apr;9:23982128251335792. DOI: [10.1177/23982128251335792](https://doi.org/10.1177/23982128251335792)

## DATA ACCESS
Data associated with this study is publicly available in processed and de-identified form via the Edinburgh DataShare service at [this link](https://doi.org/10.7488/ds/7843). Due to the small sample size of the study, this publicly available data set does not include certain sensitive demographic variables. Researchers can request to access the full de-identified data set (including demographic variables) via the Edinburgh DataVault service at [this link](https://doi.org/10.7488/5d6256a4-2132-400d-84f8-f6445cc624f8). 

## CONTACT
For any questions about the code or dataset, please contact the lead investigator: Dr Arish Mudra Rakshasa-Loots (arish.mrl@ed.ac.uk).


## SCRIPTS

Data analysis is carried out using the following scripts:

### Section 1: Preamble
#### Always run these when starting the analysis
**01-packages**
Load (and install, if not already installed) the required packages for analysis

**02-load-data**
Load the dataset (which can be downloaded from the DOI provided above)

### Section 2: Demographics
**03-demos**
Summarise participant characteristics overall and by depression group

**04-summarystats**
Summarise depression scores and biomarker concentrations overall and by depression group

### Section 3: Neuroimaging biomarkers
**05-plots-by-score**
Visualise neuroimaging biomarker data by PHQ-9 scores

**06-correlations**
Calculate correlations of neuroimaging biomarkers with PHQ-9 scores

**06-adjusted-corr**
Calculate partial correlations of neuroimaging biomarkers with PHQ-9 scores, adjusted for relevant sociodemographic factors

**07-plots-by-group**
Visualise neuroimaging biomarker data by PHQ-9 group (high vs low depression)

**08-group-comparisons**
Statistical comparison of neuroimaging biomarker data between PHQ-9 groups (high vs low depression)

### Section 4: Blood biomarkers
**09-blood-plots-by-score**
Visualise blood biomarker data by PHQ-9 scores

**10-blood-correlations**
Calculate correlations of blood biomarkers with PHQ-9 scores

**11-blood-plots-by-group**
Visualise blood biomarker data by PHQ-9 group (high vs low depression)

**12-blood-group-comparisons**
Statistical comparison of blood biomarker data between PHQ-9 groups (high vs low depression)

**12-blood-imaging**
Calculate and visualise correlations of blood biomarkers with neuroimaging biomarkers


## VARIABLES

Variable name | Description
--- | ---
_Basic_ | 
`PID` | Participant ID
`Group` | Participant group (high or low depression severity)
`PHQ` | Total PHQ-9 score (calculated as sum of PHQ1:PHQ9)
_Demographic_ | 
`age` | Participant age in years
`sex` | Self-identified sex
`gender` | Self-identified gender (multiple selections possible)
`ethnicity` | Self-identified ethnicity (multiple selections possible)
`edu` | Self-reported number of years of education
`smoking` | Have you smoked cigarettes in the past 3 months?
`alcohol` | Have you drank alcohol in the past 3 months?
`recdrugs` | Have you used any recreational drugs in the past 3 months?
_Clinical_ | 
`curr.antidep` | Are you currently taking antidepressant medication?
`ever.antidep` | Have you ever taken antidepressant medication, now or in the past?
`curr.therapy` | Are you currently in therapy or counselling for depression?
`ever.therapy` | Have you ever been in therapy or counselling for depression, now or in the past?
`ViralLoad` | Latest documented HIV RNA viral load
`CD4` | Latest document absolute CD4 cell count
_Depression_ | 
`PHQ1`:`PHQ9` | Items on the PHQ-9 (asking about frequency IN THE PAST TWO WEEKS)
`PHQ1` | Little interest or pleasure in doing things
`PHQ2` | Feeling down, depressed, or hopeless
`PHQ3` | Trouble falling or staying asleep, or sleeping too much
`PHQ4` | Feeling tired or having little energy
`PHQ5` | Poor appetite or overeating
`PHQ6` | Feeling bad about yourself
`PHQ7` | Trouble concentrating on things
`PHQ8` | Moving or speaking so slowly that other people could have noticed? Or the opposite
`PHQ9` | Thoughts that you would be better off dead or of hurting yourself in some way
_MRS biomarkers_ | 
`Asp`:`Cho.Cr` | MRS biomarkers - those with suffix "`.Cr`" in the variable name are creatine referenced, those without the suffix are water referenced
`Asp` | Aspartate
`Cr` | Creatine
`Glc` | Glucose
`Glu.Gln` | Glutamate + Glutamine
`Ins` | Myo-inositol
`Lac` | Lactate
`NAA` | N-Acetyl Aspartate
`fGM` | fractional Gray Matter volume
`fWM` | fractional White Matter volume
`fCSF` | fractional CerebroSpinal Fluid volume
`Scandate` | Date of neuroimaging scan
_DW-MRS biomarkers_ | 
`ADC.NAA` | Apparent Diffusion Coefficient for NAA (DW-MRS)
`ADC.Cr` | Apparent Diffusion Coefficient for creatine (DW-MRS)
`ADC.Cho` | Apparent Diffusion Coefficient for choline (DW-MRS)
_DCE-MRI biomarkers_ | 
`K.trans_WB` | Volume transfer constant in whole brain (DCE-MRI)
`K.trans_thal` | Volume transfer constant in thalamus (DCE-MRI)
_Blood biomarkers_ | 
`BDNF`:`YKL40` | Inflammatory and other soluble biomarkers measured in blood serum, each variable name is the commonly used acronym for the respective biomarker


## NOTES

- My code is set up so that the 01-packages and 02-load-data scripts need to be run once each time when the project is opened, but not necessarily again.
- I export figures and results files to a /Figures and /Results sub-folder within my working directory, respectively. 
