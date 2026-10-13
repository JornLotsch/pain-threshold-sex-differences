# Quantitative Sensory Testing Dataset - Preprocessing Pipeline

## Dataset: Quantitative Sensory Testing Dataset with Capsaicin and Menthol Sensitization Assessing Modality-Specific Sex Differences and Sensitization Effect Sizes in 125 Healthy Volunteers

Medeley dataset: "Quantitative sensory testing dataset assessing sex differences and sensitization in 125 healthy volunteers" (links and details, see below).

This repository contains the complete preprocessing pipeline for the quantitative sensory testing (QST) dataset accompanying the publication:

> **Quantitative sensory testing dataset assessing sex differences and sensitization in 125 healthy volunteers**
>
> [Citation to be added upon publication]
>
> **Repository:** https://github.com/JornLotsch/pain-threshold-sex-differences

## Overview

This repository contains the preprocessing pipeline to prepare the quantitative sensory testing dataset for machine learning applications. The original pipeline implements two stages: (1) natural-logarithmic transformation and sensitization effect calculation in Python, and (2) censored data imputation, standardization, and train-test splitting in R. A final R script (3) assembles the files deposited at Mendeley Data and produces quality-control output. The result is a fully processed dataset with standardized features ready for statistical modeling or machine learning analysis.

In addition, `R/verify_reproduction_from_deposited_data.R` regenerates `exp_pain_data_transformed.csv` and the training/validation assignment **entirely in R and directly from the deposited Mendeley files** (see [Verification of Reproducibility](#verification-of-reproducibility)).

## Dataset Description

This dataset contains experimental pain threshold measurements collected from 125 healthy adult volunteers under both baseline and chemically-sensitized conditions. The data were generated using established quantitative sensory testing (QST) protocols, including mechanical (punctate and blunt pressure), thermal (heat and cold), and electrical pain modalities. Sensitization was induced using topical capsaicin (for heat and punctate pressure) and menthol (for cold) to assess pain plasticity and sensitization mechanisms.

**Dataset Citation:**
> Lotsch, Jorn; Flühr, Karin; Neddermeyer, Till J. (2026), "Quantitative sensory testing dataset with capsaicin and menthol sensitization assessing modality-specific sex differences and sensitization effect sizes in 125 healthy volunteers", Mendeley Data, V1, doi: 10.17632/y6v3rtpgps.1

### Input Data

**Input:** Participant-level pain threshold estimates from 125 participants. Each value is a threshold estimate derived from repeated measurements or logistic regression, not a trial-level observation; trial-level observations are not part of the deposited dataset.
- 5 non-sensitized pain modalities: von Frey (punctate pressure), heat, cold, blunt pressure, electrical stimulation
- 3 sensitized measurements: von Frey + capsaicin, heat + capsaicin, cold + menthol
- Metadata: participant ID, age, sex, smoking status, experimenter identifier

**Units of the untransformed thresholds** (`exp_pain_data_orig.csv`):

| Variable | Unit | Theoretical range |
|---|---|---|
| vonFrey, vonFrey_Capsaicin | g | 0.008 – 300 g (filament set) |
| Heat, Heat_Capsaicin | °C | 32 – 50 °C |
| Cold, Cold_Menthol | °C | 0 – 32 °C |
| Pressure | N/cm² | > 0 |
| Current | mA | 0 – 20 mA |

### Output Data

**Output of the preprocessing (Stage 2):** Two files ready for machine learning
- `PainThresholds_scaled_Training.csv` - Training set (80% of data, n = 100; 45 female, 55 male)
- `PainThresholds_scaled_Test.csv` - Validation set (20% of data, n = 25; 11 female, 14 male)

Each output file contains:
- 11 standardized pain threshold features (8 transformed thresholds + 3 sensitization effects)
- Target variable (sex: 0=male, 1=female)

**Deposited files (Stage 3):** These two files are combined into `exp_pain_data_transformed.csv` (without the target), and the assignment is stored in the column `TrainingValidation` of `exp_pain_metadata.csv`. Variable names map as follows:

| Original analysis | Deposited files |
|---|---|
| vonFrey, vonFrey_C | vonFrey, vonFrey_Capsaicin |
| Hitze, Hitze_C | Heat, Heat_Capsaicin |
| Kaelte, Kaelte_M | Cold, Cold_Menthol |
| Druck | Pressure |
| Strom | Current |
| CapsHeat, CapsvFrey, MenthCold | Capsaicin_Effect_Heat, Capsaicin_Effect_vonFrey, Menthol_Effect_Cold |

## Pipeline Structure

### Stage 1: Data Transformation (Python)

**Script:** `Python/Read_explore_PainThresholds.py`

**Purpose:** Apply signed logarithmic transformation and calculate sensitization effects

**Steps:**
1. Read participant-level pain threshold estimates (`Phaeno_125.csv`)
2. Apply signed logarithmic transformation: `sign(x) × ln(|x + 1|)` (natural logarithm, `np.log`)
   - As all thresholds are ≥ 0, this equals `ln(x + 1)`; no sign reversal is applied
   - Handles zero values
   - Reduces skewness while preserving direction
3. Calculate sensitization effects (difference scores):
   - `CapsHeat` = Heat - Heat_Capsaicin
   - `CapsvFrey` = vonFrey - vonFrey_Capsaicin
   - `MenthCold` = Cold - Cold_Menthol
4. Write transformed data: `dfPainThresholdsAnalyzed_log_eff.csv`

The modules `explore_tukey_lop` and `box_and_heatplots` imported by the script are only used for visual exploration of distributions and are not required to produce the output file.

**Key transformations:**
- Signed log transformation normalizes skewed distributions
- Sensitization effects quantify how chemical sensitizers alter pain perception

### Stage 2: Data Preprocessing (R)

**Script:** `R/preprocessing_pain_Mendeley_mf_dataset.R`

**Purpose:** Handle censored data, standardize features, and create train/test split

**Steps:**

1. **Data Loading**
   - Load transformed pain thresholds from Stage 1
   - Combine with sex variable from the phenotype source file `Phaeno_125.csv`

2. **Censored Data Imputation**
   - **Problem:** Some measurements hit instrument limits
     - Cold thresholds: No pain at 0°C (left-censored, coded as 0; transformed: ln(1) = 0)
     - von Frey thresholds: No pain at 300 g maximum filament force (right-censored, coded as 300; transformed: ln(301))
   - **Solution:** Censored likelihood multiple imputation (`clmi`)
     - Identify predictor variables with |Spearman r| ≥ 0.4 (resulting predictors: Heat for cold; Heat and Pressure for von Frey)
     - Generate 5 imputed datasets per censored variable (`seed = 42`)
     - Take median across imputations
     - Recalculate sensitization effects using imputed values
   - **Sign reversal:** As `clmi()` handles left-censoring, the right-censored von Frey values are negated before imputation and negated back afterwards. The deposited von Frey variables are therefore **not** sign-reversed (higher values = higher thresholds).
   - **Implementation:** `clmi()` function from `lodi` package

3. **Z-score Standardization**
   - Transform all 11 features to mean=0, SD=1
   - Ensures all variables on comparable scales for machine learning

4. **Train/Test Split**
   - **Method:** Optimal Distribution-Preserving Downsampling (OPDIS)
   - **Rationale:** Maintains distributional properties of full dataset
   - **Parameters:**
     - 80% training, 20% validation
     - 1,000,000 optimization trials
     - Seed = 42 (reproducibility)
   - **Implementation:** `opdisDownsampling()` function

5. **Output**
   - Training set: 80% of participants, distribution-matched to full dataset
   - Validation set: Remaining 20%

**Exact order of preprocessing steps and implication for independence:**

1. ln(x + 1) transformation and sensitization effects (Python, deterministic, per subject)
2. Imputation of censored values, estimated on **all 125 subjects**
3. Z-score standardization, mean and SD estimated on **all 125 subjects**
4. Training/validation split (opdisDownsampling)

Because the imputation models and scaling parameters were estimated on the complete dataset before the split, the validation subset is **not fully independent** of the training subset. Users requiring a strictly independent validation should re-estimate imputation and scaling on the training subset only, starting from `exp_pain_data_orig.csv` and the `TrainingValidation` column of `exp_pain_metadata.csv`.

### Stage 3: Assembly of the Deposited Files and Quality Control (R)

**Script:** `R/prepare_pain_Mendeley_mf_dataset.R`

**Purpose:** Write the files deposited at Mendeley Data, create Figure 1 and the metadata plot, and report quality-control information

**Steps:**
1. Read `Phaeno_125.csv`, select and rename the threshold and metadata variables (English names)
2. Read `PainThresholds_scaled_Training.csv` and `PainThresholds_scaled_Test.csv`, rename, and sort by participant ID
3. Write `exp_pain_data_orig.csv`, `exp_pain_data_transformed.csv`, `exp_pain_metadata.csv`
4. Plot untransformed data (Figure 1, `plot_orig_data.svg/.png`) and metadata (`plot_metadata.svg/.png`)
5. Quality control: ID consistency across files, observed ranges, and counts of values at the theoretical limits (`qc_ranges_censored_counts.csv`), and comparison of training and validation sets (`qc_training_validation_comparison.csv`)

## Quality Control

Output of the quality-control section of `R/prepare_pain_Mendeley_mf_dataset.R`:

- **ID consistency:** 125 unique, non-missing IDs, identical and in identical order in all three deposited files
- **Training/validation:** 100 training (45 female, 55 male), 25 validation (11 female, 14 male)
- **Observed ranges and values at theoretical limits** (n = 125 per variable, no missing values):

| Variable | Unit | Observed min | Observed max | n at min | n at max |
|---|---|---|---|---|---|
| vonFrey | g | 2.60 | 300.0 | 0 | 38 (30.4%) |
| Heat | °C | 39.4 | 49.8 | 0 | 0 |
| Cold | °C | 0.0 | 28.1 | 30 (24.0%) | 0 |
| Pressure | N/cm² | 15.4 | 103.0 | 0 | 0 |
| Current | mA | 0.6 | 7.0 | 0 | 0 |
| vonFrey_Capsaicin | g | 1.61 | 300.0 | 0 | 21 (16.8%) |
| Heat_Capsaicin | °C | 32.0 | 48.3 | 1 (0.8%) | 0 |
| Cold_Menthol | °C | 0.0 | 30.3 | 19 (15.2%) | 0 |

The censored values of von Frey (at 300 g) and cold (at 0 °C) were imputed (Stage 2). The single Heat_Capsaicin value at the 32 °C baseline temperature was not imputed.

- **Dependence between training and validation sets** (`qc_training_validation_comparison.csv`): Imputation and standardization preceded the split, which was overlooked in the original analyses. Neither step used the class variable (sex), and imputation affected only 108 of 1,000 threshold values. Means of the z-scaled variables in the training set range from −0.04 to 0.01 and SDs from 0.99 to 1.04. Standardizing the validation set with training-only parameters would change its values by a median of 1.5% and at most 10.5% of one SD (Menthol_Effect_Cold; ≤ 6.0% for all other variables), without changing their rank order. No variable differs in distribution between the subsets (Wilcoxon and Kolmogorov-Smirnov tests, all p ≥ 0.54).
- **Strictly independent validation:** impute the censored values and estimate the standardization parameters in the training set only, then apply them to the validation set, starting from `exp_pain_data_orig.csv` and the `TrainingValidation` column of `exp_pain_metadata.csv`. `R/verify_reproduction_from_deposited_data.R` can be adapted for this purpose.

## Verification of Reproducibility

**Script:** `R/verify_reproduction_from_deposited_data.R`

This script uses only the three deposited Mendeley files and R. It re-implements the Python transformation (ln(x + 1), sensitization effects), the predictor selection, the `clmi` imputation, the z-scaling and the opdisDownsampling split with the settings of the original analysis, and compares the results with the deposited files.

Result (R 4.6.1, lodi 0.9.2, opdisDownsampling 1.7):
- `exp_pain_data_transformed.csv`: all 11 variables reproduced (maximum absolute difference < 1e-14)
- Training/validation assignment in `exp_pain_metadata.csv`: reproduced identically

Running `R/prepare_pain_Mendeley_mf_dataset.R` on the original analysis files reproduces the three deposited files byte-identically.

### Helper Functions

**Script:** `R/FindMostCorrelatedVariables2.R`

**Purpose:** Identify highly correlated variables for imputation predictor selection

**Features:**
- Computes all pairwise Spearman correlations
- Generates visualization plots with confidence intervals
- Highlights correlations involving target variables
- Used to select predictors (|r| ≥ 0.4) for censored data imputation

## Software Requirements and Installation

### Software Versions

The Python and R preprocessing stages were originally run in 2022; the exact software versions used at that time were not recorded. All stages were re-run with the versions listed below: the Python stage reproduces the original intermediate file `dfPainThresholdsAnalyzed_log_eff.csv` (maximum absolute deviation < 1e-15), the R preprocessing stage reproduces `PainThresholds_scaled_Training.csv` and `PainThresholds_scaled_Test.csv` byte-identically, and the assembly stage reproduces the three deposited files byte-identically.

### Python
- **Version:** re-run and verified with Python 3.12.12
- **Packages:** `numpy` 1.26.4, `pandas` 3.0.5
- **Install:** `pip install numpy pandas`

### R
- **Version:** re-run and verified with R 4.6.1, lodi 0.9.2, opdisDownsampling 1.7, inspectdf 0.0.12.1, ggplot2 4.0.3, dplyr 1.2.1
- **Required packages (Stage 2 and verification):** `lodi`, `opdisDownsampling`, `dplyr`, `parallel`, `inspectdf`, `ggplot2` (`inspectdf` and `ggplot2` are used by `FindMostCorrelatedVariables2()` for the predictor selection)
- **Required packages (Stage 3):** `ggplot2`, `reshape2`, `ggthemes`, `dplyr`, `tidyr`, `ggnewscale`
- **Install:**
```r
install.packages(c("lodi", "opdisDownsampling", "dplyr", "inspectdf", "ggplot2",
                   "reshape2", "ggthemes", "tidyr", "ggnewscale"))
```

## Usage

### Quick Start

Before running the pipeline:

1. **Download the dataset from Mendeley Data:** https://doi.org/10.17632/y6v3rtpgps.1
   - Download `exp_pain_data_orig.csv`, `exp_pain_data_transformed.csv`, `exp_pain_metadata.csv`

2. **Install dependencies:** Python 3.12.3 and R 4.5.2 with required packages (see Installation section)

3. **Update paths in scripts** to point to your downloaded data directory

### Regenerating the Transformed Data from the Deposited Files (recommended)

```r
setwd("R")
source("verify_reproduction_from_deposited_data.R")
```

Set `DATA_DIR` at the top of the script to the directory containing the downloaded files:
- `exp_pain_data_orig.csv`
- `exp_pain_data_transformed.csv`
- `exp_pain_metadata.csv`

The training/validation split (10⁶ trials) takes about one minute with 30 cores; reduce `MAX_CORES` if fewer are available.

### Running the Original Pipeline

The original pipeline reads the single-file source `Phaeno_125.csv` (see [Data Not Included](#data-not-included-in-this-repository)), which contains the same values as `exp_pain_data_orig.csv` and `exp_pain_metadata.csv` with German variable names.

1. **Python Stage:**
   ```bash
   python3 Python/Read_explore_PainThresholds.py
   ```

2. **R Stage:**
   ```r
   source("R/preprocessing_pain_Mendeley_mf_dataset.R")
   ```

3. **R Assembly and QC:**
   ```r
   source("R/prepare_pain_Mendeley_mf_dataset.R")
   ```

### Path Configuration

The paths are specific to the original analysis environment and must be adapted:

- **Python:** `os.chdir(...)`, `pfad_o` and `pfad_u1` in `Python/Read_explore_PainThresholds.py`
- **R, Stage 2:** `pfad_o` and `pfad_u` in `R/preprocessing_pain_Mendeley_mf_dataset.R`
- **R, Stage 3:** `INPUT_FILE_PATH`, `OUTPUT_DIR` and `BASE_PATH` in `R/prepare_pain_Mendeley_mf_dataset.R`
- **Verification:** `DATA_DIR` in `R/verify_reproduction_from_deposited_data.R`

### Input and Output Files per Script

| Script | Input | Output |
|---|---|---|
| `Python/Read_explore_PainThresholds.py` | `Phaeno_125.csv` | `dfPainThresholdsAnalyzed_log_eff.csv` |
| `R/preprocessing_pain_Mendeley_mf_dataset.R` | `Phaeno_125.csv`, `dfPainThresholdsAnalyzed_log_eff.csv` | `PainThresholds_scaled_Training.csv`, `PainThresholds_scaled_Test.csv` |
| `R/prepare_pain_Mendeley_mf_dataset.R` | `Phaeno_125.csv`, `PainThresholds_scaled_Training.csv`, `PainThresholds_scaled_Test.csv` | `exp_pain_data_orig.csv`, `exp_pain_data_transformed.csv`, `exp_pain_metadata.csv`, `plot_orig_data.*`, `plot_metadata.*`, `qc_ranges_censored_counts.csv`, `qc_training_validation_comparison.csv` |
| `R/verify_reproduction_from_deposited_data.R` | `exp_pain_data_orig.csv`, `exp_pain_data_transformed.csv`, `exp_pain_metadata.csv` | Console report of the comparison |

## Key Methodological Decisions

### Why Signed Logarithmic Transformation?
- Pain thresholds are highly right-skewed
- Log transformation normalizes distributions
- Signed version preserves information about direction (important for sensitization effects)

### Why Censored Likelihood Multiple Imputation?
- Censored data cannot simply be discarded (biases results)
- `clmi` uses correlated variables to predict censored values
- Accounts for uncertainty through multiple imputations
- Specifically designed for limit-of-detection censoring

### Why Correlation Threshold ≥ 0.4?
- Balances predictive power with model complexity
- Strong enough correlations to inform imputation
- Avoids overfitting by excluding weak predictors

### Why Optimal Distribution-Preserving Downsampling?
- Standard random splitting may not preserve distributional properties in small datasets
- OPDIS minimizes distributional divergence between full dataset and subset
- Ensures training set is representative of full population
- Particularly important for datasets with n < 200

## Reproducibility

All random processes use fixed seeds:
- **Python transformations:** Deterministic (no random processes)
- **R imputation:** `seed = 42`
- **R train/test split:** `Seed = 42`

Running the pipeline with the same input data will produce identical results.

## Code Availability

This code is freely available on GitHub at https://github.com/JornLotsch/pain-threshold-sex-differences. All code is provided under the Creative Commons Attribution 4.0 International (CC-BY 4.0) License, allowing free use, modification, and distribution with appropriate attribution. The pipeline is fully documented and reproducible; all scripts include inline comments explaining key steps. To reproduce the analysis:

1. Clone this repository
2. Install dependencies as described in the Installation section
3. Update paths in both scripts
4. Follow the Usage instructions

For detailed methodology and validation, refer to the accompanying manuscript (citation to be added upon publication).

## License

This project is licensed under the **Creative Commons Attribution 4.0 International License (CC-BY 4.0)**. See [LICENSE](LICENSE) file for details. This license applies to both code and data in this repository.

**What CC-BY 4.0 means:**
- ✅ You can use, modify, and distribute this work for any purpose
- ✅ You can use it for commercial or non-commercial work
- ✅ You must give appropriate credit and indicate changes made
- ✅ You must provide a link to the license

## Citation

If you use this code or preprocessing pipeline, please cite:

```
Lötsch, J., et al. (2025). Quantitative sensory testing dataset assessing sex differences 
and sensitization in 125 healthy volunteers. Nature Scientific Data, [volume], article [number].
https://github.com/JornLotsch/pain-threshold-sex-differences
```

(Full citation to be completed upon publication)

### Package Citations

When using output from this pipeline, also cite the key R packages:

**lodi:**
> Boss, J., & Rix, A. (2020). lodi: Limit of Detection Imputation for Single-Pollutant Models. R package version 0.9.2. https://CRAN.R-project.org/package=lodi

**opdisDownsampling:**
> Lötsch, J., & Ultsch, A. (2021). opdisDownsampling: Optimal Distribution Preserving Down-Sampling of Bio-Medical Data. R package version 1.5. https://CRAN.R-project.org/package=opdisDownsampling

## Contact

For questions about this dataset or preprocessing pipeline, please open an issue on GitHub:
https://github.com/JornLotsch/pain-threshold-sex-differences/issues

Or contact the corresponding author (to be specified upon publication).

## File Structure

```
.
├── README.md                                      # This file
├── Python/
│   ├── Read_explore_PainThresholds.py            # Stage 1: Transformation
│   └── transformation_example.py                 # Condensed code excerpt
├── R/
│   ├── preprocessing_pain_Mendeley_mf_dataset.R  # Stage 2: Preprocessing
│   ├── prepare_pain_Mendeley_mf_dataset.R        # Stage 3: Deposited files, figures, QC
│   ├── verify_reproduction_from_deposited_data.R # Regeneration from deposited files
│   ├── FindMostCorrelatedVariables2.R            # Helper: Correlation analysis
│   ├── imputation_example.R                      # Condensed code excerpt
│   └── training_test_splitting_example.R         # Condensed code excerpt
└── data/                                          # Not included in repository
    └── README.md                                  # Data availability information
```

## Data Availability

The complete participant-level pain threshold dataset and metadata are publicly available through Mendeley Data:

**Mendeley Data Repository:**
> Lotsch, Jorn; Flühr, Karin; Neddermeyer, Till J. (2026), "Quantitative sensory testing dataset with capsaicin and menthol sensitization assessing modality-specific sex differences and sensitization effect sizes in 125 healthy volunteers", Mendeley Data, V1, doi: 10.17632/y6v3rtpgps.1

### Deposited Data Files (Available on Mendeley Data)

The complete dataset includes:
- **exp_pain_data_orig.csv** - Untransformed participant-level pain threshold estimates in original units (125 subjects × 9 variables: ID + 8 thresholds). Trial-level observations are not included.
- **exp_pain_data_transformed.csv** - ln(x + 1) transformed, imputed and z-standardized pain thresholds and sensitization effects, ready for machine learning (125 subjects × 12 variables: ID + 11 features)
- **exp_pain_metadata.csv** - Participant metadata including age, sex, smoking status, experimenter, and train-validation assignment (125 subjects × 6 variables)

### Generated Output Files (From This Preprocessing Pipeline)

After running the preprocessing pipeline:
- `PainThresholds_scaled_Training.csv` - Final training dataset (80%, n = 100) with z-score standardized features
- `PainThresholds_scaled_Test.csv` - Final validation dataset (20%, n = 25) with z-score standardized features

### Data Not Included in This Repository

The source file `Phaeno_125.csv` (original single-file format, German variable names) is not included in this repository; use the Mendeley Data repository files instead. However, the original `Phaeno_125.csv` file is available from the corresponding author upon reasonable request.


---

**Last updated:** October 13, 2026 (revision)  
**Repository:** https://github.com/JornLotsch/pain-threshold-sex-differences  
**License:** CC-BY 4.0  
**DOI:** [To be assigned upon publication]
