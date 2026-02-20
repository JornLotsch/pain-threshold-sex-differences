# Quantitative Sensory Testing Dataset - Preprocessing Pipeline

## Dataset: Quantitative Sensory Testing Dataset with Capsaicin and Menthol Sensitization Assessing Modality-Specific Sex Differences and Sensitization Effect Sizes in 125 Healthy Volunteers

## Quantitative sensory testing dataset assessing sex differences and sensitization in 125 healthy volunteers

This repository contains the complete preprocessing pipeline for the quantitative sensory testing (QST) dataset accompanying the publication:

> **Quantitative sensory testing dataset assessing sex differences and sensitization in 125 healthy volunteers**
>
> [Citation to be added upon publication]
>
> **Repository:** https://github.com/JornLotsch/pain-threshold-sex-differences

## Overview

This repository contains the preprocessing pipeline to prepare the quantitative sensory testing dataset for machine learning applications. The pipeline implements two stages: (1) signed logarithmic transformation and sensitization effect calculation in Python, and (2) censored data imputation, standardization, and train-test splitting in R. The result is a fully processed dataset with standardized features ready for statistical modeling or machine learning analysis.

## Dataset Description

This dataset contains experimental pain threshold measurements collected from 125 healthy adult volunteers under both baseline and chemically-sensitized conditions. The data were generated using established quantitative sensory testing (QST) protocols, including mechanical (punctate and blunt pressure), thermal (heat and cold), and electrical pain modalities. Sensitization was induced using topical capsaicin (for heat and punctate pressure) and menthol (for cold) to assess pain plasticity and sensitization mechanisms.

**Dataset Citation:**
> Lotsch, Jorn; Flühr, Karin; Neddermeyer, Till J. (2026), "Quantitative sensory testing dataset with capsaicin and menthol sensitization assessing modality-specific sex differences and sensitization effect sizes in 125 healthy volunteers", Mendeley Data, V1, doi: 10.17632/y6v3rtpgps.1

### Input Data

**Input:** Raw pain threshold measurements from 125 participants
- 5 non-sensitized pain modalities: von Frey (punctate pressure), heat, cold, blunt pressure, electrical stimulation
- 3 sensitized measurements: von Frey + capsaicin, heat + capsaicin, cold + menthol
- Metadata: participant ID, age, sex, smoking status, experimenter identifier

### Output Data

**Output:** Two files ready for machine learning
- `PainThresholds_scaled_Training.csv` - Training set (80% of data, n≈100)
- `PainThresholds_scaled_Test.csv` - Validation set (20% of data, n≈25)

Each output file contains:
- 11 standardized pain threshold features (8 transformed thresholds + 3 sensitization effects)
- Target variable (sex: 0=male, 1=female)

## Pipeline Structure

### Stage 1: Data Transformation (Python)

**Script:** `Python/Read_explore_PainThresholds.py`

**Purpose:** Apply signed logarithmic transformation and calculate sensitization effects

**Steps:**
1. Read raw pain threshold data (`Phaeno_125.csv`)
2. Apply signed logarithmic transformation: `sign(x) × log(|x| + 1)`
   - Handles zero and negative values
   - Reduces skewness while preserving direction
3. Calculate sensitization effects (difference scores):
   - `CapsHeat` = Heat - Heat_Capsaicin
   - `CapsvFrey` = vonFrey - vonFrey_Capsaicin
   - `MenthCold` = Cold - Cold_Menthol
4. Write transformed data: `dfPainThresholdsAnalyzed_log_eff.csv`

**Key transformations:**
- Signed log transformation normalizes skewed distributions
- Sensitization effects quantify how chemical sensitizers alter pain perception

### Stage 2: Data Preprocessing (R)

**Script:** `R/preprocessing_pain_Mendeley_mf_dataset.R`

**Purpose:** Handle censored data, standardize features, and create train/test split

**Steps:**

1. **Data Loading**
   - Load transformed pain thresholds from Stage 1
   - Combine with sex variable from raw phenotype data

2. **Censored Data Imputation**
   - **Problem:** Some measurements hit instrument limits
     - Cold thresholds: No pain at 0°C (left-censored, coded as 0)
     - von Frey thresholds: No pain at 300mN max force (right-censored, coded as log(301))
   - **Solution:** Conditional Linear Multiple Imputation (CLMI)
     - Identify predictor variables with |Spearman r| ≥ 0.4
     - Generate 5 imputed datasets per censored variable
     - Take median across imputations
     - Recalculate sensitization effects using imputed values
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

### Helper Functions

**Script:** `R/FindMostCorrelatedVariables2.R`

**Purpose:** Identify highly correlated variables for imputation predictor selection

**Features:**
- Computes all pairwise Spearman correlations
- Generates visualization plots with confidence intervals
- Highlights correlations involving target variables
- Used to select predictors (|r| ≥ 0.4) for censored data imputation

## Software Requirements and Installation

### Python
- **Version:** 3.12.3
- **Packages:** `numpy` (≥1.26.4), `pandas`
- **Install:** `pip install numpy pandas`

### R
- **Version:** 4.5.2
- **Required packages:** `lodi`, `opdisDownsampling`, `dplyr`, `parallel`
- **Optional packages:** `inspectdf`, `ggplot2` (for correlation visualization in helper script)
- **Install:**
```r
install.packages(c("lodi", "opdisDownsampling", "dplyr"))
# Optional for correlation analysis:
install.packages(c("inspectdf", "ggplot2"))
```

## Usage

### Quick Start

Before running the pipeline, ensure you have:
1. **Python 3.12.3** and **R 4.5.2** installed
2. **All required packages** installed (see Installation section)
3. **Raw data file** `Phaeno_125.csv` placed in the appropriate directory
4. **Updated paths** in both scripts (see Path Configuration below)

### Complete Pipeline

**Note:** You must have the raw data file `Phaeno_125.csv` to run the complete pipeline.

1. **Run Python transformation:**
   ```bash
   cd Python/
   python3 Read_explore_PainThresholds.py
   ```
   Output: `dfPainThresholdsAnalyzed_log_eff.csv`

2. **Run R preprocessing:**
   ```r
   setwd("R/")
   source("preprocessing_pain_Mendeley_mf_dataset.R")
   ```
   Output:
   - `PainThresholds_scaled_Training.csv`
   - `PainThresholds_scaled_Test.csv`

### Path Configuration

**Important:** Update the paths in both scripts to match your local directory structure:

**Python script** (`Python/Read_explore_PainThresholds.py`):
```python
pfad_o = "/path/to/your/data/"  # Absolute path to parent directory containing raw data
pfad_u1 = "09Originale/"         # Relative path to subfolder within pfad_o
```

**R script** (`R/preprocessing_pain_Mendeley_mf_dataset.R`):
```r
pfad_o  <- "/path/to/your/data/"  # Absolute path to parent directory (same as Python)
pfad_u  <- "09Originale/"          # Raw data folder (relative to pfad_o)
pfad_u1 <- "08AnalyseProgramme/"   # Analysis folder (relative to pfad_o)
```

**Expected directory structure:**
```
/path/to/your/data/
├── 09Originale/
│   ├── Phaeno_125.csv              # Raw pain threshold data (required)
│   └── [other phenotype files]
└── 08AnalyseProgramme/             # This repository
    ├── Python/
    ├── R/
    └── README.md
```

## Key Methodological Decisions

### Why Signed Logarithmic Transformation?
- Pain thresholds are highly right-skewed
- Log transformation normalizes distributions
- Signed version preserves information about direction (important for sensitization effects)

### Why Conditional Linear Multiple Imputation?
- Censored data cannot simply be discarded (biases results)
- CLMI uses correlated variables to predict censored values
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
> Boss, J., Wetzels, M., & Muggeo, V. M. R. (2020). lodi: Limit of Detection Imputation for Single-Pollutant Models. R package version 0.9.2. https://CRAN.R-project.org/package=lodi

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
│   └── Read_explore_PainThresholds.py            # Stage 1: Transformation
├── R/
│   ├── preprocessing_pain_Mendeley_mf_dataset.R  # Stage 2: Preprocessing
│   └── FindMostCorrelatedVariables2.R            # Helper: Correlation analysis
└── data/                                          # Not included in repository
    └── README.md                                  # Data availability information
```

## Data Availability

The complete raw pain threshold dataset and metadata are publicly available through Mendeley Data:

**Mendeley Data Repository:**
> Lotsch, Jorn; Flühr, Karin; Neddermeyer, Till J. (2026), "Quantitative sensory testing dataset with capsaicin and menthol sensitization assessing modality-specific sex differences and sensitization effect sizes in 125 healthy volunteers", Mendeley Data, V1, doi: 10.17632/y6v3rtpgps.1

### Raw Data Files (Available on Mendeley Data)

The complete dataset includes:
- **exp_pain_data_orig.csv** - Original pain threshold measurements (125 subjects × 9 variables)
- **exp_pain_data_transformed.csv** - Transformed pain thresholds ready for machine learning (125 subjects × 12 variables)
- **exp_pain_metadata.csv** - Participant metadata including age, sex, smoking status, experimenter, and train-validation assignment (125 subjects × 6 variables)

### Generated Output Files (From This Preprocessing Pipeline)

After running the preprocessing pipeline:
- `PainThresholds_scaled_Training.csv` - Final training dataset (80%, n≈100) with z-score standardized features
- `PainThresholds_scaled_Test.csv` - Final validation dataset (20%, n≈25) with z-score standardized features

### Data Not Included in This Repository

The raw data file `Phaeno_125.csv` (original single-file format) is not included; use the Mendeley Data repository files instead.


---

**Last updated:** February 20, 2025  
**Repository:** https://github.com/JornLotsch/pain-threshold-sex-differences  
**License:** CC-BY 4.0  
**DOI:** [To be assigned upon publication]
