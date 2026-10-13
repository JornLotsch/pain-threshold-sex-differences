################################################################################
# Verify that the deposited transformed data and the training/validation
# assignment can be regenerated directly from the deposited Mendeley files
################################################################################
#
# Input (Mendeley Data, doi: 10.17632/y6v3rtpgps.1):
#   exp_pain_data_orig.csv        - participant-level threshold estimates
#   exp_pain_data_transformed.csv - deposited transformed data (to be verified)
#   exp_pain_metadata.csv         - metadata incl. Sex and TrainingValidation
#
# The script re-implements, in R only, the steps of
#   Python/Read_explore_PainThresholds.py     (ln(x + 1), sensitization effects)
#   R/preprocessing_pain_Mendeley_mf_dataset.R (imputation, scaling, split)
# using identical settings (predictor selection |rho| >= 0.4, clmi with
# n.imps = 5 and seed = 42, opdisDownsampling with Seed = 42, nTrials = 1e6),
# and compares the result with the deposited files.
#
# Order of steps (as in the original pipeline): transformation -> imputation on
# all 125 subjects -> z-scaling on all 125 subjects -> training/validation split.
#
################################################################################

############### Paths ##############################
# Users should set this to the directory containing the downloaded Mendeley files
DATA_DIR <- "/home/joern/Aktuell/QST_Mendeley_MF/DataSetPublished"

############### Libraries ##############################
library(lodi)
library(dplyr)
library(ggplot2)             # Required by FindMostCorrelatedVariables2()
library(opdisDownsampling)

source("FindMostCorrelatedVariables2.R")

N_TRIALS <- 1000000          # As in the original analysis
MAX_CORES <- 30

############### Read deposited data ##############################
orig_data        <- read.csv(file.path(DATA_DIR, "exp_pain_data_orig.csv"))
transformed_data <- read.csv(file.path(DATA_DIR, "exp_pain_data_transformed.csv"))
metadata         <- read.csv(file.path(DATA_DIR, "exp_pain_metadata.csv"))

stopifnot(identical(orig_data$ID, metadata$ID),
          identical(as.character(orig_data$ID), as.character(transformed_data$ID)))

############### Step 1: ln(x + 1) transformation and sensitization effects ##############
# Variable names of the original analysis are used, as the predictor selection
# depends on them only through their identity, not their order
pain_data_log <- data.frame(
  vonFrey   = orig_data$vonFrey,
  Hitze     = orig_data$Heat,
  Kaelte    = orig_data$Cold,
  Druck     = orig_data$Pressure,
  Strom     = orig_data$Current,
  vonFrey_C = orig_data$vonFrey_Capsaicin,
  Hitze_C   = orig_data$Heat_Capsaicin,
  Kaelte_M  = orig_data$Cold_Menthol
)
# Natural logarithm, as np.log in the Python stage
pain_data_log <- as.data.frame(lapply(pain_data_log, function(x) sign(x) * log(abs(x + 1))))
pain_data_log$CapsHeat  <- pain_data_log$Hitze - pain_data_log$Hitze_C
pain_data_log$CapsvFrey <- pain_data_log$vonFrey - pain_data_log$vonFrey_C
pain_data_log$MenthCold <- pain_data_log$Kaelte - pain_data_log$Kaelte_M
rownames(pain_data_log) <- orig_data$ID

############### Step 2: Predictor selection for imputation ##############
pCorrRanked <- FindMostCorrelatedVariables2(Data = pain_data_log, TargetVariable = "Kaelte")
corr_tab <- pCorrRanked$CorrTab$data
excluded <- c("Kaelte", "Kaelte_M", "vonFrey", "vonFrey_C", "CapsHeat", "CapsvFrey", "MenthCold")

select_predictors <- function(targets) {
  x <- corr_tab[abs(corr_tab$corr) >= 0.4 & (corr_tab$col_1 %in% targets | corr_tab$col_2 %in% targets), ]
  setdiff(unique(c(x$col_1, x$col_2)), excluded)
}
CorrWkaelte <- select_predictors(c("Kaelte", "Kaelte_M"))
CorrWFrey   <- select_predictors(c("vonFrey", "vonFrey_C"))
cat("Predictors for cold:", CorrWkaelte, "\n")
cat("Predictors for von Frey:", CorrWFrey, "\n")

############### Step 3: Imputation of censored values (all 125 subjects) ##############
df_to_impute <- pain_data_log
df_to_impute$Kaelte_NA    <- ifelse(df_to_impute$Kaelte == 0, NA, df_to_impute$Kaelte)
df_to_impute$Kaelte_M_NA  <- ifelse(df_to_impute$Kaelte_M == 0, NA, df_to_impute$Kaelte_M)
# Sign reversal for clmi() only (right- to left-censoring); reversed after imputation
df_to_impute$vonFrey_NA   <- ifelse(df_to_impute$vonFrey >= log(300 + 1), NA, -df_to_impute$vonFrey)
df_to_impute$vonFrey_C_NA <- ifelse(df_to_impute$vonFrey_C >= log(300 + 1), NA, -df_to_impute$vonFrey_C)
df_to_impute$lod  <- 0
df_to_impute$lodF <- -log(300 + 1)

impute_censored <- function(variable, predictors, lod_column) {
  # clmi() expects the lod column name unquoted, hence the constructed call
  clmi.out <- eval(bquote(clmi(
    formula = .(as.formula(paste(variable, "~", paste(predictors, collapse = " + ")))),
    df = df_to_impute, lod = .(as.name(lod_column)), seed = 42, n.imps = 5
  )))
  apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x[[paste0(variable, "_transform_imputed")]])), 1, median)
}

Kaelte_imputed    <- impute_censored("Kaelte_NA", CorrWkaelte, "lod")
Kaelte_M_imputed  <- impute_censored("Kaelte_M_NA", CorrWkaelte, "lod")
vonFrey_imputed   <- -impute_censored("vonFrey_NA", CorrWFrey, "lodF")
vonFrey_C_imputed <- -impute_censored("vonFrey_C_NA", CorrWFrey, "lodF")

############### Step 4: Assemble and z-scale (all 125 subjects) ##############
regenerated <- data.frame(
  Heat                     = pain_data_log$Hitze,
  Pressure                 = pain_data_log$Druck,
  Current                  = pain_data_log$Strom,
  Heat_Capsaicin           = pain_data_log$Hitze_C,
  Capsaicin_Effect_Heat    = pain_data_log$CapsHeat,
  Cold                     = Kaelte_imputed,
  Cold_Menthol             = Kaelte_M_imputed,
  Menthol_Effect_Cold      = Kaelte_imputed - Kaelte_M_imputed,
  vonFrey                  = vonFrey_imputed,
  vonFrey_Capsaicin        = vonFrey_C_imputed,
  Capsaicin_Effect_vonFrey = vonFrey_imputed - vonFrey_C_imputed
)
regenerated <- as.data.frame(scale(regenerated))
rownames(regenerated) <- orig_data$ID

############### Verification 1: transformed data ##############
max_abs_diff <- sapply(names(regenerated), function(v) max(abs(regenerated[[v]] - transformed_data[[v]])))
cat("Maximum absolute difference to exp_pain_data_transformed.csv per variable:\n")
print(max_abs_diff)
cat("Transformed data reproduced:",
    isTRUE(all.equal(regenerated, transformed_data[, names(regenerated)], check.attributes = FALSE)), "\n")

############### Step 5 / Verification 2: training/validation split ##############
# Target: sex (0 = male, 1 = female)
split <- opdisDownsampling::opdisDownsampling(
  Data = regenerated, Cls = ifelse(metadata$Sex == "Female", 1, 0),
  Size = 0.8 * nrow(regenerated), Seed = 42,
  nTrials = N_TRIALS, MaxCores = MAX_CORES, PCAimportance = FALSE
)
ids_training_regenerated <- rownames(split$ReducedData)
ids_training_deposited   <- as.character(metadata$ID[metadata$TrainingValidation == "Training"])
cat("Training/validation assignment reproduced:",
    setequal(ids_training_regenerated, ids_training_deposited), "\n")

sessionInfo()
