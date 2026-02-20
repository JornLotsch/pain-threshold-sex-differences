################################################################################
# Pain Threshold Data Preprocessing Pipeline
################################################################################
#
# This script performs preprocessing of pain threshold data including:
# 1. Data loading and preparation
# 2. Imputation of censored data using correlated variables
# 3. Scaling and train/test split generation
#
# Input:  Phaeno_125.csv (raw phenotype data)
#         dfPainThresholdsAnalyzed_log_eff.csv (log-transformed pain thresholds)
# Output: PainThresholds_scaled_Training.csv (80% training set, scaled)
#         PainThresholds_scaled_Test.csv (20% test set, scaled)
#
################################################################################

############################### Paths #########################################
# Note: These paths are specific to the original analysis environment
# Users should modify these to match their local directory structure
pfad_o <- "/home/joern/Aktuell/PainGenesDrugs/"
pfad_u <- "09Originale/"
pfad_u1 <- "08AnalyseProgramme/"
pfad_u2 <- "04Umatrix/"

#################### Libraries ################################################
library(parallel)
library(lodi)                # For censored data imputation (clmi function)
library(dplyr)

# Source custom function for correlation analysis
source("FindMostCorrelatedVariables2.R")


nProc <- detectCores() - 1 # Number of cores for parallel processing

################################################################################
######################### STEP 1: Read and prepare data ########################
################################################################################
# Load phenotype data (125 subjects)
Phex2_125 <- read.csv(paste0(pfad_o, pfad_u, "Phaeno_125.csv"))
Phex2_125 <- data.frame(sapply(Phex2_125, function(x) as.numeric(as.character(x))))
names(Phex2_125)
rownames(Phex2_125) <- Phex2_125$Probandennummer

# Define pain threshold variable names
# Non-sensitized: baseline pain thresholds to different stimuli
PainThresholdVarnames_nonsensitized <- c("vonFrey", "Hitze", "Kaelte", "Druck", "Strom")
# Sensitized: pain thresholds after application of sensitizing agent
PainThresholdVarnames_sensitized <- c("vonFrey_C", "Hitze_C", "Kaelte_M")

# Load pre-transformed pain threshold data (log-transformed with effect calculations)
PainThresholdsData_transformed <- read.csv(paste0(pfad_o, "/08AnalyseProgramme/Python/", "dfPainThresholdsAnalyzed_log_eff.csv"), row.names = 1)
# Combine sex variable (0=male, 1=female) with pain threshold data
PainThresholdsData_transformed_all <- cbind.data.frame(Sex = Phex2_125$GeschlechtN, PainThresholdsData_transformed)
names(PainThresholdsData_transformed_all)
dim(PainThresholdsData_transformed)

################################################################################
################ STEP 2: Identify variables for imputation #####################
################################################################################
# Pain threshold measurements have censoring due to instrument limits:
# - Cold (Kaelte): Some subjects did not feel pain even at 0°C (coded as 0)
# - vonFrey: Some subjects did not feel pain even at max 300mN (coded as log(301))
#
# Strategy: Use correlated variables to impute censored values via clmi()
# (Censored data Left-truncated Multiple Imputation)
# Correlation threshold: |r| >= 0.4 was used to select predictor variables
################################################################################

# Find and rank correlations to identify predictor variables for imputation
pCorrRanked <- FindMostCorrelatedVariables2(Data = subset(PainThresholdsData_transformed_all, select = !names(PainThresholdsData_transformed_all) %in% "Sex"), TargetVariable = "Kaelte")

# Identify variables with correlation >= 0.4 to cold thresholds
pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4, ]
CorrKalt <- pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4 & (pCorrRanked$CorrTab$data$col_1 %in% c("Kaelte", "Kaelte_M") | pCorrRanked$CorrTab$data$col_2 %in% c("Kaelte", "Kaelte_M")), ]
# Exclude derived variables (CapsHeat, CapsvFrey, MenthCold) from predictors
# as they are difference scores and would create circular dependencies
CorrWkaelte <- setdiff(unique(c(CorrKalt$col_1, CorrKalt$col_2)), c("Kaelte", "Kaelte_M", "vonFrey", "vonFrey_C",  "CapsHeat",   "CapsvFrey",    "MenthCold"))

# Identify variables with correlation >= 0.4 to vonFrey thresholds
CorrFrey <- pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4 & (pCorrRanked$CorrTab$data$col_1 %in% c("vonFrey", "vonFrey_C") | pCorrRanked$CorrTab$data$col_2 %in% c("vonFrey", "vonFrey_C")), ]
CorrWFrey <- setdiff(unique(c(CorrFrey$col_1, CorrFrey$col_2)), c("Kaelte", "Kaelte_M", "vonFrey", "vonFrey_C",  "CapsHeat",   "CapsvFrey",    "MenthCold"))

################################################################################
#################### STEP 3: Impute censored data ##############################
################################################################################
# Prepare data for imputation
# Reference: https://cran.r-project.org/web/packages/lodi/vignettes/lodi.html
PainThresholdsData_transformed_all_toImpute <- PainThresholdsData_transformed_all

# Count censored observations
sum(PainThresholdsData_transformed_all_toImpute$Kaelte == 0)
sum(PainThresholdsData_transformed_all_toImpute$Kaelte_M == 0)
sum(PainThresholdsData_transformed_all_toImpute$vonFrey >= log(300+1))
sum(PainThresholdsData_transformed_all_toImpute$vonFrey_C >= log(300+1))

# Create NA versions of censored variables for imputation
# Kaelte: 0 indicates no pain at minimum temperature (left-censored)
PainThresholdsData_transformed_all_toImpute$Kaelte_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$Kaelte == 0, NA, PainThresholdsData_transformed_all_toImpute$Kaelte)
PainThresholdsData_transformed_all_toImpute$Kaelte_M_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$Kaelte_M == 0, NA, PainThresholdsData_transformed_all_toImpute$Kaelte_M)
# vonFrey: log(301) indicates no pain at maximum pressure (right-censored, so negate for left-censoring)
PainThresholdsData_transformed_all_toImpute$vonFrey_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$vonFrey >=  log(300+1), NA, -PainThresholdsData_transformed_all_toImpute$vonFrey)
PainThresholdsData_transformed_all_toImpute$vonFrey_C_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$vonFrey_C >=  log(300+1), NA, -PainThresholdsData_transformed_all_toImpute$vonFrey_C)
PainThresholdsData_transformed_all_toImpute
# Set limits of detection (LOD) for imputation algorithm
PainThresholdsData_transformed_all_toImpute$lod <- 0           # For cold (lower limit)
PainThresholdsData_transformed_all_toImpute$lodF <- -log(300+1) # For vonFrey (upper limit, negated)


# Perform censored left multiple imputation (clmi) for each censored variable
# n.imps=5 generates 5 imputed datasets; median is taken across imputations

# Impute Kaelte (cold threshold) using correlated predictors
clmi.out <- clmi(
  formula = as.formula(paste("Kaelte_NA ~ ", paste(CorrWkaelte, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lod, seed = 42, n.imps = 5
)
Kaelte_imputed <- apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$Kaelte_NA_transform_imputed)), 1, median)

# Impute Kaelte_M (cold + menthol threshold)
clmi.out <- clmi(
  formula = as.formula(paste("Kaelte_M_NA ~ ", paste(CorrWkaelte, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lod, seed = 42, n.imps = 5
)
Kaelte_M_imputed <- apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$Kaelte_M_NA_transform_imputed)), 1, median)

# Impute vonFrey (punctate pressure threshold)
clmi.out <- clmi(
  formula = as.formula(paste("vonFrey_NA ~ ", paste(CorrWFrey, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lodF, seed = 42, n.imps = 5
)
vonFrey_imputed <- -apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$vonFrey_NA_transform_imputed)), 1, median)

# Impute vonFrey_C (vonFrey + capsaicin threshold)
clmi.out <- clmi(
  formula = as.formula(paste("vonFrey_C_NA ~ ", paste(CorrWFrey, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lodF, seed = 42, n.imps = 5
)
vonFrey_C_imputed <- -apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$vonFrey_C_NA_transform_imputed)), 1, median)

################################################################################
############# STEP 4: Construct final imputed dataset #########################
################################################################################
# Remove original censored variables and replace with imputed versions
# Also recalculate derived variables (difference scores) using imputed values

# Remove original censored variables and derived variables
PainThresholdsData_transformed_imputed <- PainThresholdsData_transformed %>% select(-c(Kaelte, Kaelte_M, MenthCold, vonFrey, vonFrey_C, CapsvFrey))

# Add imputed variables with "_I" suffix
PainThresholdsData_transformed_imputed$Kaelte_I <- Kaelte_imputed
PainThresholdsData_transformed_imputed$Kaelte_M_I <- Kaelte_M_imputed
# MenthCold_I: effect of menthol on cold threshold (difference score)
PainThresholdsData_transformed_imputed$MenthCold_I <- Kaelte_imputed - Kaelte_M_imputed
PainThresholdsData_transformed_imputed$vonFrey_I <- vonFrey_imputed
PainThresholdsData_transformed_imputed$vonFrey_C_I <- vonFrey_C_imputed
# CapsvFrey_I: effect of capsaicin on vonFrey threshold (difference score)
PainThresholdsData_transformed_imputed$CapsvFrey_I <- vonFrey_imputed - vonFrey_C_imputed

################################################################################
########### STEP 5: Add target variable (sex) for supervised learning #########
################################################################################
PainThresholdsData_transformed_imputed_proj <- cbind.data.frame(Target = PainThresholdsData_transformed_all$Sex, PainThresholdsData_transformed_imputed)
PainThresholds <- PainThresholdsData_transformed_imputed_proj

################################################################################
########### STEP 6: Scale and generate train/test split (FINAL OUTPUT) ########
################################################################################
# 1. Z-score standardization: Mean=0, SD=1 for all features
# 2. Optimal distribution downsampling: Creates balanced 80/20 train/test split
#    that preserves distributional properties of the full dataset

# Z-score standardization of all features (excluding Target column)
PainThresholds_scaled <- PainThresholds
PainThresholds_scaled[,2:ncol(PainThresholds_scaled)] <- scale(PainThresholds_scaled[,2:ncol(PainThresholds_scaled)])

# opdisDownsampling: Optimal distribution-preserving downsampling
# Creates train/test split that maintains distributional characteristics
# Size = 0.8 means 80% training, 20% test
# nTrials = 1000000 ensures thorough search for optimal split
# Seed = 42 for reproducibility
PainThresholdsDown_scaled <- opdisDownsampling::opdisDownsampling(Data = within(PainThresholds_scaled, rm(Target)),
                                                           Cls = PainThresholds_scaled$Target, Size = 0.8*nrow(PainThresholds_scaled), Seed = 42,
                                                           nTrials = 1000000, MaxCores = 30, PCAimportance = F)

################################################################################
########################### WRITE FINAL OUTPUT FILES ###########################
################################################################################
write.csv(PainThresholdsDown_scaled$ReducedData, "PainThresholds_scaled_Training.csv")
write.csv(PainThresholdsDown_scaled$RemovedData, "PainThresholds_scaled_Test.csv")
################################################################################
