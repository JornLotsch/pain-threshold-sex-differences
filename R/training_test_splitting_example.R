# Optimal distribution-preserving downsampling for 80/20 train/test split

# 1. Add binary sex classification target (0=male, 1=female)
PainThresholds <- cbind.data.frame(
  Target = Sex,
  PainThresholdsData_transformed_imputed
)

# 2. Z-score standardization of all 11 pain threshold variables
PainThresholds_scaled <- PainThresholds
PainThresholds_scaled[, 2:ncol(PainThresholds_scaled)] <-
  scale(PainThresholds_scaled[, 2:ncol(PainThresholds_scaled)])

# 3. Iterative optimization to identify 80% subset that minimizes
#    distributional divergence from full dataset
PainThresholdsDown_scaled <- opdisDownsampling(
  Data = within(PainThresholds_scaled, rm(Target)),
  Cls = PainThresholds_scaled$Target,
  Size = 0.8 * nrow(PainThresholds_scaled),  # 80% training
  Seed = 42,
  nTrials = 1000000,  # optimization iterations
  MaxCores = 30,
  PCAimportance = FALSE
)

# Output files (n_train ≈ 100, n_test ≈ 25):
# - PainThresholds_scaled_Training.csv: optimally selected 80%
# - PainThresholds_scaled_Test.csv: remaining 20% for validation