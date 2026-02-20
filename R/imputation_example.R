# Conditional linear multiple imputation (clmi) for censored measurements
# Left-censored: Cold thresholds = 0 (below detection limit)
# Right-censored: von Frey thresholds ≥ log(301) (above 300 mN maximum force)

# 1. Identify censored values
PainThresholdsData_transformed_all_toImpute$Kaelte_NA <-
  ifelse(Kaelte == 0, NA, Kaelte)
PainThresholdsData_transformed_all_toImpute$vonFrey_NA <-
  ifelse(vonFrey >= log(301), NA, -vonFrey)

# 2. Select predictors: variables with |Spearman correlation| ≥ 0.4
CorrKalt <- CorrTab[abs(corr) >= 0.4 &
                      (col_1 %in% c("Kaelte","Kaelte_M") |
                        col_2 %in% c("Kaelte","Kaelte_M")), ]
CorrWkaelte <- setdiff(unique(c(CorrKalt$col_1, CorrKalt$col_2)),
                       c("Kaelte","Kaelte_M","vonFrey","vonFrey_C",
                         "CapsHeat","CapsvFrey","MenthCold"))

# 3. Impute with 5 iterations, extract median
clmi.out <- clmi(
  formula = as.formula(paste("Kaelte_NA ~", paste(CorrWkaelte, collapse=" + "))),
  df = PainThresholdsData_transformed_all_toImpute,
  lod = 0,  # limit of detection
  seed = 42,
  n.imps = 5
)
Kaelte_imputed <- apply(
  data.frame(lapply(clmi.out$imputed.dfs,
                    function(x) x$Kaelte_NA_transform_imputed)),
  1, median
)

# 4. Recalculate sensitization effects with imputed values
MenthCold_I <- Kaelte_I - Kaelte_M_I
CapsvFrey_I <- vonFrey_I - vonFrey_C_I