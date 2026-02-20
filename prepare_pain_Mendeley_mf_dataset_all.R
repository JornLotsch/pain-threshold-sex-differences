
############### Libraries ##############################

############### Constants ##############################
# File paths
INPUT_FILE_PATH <- "/home/joern/Dokumente/PainGenesDrugs/09Originale/Phaeno_125.csv"
OUTPUT_DIR <- "/home/joern/Aktuell/QST_Mendeley_MF/DataSetPublished"

# Variable definitions
PAIN_TESTS_NAMES <- c(
  PainThresholdVarnames_nonsensitized = c("vonFrey", "Hitze", "Kaelte", "Druck", "Strom"),
  PainThresholdVarnames_sensitized = c("vonFrey_C", "Hitze_C", "Kaelte_M")
)

# Variable definitions
PAIN_TESTS_NAMES_ENGLISH <- c(
  PainThresholdVarnames_nonsensitized = c("vonFrey", "Heat", "Cold", "Pressure", "Current"),
  PainThresholdVarnames_sensitized = c("vonFrey_Capsaicin", "Heat_Capsaicin", "Cold_Menthol")
)


METADATA <- c("Probandennummer", "Alter", "Raucher", "Geschlecht", "Pruefer")
METADATA_ENGLISH <-c("ID", "Age", "Smoker", "Sex", "Tester")


rename_pain_data_columns <- function(data, new_names) {
  # Assign new column names for consistency.
  names(data) <- new_names
  data
}



############### Load Data ##############################

Phaeno_125 <- read.csv(INPUT_FILE_PATH, )
Phaeno_125_orig_pain_data <- Phaeno_125[, c("Probandennummer", PAIN_TESTS_NAMES)]
Phaeno_125_orig_pain_data_renamed  <- rename_pain_data_columns(Phaeno_125_orig_pain_data, c("ID", PAIN_TESTS_NAMES_ENGLISH))

Phaeno_125_metadata <- Phaeno_125[,METADATA]
Phaeno_125_metadata_renamed <- rename_pain_data_columns(Phaeno_125_metadata, METADATA_ENGLISH)
Phaeno_125_metadata_renamed$Tester <- as.integer(as.factor(Phaeno_125_metadata_renamed$Tester))
Phaeno_125_metadata_renamed$Smoker[Phaeno_125_metadata_renamed$Smoker == "J"] <- "Y"
Phaeno_125_metadata_renamed$Sex[Phaeno_125_metadata_renamed$Sex == "W"] <- "F"

# Dataset paths (now only training + validation CSVs)
base_path <- "/home/joern/Dokumente/PainGenesDrugs"
r_path <- "08AnalyseProgramme/R"

train_file <- file.path(base_path, r_path, "PainThresholds_scaled_Training.csv")
val_file <- file.path(base_path, r_path, "PainThresholds_scaled_Test.csv")

# Dataset-specific configuration
DATASET_COLUMN_NAMES <- c(
  "Heat", "Pressure", "Current", "Heat_Capsaicin",
  "Capsaicin_Effect_Heat", "Cold", "Cold_Menthol", "Menthol_Effect_Cold",
  "vonFrey", "vonFrey_Capsaicin", "Capsaicin_Effect_vonFrey"
)

CURATED_COLUMN_NAMES <- c(
  "Heat", "Pressure", "Current", "Heat_Capsaicin",
  "Cold", "Cold_Menthol",
  "vonFrey", "vonFrey_Capsaicin"
)

COLUMNS_COLINEAR <- NULL # Will be determined later

# Noise addition parameters
ADD_NOISE_COLUMN <- TRUE
NOISE_COLUMN_SOURCE <- "Pressure"
NOISE_COLUMN_NAME <- "Pressure2"


###############################################################################
# Check data files existence
###############################################################################
if (!file.exists(train_file)) {
  stop(paste("Training data file not found:", train_file))
}
if (!file.exists(val_file)) {
  stop(paste("Validation data file not found:", val_file))
}


###############################################################################
# Load data and modify files when needed
###############################################################################

rename_pain_data_columns <- function(data, new_names) {
  # Assign new column names for consistency.
  names(data) <- new_names
  data
}


# --- Define paths ---
base_path <- "/home/joern/Dokumente/PainGenesDrugs"
r_path <- "08AnalyseProgramme/R"

train_file <- file.path(base_path, r_path, "PainThresholds_scaled_Training.csv")
val_file <- file.path(base_path, r_path, "PainThresholds_scaled_Test.csv")

# --- Load datasets ---
train_df <- read.csv(train_file, row.names = 1)
val_df <- read.csv(val_file, row.names = 1)

# --- Sizes to later reseparate ---
n_train <- nrow(train_df)
n_val <- nrow(val_df)

# --- Extract features + targets ---
train_features <- train_df[, - ncol(train_df)]
train_target <- train_df[, ncol(train_df)]

val_features <- val_df[, - ncol(val_df)]
val_target <- val_df[, ncol(val_df)]

# --- Combine into pain_data + target_data ---
pain_data <- rbind(train_features, val_features)
target_data <- c(train_target, val_target)

Phaeno_125_metadata_renamed$TrainingValidation <- ifelse(Phaeno_125_metadata_renamed$ID %in% as.integer(rownames(train_features)), "Training", "Validation")

# --- Rename columns once, for the combined dataset ---
pain_data <- rename_pain_data_columns(pain_data, DATASET_COLUMN_NAMES)

# --- Reseparate into train/validation, consistent with original splits ---
training_data_original <- pain_data[1:n_train,]
validation_data_original <- pain_data[(n_train + 1):(n_train + n_val),]

training_target <- target_data[1:n_train]
validation_target <- target_data[(n_train + 1):(n_train + n_val)]


# Save data file for use in other analyses
Pheno_125_prepared_data <- cbind.data.frame(Target = target_data, pain_data)
Pheno_125_prepared_data_1 <- cbind.data.frame(ID = rownames(Pheno_125_prepared_data), Pheno_125_prepared_data)

Pheno_125_prepared_data_1_sorted <- Pheno_125_prepared_data_1[match(Phaeno_125_metadata_renamed$ID, Pheno_125_prepared_data_1$ID), ]
plot(Pheno_125_prepared_data_1_sorted$Target ~ as.integer(as.factor(Phaeno_125_metadata_renamed$Sex))-1)

Pheno_125_prepared_data_1_sorted_2 <- Pheno_125_prepared_data_1_sorted[,!names(Pheno_125_prepared_data_1_sorted) %in% c("Target")]

# Write final data
# Original data
write.csv(Phaeno_125_orig_pain_data_renamed, file = paste0(OUTPUT_DIR, "/", "exp_pain_data_orig.csv"))
# Transfored data
write.csv(Pheno_125_prepared_data_1_sorted_2, file = paste0(OUTPUT_DIR, "/", "exp_pain_data_transformed.csv"))
# Metadata
write.csv(Phaeno_125_metadata_renamed, file = paste0(OUTPUT_DIR, "/", "exp_pain_metadata.csv"))

########################################################################################################

############################### Paths #############################################################
pfad_o <- "/home/joern/Aktuell/PainGenesDrugs/"
pfad_u <- "09Originale/"
pfad_u1 <- "08AnalyseProgramme/"
pfad_u2 <- "04Umatrix/"

#################### Libraries ############################################################
library(FactoMineR)
library(readxl)
library(opGMMassessment)
library(parallel)
library(pbmcapply)
library(cowplot)
library(ggplot2)
library(ggExtra)
library(GGally)
library(ggthemes)
library(dplyr)
library(cowplot)
library(effsize)
library(stringr)

source("/home/joern/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/ProjectionsBiomed_MainFunctions.R")
source("/home/joern/Aktuell/RFABCFeatureSelection2/08AnalyseProgramme/R/FindMostCorrelatedVariables2.R")


nProc <- detectCores() - 1 # 20#

######################### Read data ######################################################
Phex2_125 <- read.csv(paste0(pfad_o, pfad_u, "Phaeno_125.csv"))
Phex2_125 <- data.frame(sapply(Phex2_125, function(x) as.numeric(as.character(x))))
names(Phex2_125)
rownames(Phex2_125) <- Phex2_125$Probandennummer
PainThresholdVarnames_nonsensitized <- c("vonFrey", "Hitze", "Kaelte", "Druck", "Strom")
PainThresholdVarnames_sensitized <- c("vonFrey_C", "Hitze_C", "Kaelte_M")

PainThresholdsData_transformed <- read.csv(paste0(pfad_o, "/08AnalyseProgramme/Python/", "dfPainThresholdsAnalyzed_log_eff.csv"), row.names = 1)
PainThresholdsData_transformed_all <- cbind.data.frame(Sex = Phex2_125$GeschlechtN, PainThresholdsData_transformed)
names(PainThresholdsData_transformed_all)
dim(PainThresholdsData_transformed)

######################## Check correlation #################################################

corrThresholds <- PainThresholdsData_transformed_all %>%
  select(-Sex) %>%
  ggpairs(.,
          mapping = ggplot2::aes(
            colour = factor(PainThresholdsData_transformed_all$Sex),
            # fill=interaction(OlfDiag,Gender),
            alpha = .3
          ),
          lower = list(continuous = wrap("smooth", size = 1)),
          upper = list(continuous = wrap("cor", method = "spearman")),
          # diag = list(continuous = my_dens)
  ) +
  theme_linedraw() +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(color = "black", fill = "cornsilk")
  ) +
  scale_color_manual(values = colorblind_pal()(6)) +
  scale_fill_manual(values = colorblind_pal()(6))


pCorrMat <- plotCorrMatrix(X = within(PainThresholdsData_transformed_all, rm(Sex)))
pCorrRanked <- FindMostCorrelatedVariables2(Data = subset(PainThresholdsData_transformed_all, select = !names(PainThresholdsData_transformed_all) %in% "Sex"), TargetVariable = "Kaelte")
pClasswise <- plotClasswiseData(
  X = subset(PainThresholdsData_transformed_all, select = !names(PainThresholdsData_transformed_all) %in% "Sex"),
  y = PainThresholdsData_transformed_all$Sex, rows = 1
)

CohenD <- data.frame(EffectSize = apply(PainThresholdsData_transformed_all[2:ncol(PainThresholdsData_transformed_all)], 2,
                                        function(x) cohen.d(x ~ PainThresholdsData_transformed_all$Sex )$estimate))
CohenD$Stimulus = rownames(CohenD)
CohenD$color <- ifelse(CohenD$EffectSize >= 0, "F>M", "M>F")
for (i in grep("Kaelte", rownames(CohenD))) {
  if (CohenD$color[i] == "F>M") CohenD$color[i] <- "M>F"
  if (CohenD$color[i] == "M>F") CohenD$color[i] <- "F>M"
}

CohenD$Number <- 1:nrow(CohenD)

pCohenD <-
  ggplot(CohenD, aes(x = Stimulus, y = EffectSize, color = color, fill = color)) +
  geom_bar(stat = "identity", alpha = .5) +
  theme_light() +
  facet_wrap(Number~., scales = "free_x", nrow = 1) +
  theme(legend.position = c(.8, 1.1), strip.background = element_blank(), strip.text = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.7)), legend.direction="horizontal",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values = c("chartreuse3", "dodgerblue") )+
  scale_fill_manual(values = c("chartreuse3", "dodgerblue")) +
  guides(color = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "salmon") +
  labs(title = "Direction and effect size of gender difference", y = "Cohen's d", x = NULL, fill = "Pain\nsensitivity")

plot_grid(plot_grid(pClasswise + labs(title = "Classwise data", x = "Sex (0 = male, 1 = female)", y = "slog(threshold)") + guides(color = "none", fill = "none") + ylim(-3.2, 6.1),
                    pCohenD,      ncol = 1, align = "v", axis = "lr",
                    labels = LETTERS[1:2], rel_heights = c(4, 1)),
          plot_grid(pCorrMat+ labs(title = "Correlation matrix\n"),
                    labels = LETTERS[3]),
          nrow = 1, align = "hv", axis = "lr",rel_widths = c(2, 1)
)

plot_grid(plot_grid(pClasswise + labs(title = "Classwise data", x = "Sex (0 = male, 1 = female)", y = "slog(threshold)") + guides(color = "none", fill = "none") + ylim(-3.2, 6.1),
                    pCorrMat,
                    nrow = 1, align = "hv", axis = "lr",
                    labels = LETTERS[1:2], rel_widths = c(2, 1)
),
plot_grid(pCorrRanked$plotCorrsRanksAbs + labs(title = "Absolute mutual correlations") + coord_flip() +
            theme_light() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = c(.8, .8)) +
            scale_y_discrete(limits = rev),
          labels = LETTERS[3]
),
nrow = 2, align = "hv", axis = "lr"
)

######################### Impute censored data from other data

pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4, ]
CorrKalt <- pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4 & (pCorrRanked$CorrTab$data$col_1 %in% c("Kaelte", "Kaelte_M") | pCorrRanked$CorrTab$data$col_2 %in% c("Kaelte", "Kaelte_M")), ]
CorrWkaelte <- setdiff(unique(c(CorrKalt$col_1, CorrKalt$col_2)), c("Kaelte", "Kaelte_M", "vonFrey", "vonFrey_C",  "CapsHeat",   "CapsvFrey",    "MenthCold"))
CorrFrey <- pCorrRanked$CorrTab$data[abs(pCorrRanked$CorrTab$data$corr) >= 0.4 & (pCorrRanked$CorrTab$data$col_1 %in% c("vonFrey", "vonFrey_C") | pCorrRanked$CorrTab$data$col_2 %in% c("vonFrey", "vonFrey_C")), ]
CorrWFrey <- setdiff(unique(c(CorrFrey$col_1, CorrFrey$col_2)), c("Kaelte", "Kaelte_M", "vonFrey", "vonFrey_C",  "CapsHeat",   "CapsvFrey",    "MenthCold"))

# https://cran.r-project.org/web/packages/lodi/vignettes/lodi.html
PainThresholdsData_transformed_all_toImpute <- PainThresholdsData_transformed_all
sum(PainThresholdsData_transformed_all_toImpute$Kaelte == 0)
sum(PainThresholdsData_transformed_all_toImpute$Kaelte_M == 0)
sum(PainThresholdsData_transformed_all_toImpute$vonFrey >= log(300+1))
sum(PainThresholdsData_transformed_all_toImpute$vonFrey_C >= log(300+1))

PainThresholdsData_transformed_all_toImpute$Kaelte_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$Kaelte == 0, NA, PainThresholdsData_transformed_all_toImpute$Kaelte)
PainThresholdsData_transformed_all_toImpute$Kaelte_M_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$Kaelte_M == 0, NA, PainThresholdsData_transformed_all_toImpute$Kaelte_M)
PainThresholdsData_transformed_all_toImpute$vonFrey_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$vonFrey >=  log(300+1), NA, -PainThresholdsData_transformed_all_toImpute$vonFrey)
PainThresholdsData_transformed_all_toImpute$vonFrey_C_NA <- ifelse(PainThresholdsData_transformed_all_toImpute$vonFrey_C >=  log(300+1), NA, -PainThresholdsData_transformed_all_toImpute$vonFrey_C)
PainThresholdsData_transformed_all_toImpute
PainThresholdsData_transformed_all_toImpute$lod <- 0
PainThresholdsData_transformed_all_toImpute$lodF <- -log(300+1)


clmi.out <- clmi(
  formula = as.formula(paste("Kaelte_NA ~ ", paste(CorrWkaelte, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lod, seed = 42, n.imps = 5
)
Kaelte_imputed <- apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$Kaelte_NA_transform_imputed)), 1, median)
clmi.out <- clmi(
  formula = as.formula(paste("Kaelte_M_NA ~ ", paste(CorrWkaelte, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lod, seed = 42, n.imps = 5
)
Kaelte_M_imputed <- apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$Kaelte_M_NA_transform_imputed)), 1, median)
clmi.out <- clmi(
  formula = as.formula(paste("vonFrey_NA ~ ", paste(CorrWFrey, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lodF, seed = 42, n.imps = 5
)
vonFrey_imputed <- -apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$vonFrey_NA_transform_imputed)), 1, median)
clmi.out <- clmi(
  formula = as.formula(paste("vonFrey_C_NA ~ ", paste(CorrWFrey, collapse = " + "))),
  df = PainThresholdsData_transformed_all_toImpute, lod = lodF, seed = 42, n.imps = 5
)
vonFrey_C_imputed <- -apply(data.frame(lapply(clmi.out$imputed.dfs, function(x) x$vonFrey_C_NA_transform_imputed)), 1, median)

dfImputationEffect <- cbind.data.frame(Cold_all = clmi.out$imputed.dfs[[1]]$Kaelte,
                                       Cold_noncensored = clmi.out$imputed.dfs[[1]]$Kaelte_NA,
                                       Cold_imputed = Kaelte_imputed,
                                       Cold_M_all = clmi.out$imputed.dfs[[1]]$Kaelte_M,
                                       Cold_M_noncensored = clmi.out$imputed.dfs[[1]]$Kaelte_M_NA,
                                       Cold_M_imputed = Kaelte_M_imputed,
                                       vonFrey_all = clmi.out$imputed.dfs[[1]]$vonFrey,
                                       vonFrey_all_noncensored = -clmi.out$imputed.dfs[[1]]$vonFrey_NA,
                                       vonFrey_imputed = vonFrey_imputed,
                                       vonFrey_C_all = clmi.out$imputed.dfs[[1]]$vonFrey_C,
                                       vonFrey_C_all_noncensored = -clmi.out$imputed.dfs[[1]]$vonFrey_C_NA,
                                       vonFrey_C_imputed = vonFrey_C_imputed)
dfImputationEffect_long = reshape2::melt(dfImputationEffect)
dfImputationEffect_long$Stimulus <- rep(c("Cold", "Cold_M", "vonFrey", "vonFrey_C"), each = 3 * nrow(dfImputationEffect))
dfImputationEffect_long$Censored <- rep(c("Original", "Removed", "Imputed"), each = nrow(dfImputationEffect), times = 4)

head(dfImputationEffect_long)

ggplot(data = dfImputationEffect_long) +
  geom_density(aes(x = value, color = Censored, linetype = Censored, fill = Censored), alpha = 0.1) +
  facet_wrap(Stimulus  ~ .) +
  theme_light() +
  theme(
    legend.position = c(.4, .8),
    legend.background = element_rect(colour = "transparent", fill = ggplot2::alpha("white", 0.6)),
    strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black")
  ) +
  labs(x = "sLog threshold", y = "Density [1/sLog trheshold]", color = "Censored values", linetype = "Censored values", title = "Imputation of censored data") +
  scale_color_manual(values = colorblind_pal()(8) ) +
  guides(fill = "none")



par(mfrow = c(4, 3))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$Kaelte)))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$Kaelte_NA)))
plot(density(Kaelte_imputed))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$Kaelte_M)))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$Kaelte_M_NA)))
plot(density(Kaelte_M_imputed))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$vonFrey)))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$vonFrey_NA)))
plot(density(vonFrey_imputed))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$vonFrey_C)))
plot(density(na.omit(clmi.out$imputed.dfs[[1]]$vonFrey_C_NA)))
plot(density(vonFrey_C_imputed))
par(mfrow = c(1, 1))

PainThresholdsData_transformed_imputed <- PainThresholdsData_transformed %>% select(-c(Kaelte, Kaelte_M, MenthCold, vonFrey, vonFrey_C, CapsvFrey))
PainThresholdsData_transformed_imputed$Kaelte_I <- Kaelte_imputed
PainThresholdsData_transformed_imputed$Kaelte_M_I <- Kaelte_M_imputed
PainThresholdsData_transformed_imputed$MenthCold_I <- Kaelte_imputed - Kaelte_M_imputed
PainThresholdsData_transformed_imputed$vonFrey_I <- vonFrey_imputed
PainThresholdsData_transformed_imputed$vonFrey_C_I <- vonFrey_C_imputed
PainThresholdsData_transformed_imputed$CapsvFrey_I <- vonFrey_imputed - vonFrey_C_imputed

#write.csv(PainThresholdsData_transformed_imputed, "PainThresholdsData_transformed_imputed.csv")

############################ Check modal distribution ####################################
GMMtests <- pbmcapply::pbmclapply(names(PainThresholdsData_transformed_imputed), function(column) {
  GMM <- opGMMassessment::opGMMassessment(PainThresholdsData_transformed_imputed[[column]], PlotIt = T, Seed = 42, FitAlg = "normalmixEM", Criterion = "LR")
  StatSig <- ClusterAccuracy(PainThresholdsData_transformed_all$Sex, GMM$Cls)
  if (length(unique(GMM$Cls)) > 1) {
    StatSig1 <- chisq.test(PainThresholdsData_transformed_all$Sex, GMM$Cls)$p.value
  } else
    StatSig1 <- 1
  GMM$Plot <- GMM$Plot + labs(title = paste0(column, ": p-value vs. prior: ", round(StatSig1,5)), y = "Density (unscaled)") + theme(legend.position = c(.05,.8))
  return(list(Accuracy = StatSig, ChisSquare = StatSig1, GMM = GMM))
}, mc.cores = nProc)

names(GMMtests) <- names(PainThresholdsData_transformed_imputed)
ClusterAcc <- lapply(GMMtests, "[[", "Accuracy")
Chis <- lapply(GMMtests, "[[", "ChisSquare")

XX <- lapply(lapply(GMMtests, "[[", "GMM"), function(i) i$Plot)
XX1 <- XX
rep_str = c("vonFrey" = "Punctate pressure (von Frey)",
            "Hitze" = "Heat",
            "Kaelte" = "Cold",
            "Druck" = "Blunt pressure",
            "Strom" = "Electrical",
            "vonFrey_C" = "von Frey + capsaicin",
            "Hitze_C" = "Heat + capsaicin",
            "Kaelte_M" = "Cold + menthol",
            "Cold_M" = "Cold + menthol",
            "CapsHeat" = "Capsaicin effect on heat",
            "CapsvFrey" = "Capsaicin effect on von Frey",
            "MenthCold" = "Menthol effect on cold")
for (i in 1:length(XX)) {
  TitleOld <- XX1[[i]]$labels$title
  TitleNew <- str_replace_all(TitleOld, rep_str)
  XX1[[i]]$labels$title <- TitleNew
}


glist <- lapply(XX1, ggplotGrob)
plot_grid(
  plotlist = glist,
  nrow = 2, align = "hv", axis = "tblr",
  labels = LETTERS[1:14]
)

PainThresholdsData_transformed_imputed_proj <- cbind.data.frame(Target = PainThresholdsData_transformed_all$Sex, PainThresholdsData_transformed_imputed)
PainThresholds <- PainThresholdsData_transformed_imputed_proj
#write.csv(PainThresholds, "PainThresholds.csv")
# PainThresholds <- read.csv(paste0(pfad_o,pfad_u1,"R/", "PainThresholds.csv"), row.names = 1)

pClasswise_imp <- plotClasswiseData(
  X = subset(PainThresholds, select = !names(PainThresholds) %in% "Target"),
  y = PainThresholds$Target, rows = 1
)

PainThresholds_scaled <- PainThresholds
PainThresholds_scaled[,2:ncol(PainThresholds_scaled)] <- scale(PainThresholds_scaled[,2:ncol(PainThresholds_scaled)])
PainThresholdsDown_scaled <- opdisDownsampling::opdisDownsampling(Data = within(PainThresholds_scaled, rm(Target)),
                                                           Cls = PainThresholds_scaled$Target, Size = 0.8*nrow(PainThresholds_scaled), Seed = 42,
                                                           nTrials = 1000000, MaxCores = 30, PCAimportance = F)
write.csv(PainThresholdsDown_scaled$ReducedData, "PainThresholds_scaled_Training.csv")
write.csv(PainThresholdsDown_scaled$RemovedData, "PainThresholds_scaled_Test.csv")





###############################################################################
# Generic Data Analysis Pipeline: Application to pain QST data
#
# This script loads, preprocesses, and analyzes datasets for classification.
# It handles feature selection, correlation analysis, and machine learning
# classification with multiple algorithms and confidence intervals.
###############################################################################

# --- Libraries ----------------------------------------------------------------
library(parallel)
library(randomForest)
library(caret)
library(C50)
library(partykit)
library(pbmcapply)
library(Boruta)
library(reshape2)
library(pROC)
library(dplyr)
library(glmnet)
library(car)
library(tidyr)
library(ggplot2)
library(patchwork)

###############################################################################
# Configuration Parameters
###############################################################################

# External functions
FUNCTIONS_FILE_PATH <- "/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ABCPython/08AnalyseProgramme/R/ABC2way/feature_selection_and_classification_functions.R"

# Dataset name
DATASET_NAME <- "Pheno_125"
EXPERIMENTS_DIR <- "/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ABCPython/08AnalyseProgramme/R/ABC2way/"

# Dataset paths (now only training + validation CSVs)
base_path <- "/home/joern/Dokumente/PainGenesDrugs"
r_path <- "08AnalyseProgramme/R"

train_file <- file.path(base_path, r_path, "PainThresholds_scaled_Training.csv")
val_file <- file.path(base_path, r_path, "PainThresholds_scaled_Test.csv")

# Analysis parameters
SEED <- 42
noise_factor <- 0.2
CORRELATION_METHOD <- "pearson"
CORRELATION_LIMIT <- 0.9
Boruta_tentative_in <- FALSE
use_nyt <- TRUE
tune_RF <- TRUE
tune_KNN <- TRUE
tune_SVM <- TRUE
mtry_12only <- FALSE

use_curated <- FALSE
use_roc_auc <- FALSE
training_and_validation_subsplits <- TRUE
TRAINING_PARTITION_SIZE <- 0.8
VALIDATION_PARTITION_SIZE <- 0.8

max_iterations <- 5
RUN_ONE_ADDITIONAL_ITERATION <- FALSE

# Dataset-specific configuration
DATASET_COLUMN_NAMES <- c(
  "Heat", "Pressure", "Current", "Heat_Capsaicin",
  "Capsaicin_Effect_Heat", "Cold", "Cold_Menthol", "Menthol_Effect_Cold",
  "vonFrey", "vonFrey_Capsaicin", "Capsaicin_Effect_vonFrey"
)

CURATED_COLUMN_NAMES <- c(
  "Heat", "Pressure", "Current", "Heat_Capsaicin",
  "Cold", "Cold_Menthol",
  "vonFrey", "vonFrey_Capsaicin"
)

COLUMNS_COLINEAR <- NULL # Will be determined later

# Noise addition parameters
ADD_NOISE_COLUMN <- TRUE
NOISE_COLUMN_SOURCE <- "Pressure"
NOISE_COLUMN_NAME <- "Pressure2"


###############################################################################
# Load External Functions
###############################################################################
rename_pain_data_columns <- function(data, new_names) {
  # Assign new column names for consistency.
  names(data) <- new_names
  data
}

###############################################################################
# Check data files existence
###############################################################################
if (!file.exists(train_file)) {
  stop(paste("Training data file not found:", train_file))
}
if (!file.exists(val_file)) {
  stop(paste("Validation data file not found:", val_file))
}


###############################################################################
# Load data and modify files when needed
###############################################################################

# --- Define paths ---
base_path <- "/home/joern/Dokumente/PainGenesDrugs"
r_path <- "08AnalyseProgramme/R"

train_file <- file.path(base_path, r_path, "PainThresholds_scaled_Training.csv")
val_file <- file.path(base_path, r_path, "PainThresholds_scaled_Test.csv")

# --- Load datasets ---
train_df <- read.csv(train_file, row.names = 1)
val_df <- read.csv(val_file, row.names = 1)

# --- Sizes to later reseparate ---
n_train <- nrow(train_df)
n_val <- nrow(val_df)

# --- Extract features + targets ---
train_features <- train_df[, - ncol(train_df)]
train_target <- train_df[, ncol(train_df)]

val_features <- val_df[, - ncol(val_df)]
val_target <- val_df[, ncol(val_df)]

# --- Combine into pain_data + target_data ---
pain_data <- rbind(train_features, val_features)
target_data <- c(train_target, val_target)

# --- Rename columns once, for the combined dataset ---
pain_data <- rename_pain_data_columns(pain_data, DATASET_COLUMN_NAMES)

# --- Add noise ONCE (so all subsets share the same noise values) ---
set.seed(SEED)
pain_data[[NOISE_COLUMN_NAME]] <- pain_data[[NOISE_COLUMN_SOURCE]] +
  rnorm(
    n = nrow(pain_data),
    mean = 0,
    sd = abs(pain_data[[NOISE_COLUMN_SOURCE]]) * noise_factor
  )

# --- Reseparate into train/validation, consistent with original splits ---
training_data_original <- pain_data[1:n_train,]
validation_data_original <- pain_data[(n_train + 1):(n_train + n_val),]

training_target <- target_data[1:n_train]
validation_target <- target_data[(n_train + 1):(n_train + n_val)]

# --- Optional curation after renaming + noise ---
if (use_curated) {
  training_data_original <- training_data_original[, CURATED_COLUMN_NAMES]
  validation_data_original <- validation_data_original[, CURATED_COLUMN_NAMES]
  pain_data <- pain_data[, CURATED_COLUMN_NAMES]
}

# Save data file for use in other analyses
Pheno_125_prepared_data <- cbind.data.frame(Target = target_data, pain_data)
