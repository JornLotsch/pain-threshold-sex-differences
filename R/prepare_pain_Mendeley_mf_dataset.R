################################################################################
# Assemble the files deposited at Mendeley Data (doi: 10.17632/y6v3rtpgps.1)
################################################################################
#
# Final step of the pipeline (after Python/Read_explore_PainThresholds.py and
# R/preprocessing_pain_Mendeley_mf_dataset.R). Writes:
#   exp_pain_data_orig.csv        - participant-level threshold estimates
#                                   (derived from repeated measurements or
#                                   logistic regression; not trial-level data),
#                                   untransformed, in original units
#   exp_pain_data_transformed.csv - ln(x + 1) transformed, imputed, z-scaled
#                                   data from PainThresholds_scaled_*.csv
#   exp_pain_metadata.csv         - metadata incl. training/validation assignment
# Also produces Figure 1 (plot_orig_data), the metadata plot, and quality-control
# output (ID consistency, observed ranges, censored values).
#
# Units: von Frey in g, heat and cold in °C, pressure in N/cm², current in mA.
#
################################################################################

############### Libraries ##############################

library(ggplot2)
library(reshape2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(ggnewscale)


############### Constants ##############################

# File paths
INPUT_FILE_PATH <- "/home/joern/Dokumente/PainGenesDrugs/09Originale/Phaeno_125.csv"
OUTPUT_DIR <- "/home/joern/Aktuell/QST_Mendeley_MF/DataSetPublished"

# Variable definitions
PAIN_TESTS_NAMES <- c(
  "vonFrey", "Hitze", "Kaelte", "Druck", "Strom",
  "vonFrey_C", "Hitze_C", "Kaelte_M"
)

PAIN_TESTS_NAMES_ENGLISH <- c(
  "vonFrey", "Heat", "Cold", "Pressure", "Current",
  "vonFrey_Capsaicin", "Heat_Capsaicin", "Cold_Menthol"
)

# Metadata variables
METADATA <- c("Probandennummer", "Alter", "Raucher", "Geschlecht", "Pruefer")
METADATA_ENGLISH <- c("ID", "Age", "Smoker", "Sex", "Tester")

# Dataset paths
BASE_PATH <- "/home/joern/Dokumente/PainGenesDrugs"
R_PATH <- "08AnalyseProgramme/R"

TRAIN_FILE <- file.path(BASE_PATH, R_PATH, "PainThresholds_scaled_Training.csv")
VAL_FILE <- file.path(BASE_PATH, R_PATH, "PainThresholds_scaled_Test.csv")

# Dataset-specific configuration
DATASET_COLUMN_NAMES <- c(
  "Heat", "Pressure", "Current", "Heat_Capsaicin",
  "Capsaicin_Effect_Heat", "Cold", "Cold_Menthol", "Menthol_Effect_Cold",
  "vonFrey", "vonFrey_Capsaicin", "Capsaicin_Effect_vonFrey"
)


############### Utility Functions ##############################

#' Rename columns in a dataset
#' @param data Data frame to rename
#' @param new_names Character vector with new column names
#' @return Data frame with renamed columns
rename_pain_data_columns <- function(data, new_names) {
  names(data) <- new_names
  return(data)
}


############### Data Validation ##############################

# Check presence of necessary files
if (!file.exists(TRAIN_FILE)) {
  stop(paste("Training data file not found:", TRAIN_FILE))
}

if (!file.exists(VAL_FILE)) {
  stop(paste("Validation data file not found:", VAL_FILE))
}


############### Load and Prepare Data ##############################

# Load primary dataset
cat("Loading raw phenotype data...\n")
Phaeno_125 <- read.csv(INPUT_FILE_PATH)

# Extract and rename original pain data
Phaeno_125_orig_pain_data <- Phaeno_125[, c("Probandennummer", PAIN_TESTS_NAMES)]
Phaeno_125_orig_pain_data_renamed <- rename_pain_data_columns(
  Phaeno_125_orig_pain_data, c("ID", PAIN_TESTS_NAMES_ENGLISH)
)

# Extract and rename metadata
Phaeno_125_metadata <- Phaeno_125[, METADATA]
Phaeno_125_metadata_renamed <- rename_pain_data_columns(Phaeno_125_metadata, METADATA_ENGLISH)

# Standardize metadata values
Phaeno_125_metadata_renamed$Tester <- ifelse(Phaeno_125_metadata_renamed$Tester == "Karin", "F", "M") #as.integer(as.factor(Phaeno_125_metadata_renamed$Tester))
Phaeno_125_metadata_renamed$Smoker <- ifelse(Phaeno_125_metadata_renamed$Smoker == "J", "Smoker", "Nonsmoker")
Phaeno_125_metadata_renamed$Sex <- ifelse(Phaeno_125_metadata_renamed$Sex == "W","Female", "Male")


############### Load Training and Validation Data ##############################

cat("Loading training and validation datasets...\n")
train_df <- read.csv(TRAIN_FILE, row.names = 1)
val_df <- read.csv(VAL_FILE, row.names = 1)

n_train <- nrow(train_df)
n_val <- nrow(val_df)

# Separate features and targets
train_features <- train_df[, -ncol(train_df)]
train_target <- train_df[, ncol(train_df)]

val_features <- val_df[, -ncol(val_df)]
val_target <- val_df[, ncol(val_df)]

# Combine datasets
pain_data <- rbind(train_features, val_features)
target_data <- c(train_target, val_target)

# Label participants as Training or Validation
Phaeno_125_metadata_renamed$TrainingValidation <- ifelse(
  Phaeno_125_metadata_renamed$ID %in% as.integer(rownames(train_features)),
  "Training",
  "Validation"
)


############### Data Postprocessing ##############################

# Rename columns of combined dataset
pain_data <- rename_pain_data_columns(pain_data, DATASET_COLUMN_NAMES)

# Split combined dataset back into original subsets
training_data_original <- pain_data[1:n_train, ]
validation_data_original <- pain_data[(n_train + 1):(n_train + n_val), ]

training_target <- target_data[1:n_train]
validation_target <- target_data[(n_train + 1):(n_train + n_val)]


############### Data Alignment and Export ##############################

# Align and sort data according to metadata
Pheno_125_prepared_data <- cbind.data.frame(Target = target_data, pain_data)
Pheno_125_prepared_data <- cbind.data.frame(ID = rownames(Pheno_125_prepared_data), Pheno_125_prepared_data)

Pheno_125_prepared_data_sorted <- Pheno_125_prepared_data[
  match(Phaeno_125_metadata_renamed$ID, Pheno_125_prepared_data$ID),
]

# Basic verification plot
plot(
  Pheno_125_prepared_data_sorted$Target ~
    (as.integer(as.factor(Phaeno_125_metadata_renamed$Sex)) - 1),
  main = "Target vs Sex",
  xlab = "Sex (0 = F, 1 = M)",
  ylab = "Target"
)

# Remove target before final export
Pheno_125_prepared_data_final <- Pheno_125_prepared_data_sorted[
  , !names(Pheno_125_prepared_data_sorted) %in% c("Target")
]


############### Write Output Files ##############################

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Write CSV files
write.csv(
  Phaeno_125_orig_pain_data_renamed,
  file = file.path(OUTPUT_DIR, "exp_pain_data_orig.csv"),
  row.names = FALSE
)

write.csv(
  Pheno_125_prepared_data_final,
  file = file.path(OUTPUT_DIR, "exp_pain_data_transformed.csv"),
  row.names = FALSE
)

write.csv(
  Phaeno_125_metadata_renamed,
  file = file.path(OUTPUT_DIR, "exp_pain_metadata.csv"),
  row.names = FALSE
)

cat("Data files successfully written to:", OUTPUT_DIR, "\n")
cat("Created files:\n")
cat("- exp_pain_metadata.csv\n")
print(names(Phaeno_125_metadata_renamed))
cat("- exp_pain_data_orig.csv\n")
print(names(Phaeno_125_orig_pain_data_renamed))
cat("- exp_pain_data_transformed.csv\n")
print(names(Pheno_125_prepared_data_final))


# ================= Plot data ==================

# Original data

plot_df_orig_data <- cbind.data.frame(Sex = Phaeno_125_metadata_renamed$Sex, TrainingValidation = Phaeno_125_metadata_renamed$TrainingValidation,
                                      Phaeno_125_orig_pain_data_renamed[,-1])

plot_df_orig_data_long <- reshape2::melt(plot_df_orig_data, id.vars = c("Sex", "TrainingValidation"))

# Combined Sex x TrainingValidation group, so that the legend shows all four
# combinations unambiguously (revision: previously separate color and shape legends).
# Four colors: hue codes sex, darker (training) / lighter (validation) shade codes
# the dataset; checked for deutan, protan and tritan color vision deficiency
# (colorspace::deutan etc.). Shapes additionally code the dataset.
plot_df_orig_data_long$Group <- factor(
  paste(plot_df_orig_data_long$Sex, tolower(plot_df_orig_data_long$TrainingValidation), sep = ", "),
  levels = c("Female, training", "Female, validation", "Male, training", "Male, validation")
)
group_colors <- c(
  "Female, training"   = "#C24100", # dark vermilion
  "Female, validation" = "#F39B3C", # light orange
  "Male, training"     = "#08519C", # dark blue
  "Male, validation"   = "#5AAEE6"  # light blue
)
group_shapes <- c(
  "Female, training"   = 16,
  "Female, validation" = 17,
  "Male, training"     = 16,
  "Male, validation"   = 17
)

facet_labs <- c(
  "vonFrey"           = "g",
  "Heat"              = "°C",
  "Cold"              = "°C",
  "Pressure"          = "N/cm²",
  "Current"           = "mA",
  "vonFrey_Capsaicin" = "g",
  "Heat_Capsaicin"    = "°C",
  "Cold_Menthol"      = "°C"
)

set.seed (42)
plot_orig_data <-
  ggplot(plot_df_orig_data_long,
         aes(x = variable, y = value)) +
  geom_violin(alpha = 0.2, fill = "cornsilk", color = "cornsilk4") +
  geom_boxplot(alpha = 0.5, width = .2, color = "cornsilk4") +
  geom_point(aes(color = Group, shape = Group),
             size = 2,
             position = position_jitter(width = 0.1),
             show.legend = TRUE) +
  facet_wrap(
    variable ~ .,
    scales = "free",
    nrow = 2,
    strip.position = "left",
    labeller = as_labeller(facet_labs)
  ) +
  theme_light() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(colour = "black"),
    strip.placement  = "outside",
    axis.title.y     = element_blank(),
    axis.title.x     = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal"
  ) +
  scale_color_manual(name = "Sex, dataset", values = group_colors) +
  scale_shape_manual(name = "Sex, dataset", values = group_shapes)
  # ggthemes::scale_color_colorblind()
  # scale_color_manual(values = c("dodgerblue", "chartreuse3"))

print(plot_orig_data)

ggsave(plot = plot_orig_data, filename = "plot_orig_data.svg", width = 9, height = 9)
ggsave(plot = plot_orig_data, filename = "plot_orig_data.png", width = 9, height = 9, dpi = 300)


# Metadata

plot_df_metadata <- Phaeno_125_metadata_renamed

age_long <- plot_df_metadata %>%
  mutate(ID = factor(ID)) %>%
  select(ID, Age) %>%
  mutate(variable = "Age") %>%
  rename(value = Age)

fac_long <- plot_df_metadata %>%
  mutate(ID = factor(ID)) %>%
  select(ID, Smoker, Sex, Tester, TrainingValidation) %>%
  pivot_longer(
    cols = -ID,
    names_to = "variable",
    values_to = "value"
  )

# Revision: colors consistent with Figure 1 (sex: orange/blue; training/validation:
# dark/light), tester explicitly labeled, legend ordered as the columns
plot_metadata <- ggplot() +
  geom_tile(
    data = age_long,
    aes(x = variable, y = ID, fill = value)
  ) +
  scale_fill_gradient(
    name = "Age",
    # low = ggthemes::colorblind_pal()(8)[2],
    # high = ggthemes::colorblind_pal()(8)[1]
    low = "#EFEDF5",
    high = "#54278F"
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = fac_long,
    aes(x = variable, y = ID, fill = value)
  ) +
  scale_fill_manual(
    name = NULL,
    # values = c(
    #   "Nonsmoker" = "cornsilk2",
    #   "Smoker" = "cornsilk4",
    #   "Male" = "cornsilk2",
    #   "Female" = "cornsilk4",
    #   "M" = "cornsilk2",
    #   "F" = "cornsilk4",
    #   "Training" = "cornsilk2",
    #   "Validation" = "cornsilk4"
    # )
    # values = c(
    #   "Nonsmoker" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[2], .6),
    #   "Smoker" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[3], .6),
    #   "Male" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[4], .6),
    #   "Female" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[7], .6),
    #   "M" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[5], .6),
    #   "F" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[8], .6),
    #   "Training" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[1], .5),
    #   "Validation" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[6], .6)
    # )
    values = c(
      "Female"     = ggplot2::alpha("#C24100", .75),
      "Male"       = ggplot2::alpha("#08519C", .75),
      "Nonsmoker"  = ggplot2::alpha("#009E73", .3),
      "Smoker"     = ggplot2::alpha("#009E73", .85),
      "F"          = ggplot2::alpha("#CC79A7", .7),
      "M"          = ggplot2::alpha("#F0E442", .8),
      "Training"   = "grey35",
      "Validation" = "grey80"
    ),
    breaks = c("Female", "Male", "Nonsmoker", "Smoker", "F", "M", "Training", "Validation"),
    labels = c("Female", "Male", "Nonsmoker", "Smoker", "Tester F", "Tester M", "Training", "Validation"),
    guide = guide_legend(nrow = 2, byrow = FALSE)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 3),
    panel.grid = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal"
  )

print(plot_metadata)

ggsave(plot = plot_metadata, filename = "plot_metadata.svg", width = 9, height = 9)
ggsave(plot = plot_metadata, filename = "plot_metadata.png", width = 9, height = 9, dpi = 300)


# ================= Quality control: ID consistency ==================

# All three deposited files must contain the same unique IDs in the same order
ids_orig        <- as.character(Phaeno_125_orig_pain_data_renamed$ID)
ids_transformed <- as.character(Pheno_125_prepared_data_final$ID)
ids_metadata    <- as.character(Phaeno_125_metadata_renamed$ID)

id_checks <- c(
  n_subjects                     = length(ids_orig),
  ids_unique                     = !anyDuplicated(ids_orig),
  ids_no_missing                 = !anyNA(ids_orig),
  orig_vs_metadata_identical     = identical(ids_orig, ids_metadata),
  orig_vs_transformed_identical  = identical(ids_orig, ids_transformed)
)
cat("ID consistency checks:\n")
print(id_checks)
print(table(Phaeno_125_metadata_renamed$TrainingValidation, Phaeno_125_metadata_renamed$Sex))


# ================= Count censored data ==================

# Censored values imputed in preprocessing_pain_Mendeley_mf_dataset.R:
# Cold and Cold_Menthol at 0 °C (left-censored), vonFrey and vonFrey_Capsaicin
# at 300 g (right-censored). Other values at a theoretical limit are counted
# here but were not imputed.

# Theoretical (device / protocol) limits of each pain test
PAIN_TESTS_LIMITS <- list(
  "vonFrey"           = c(min = 0.008, max = 300), # g
  "Heat"              = c(min = 32, max = 50), # °C
  "Cold"              = c(min = 0, max = 32), # °C
  "Pressure"          = c(min = 0, max = Inf), # N/cm2
  "Current"           = c(min = 0, max = 20), # mA
  "vonFrey_Capsaicin" = c(min = 0.008, max = 300), # g
  "Heat_Capsaicin"    = c(min = 32, max = 50), # °C
  "Cold_Menthol"      = c(min = 0, max = 32)  # °C
)

#' Count values at the theoretical limits of a variable
#' @param x Numeric vector of measured values
#' @param limits Named numeric vector with elements "min" and "max"
#' @param tol Tolerance for the floating point comparison
#' @return Named integer vector with counts at min and max
count_censored <- function(x, limits, tol = 1e-8) {
  x <- x[!is.na(x)]
  c(
    n_valid  = length(x),
    n_at_min = sum(abs(x - limits[["min"]]) < tol),
    n_at_max = sum(abs(x - limits[["max"]]) < tol)
  )
}

censored_counts <- mapply(
  FUN = count_censored,
  x = Phaeno_125_orig_pain_data_renamed[, PAIN_TESTS_NAMES_ENGLISH],
  limits = PAIN_TESTS_LIMITS[PAIN_TESTS_NAMES_ENGLISH]
)

censored_counts_df <- data.frame(
  Variable        = PAIN_TESTS_NAMES_ENGLISH,
  Theoretical_Min = sapply(PAIN_TESTS_LIMITS[PAIN_TESTS_NAMES_ENGLISH], `[[`, "min"),
  Theoretical_Max = sapply(PAIN_TESTS_LIMITS[PAIN_TESTS_NAMES_ENGLISH], `[[`, "max"),
  Observed_Min    = sapply(Phaeno_125_orig_pain_data_renamed[, PAIN_TESTS_NAMES_ENGLISH], min, na.rm = TRUE),
  Observed_Max    = sapply(Phaeno_125_orig_pain_data_renamed[, PAIN_TESTS_NAMES_ENGLISH], max, na.rm = TRUE),
  t(censored_counts),
  row.names = NULL
)
censored_counts_df$Percent_at_min <- round(100 * censored_counts_df$n_at_min / censored_counts_df$n_valid, 1)
censored_counts_df$Percent_at_max <- round(100 * censored_counts_df$n_at_max / censored_counts_df$n_valid, 1)

cat("Observed ranges and censored data (values at theoretical min / max):\n")
print(censored_counts_df)

write.csv(censored_counts_df, "qc_ranges_censored_counts.csv", row.names = FALSE)


# ================= Quality control: training vs. validation ==================

# Imputation and z-standardization were performed on all 125 subjects before the
# training/validation split (not considered in the original analyses). Neither step
# used the class variable (sex). Here, the size of the resulting dependence is
# quantified: (i) means and SDs of the z-scaled variables in both subsets, and
# (ii) the change of the validation values if standardization had used the
# training subset only. As the data are z-scaled, differences are in SD units,
# i.e., 0.01 corresponds to 1% of the standard deviation.
transformed_vars <- setdiff(names(Pheno_125_prepared_data_final), "ID")
is_training <- Phaeno_125_metadata_renamed$TrainingValidation == "Training"
X_training   <- Pheno_125_prepared_data_final[is_training, transformed_vars]
X_validation <- Pheno_125_prepared_data_final[!is_training, transformed_vars]

split_stats_df <- data.frame(
  Variable        = transformed_vars,
  Mean_training   = colMeans(X_training),
  SD_training     = apply(X_training, 2, sd),
  Mean_validation = colMeans(X_validation),
  SD_validation   = apply(X_validation, 2, sd),
  p_Wilcoxon      = sapply(transformed_vars, function(v) wilcox.test(X_training[[v]], X_validation[[v]])$p.value),
  p_KS            = sapply(transformed_vars, function(v) suppressWarnings(ks.test(X_training[[v]], X_validation[[v]])$p.value)),
  row.names = NULL
)

# Validation values standardized with training-only parameters
X_validation_rescaled <- sweep(sweep(X_validation, 2, split_stats_df$Mean_training), 2, split_stats_df$SD_training, "/")
rescaling_change <- abs(as.matrix(X_validation_rescaled) - as.matrix(X_validation))
split_stats_df$Max_change_validation    <- apply(rescaling_change, 2, max)
split_stats_df$Median_change_validation <- apply(rescaling_change, 2, median)

cat("Training vs. validation (z-scaled data):\n")
print(split_stats_df, digits = 3)
cat("Change of validation z-values with training-only standardization: median",
    round(100 * median(rescaling_change), 1), "%, maximum", round(100 * max(rescaling_change), 1), "% of one SD\n")
cat("Rank order of validation values unchanged:",
    all(sapply(transformed_vars, function(v) identical(rank(X_validation_rescaled[[v]]), rank(X_validation[[v]])))), "\n")

write.csv(split_stats_df, "qc_training_validation_comparison.csv", row.names = FALSE)

# Verification that the deposited transformed data and training/validation
# assignment can be regenerated from the deposited files:
# see verify_reproduction_from_deposited_data.R
