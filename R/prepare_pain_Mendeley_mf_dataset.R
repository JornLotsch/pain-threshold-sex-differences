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

facet_labs <- c(
  "vonFrey"           = "g",
  "Heat"              = "°C",
  "Cold"              = "°C",
  "Pressure"          = "N/cm2",
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
  geom_point(aes(color = Sex, shape = TrainingValidation),
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
    legend.position  = "bottom",
    legend.direction = "horizontal"
  ) +
  ggthemes::scale_color_colorblind()
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

plot_metadata <- ggplot() +
  geom_tile(
    data = age_long,
    aes(x = variable, y = ID, fill = value),
    alpha = .6
  ) +
  scale_fill_gradient(
    name = "Age",
    # low = "cornsilk",
    # high = "cornsilk4"
    low = ggthemes::colorblind_pal()(8)[2],
    high = ggthemes::colorblind_pal()(8)[1]
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
    values = c(
      "Nonsmoker" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[2], .6),
      "Smoker" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[3], .6),
      "Male" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[4], .6),
      "Female" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[7], .6),
      "M" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[5], .6),
      "F" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[8], .6),
      "Training" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[1], .5),
      "Validation" = ggplot2::alpha(ggthemes::colorblind_pal()(8)[6], .6)
    )
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 3),
    panel.grid = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal", legend.byrow = TRUE
  )

print(plot_metadata)

ggsave(plot = plot_metadata, filename = "plot_metadata.svg", width = 9, height = 9)
ggsave(plot = plot_metadata, filename = "plot_metadata.png", width = 9, height = 9, dpi = 300)
