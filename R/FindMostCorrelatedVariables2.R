################################################################################
# Find Most Correlated Variables and Produce Visualizations
################################################################################
#
# This function computes pairwise correlations between all variables in a dataset
# and creates multiple visualizations to explore correlation patterns.
#
# Main use case: Identify which variables are highly correlated with each other,
# especially useful for selecting predictor variables for imputation or feature
# selection tasks.
#
# Parameters:
#   Data          - Data frame with numeric variables
#   Corrmethod    - Correlation method: "spearman" (default), "pearson", or "kendall"
#   TargetVariable- Optional: Name of variable to highlight in plots
#   Corrlimit     - Minimum absolute correlation to include in plots (default: 0)
#   Flipped       - If TRUE, flip coordinates of correlation plot (default: TRUE)
#
# Returns:
#   List containing:
#   - CorrTab: Correlation table with confidence intervals
#   - plotCoors: Basic correlation plot
#   - plotCorrsRanks: Ranked correlations plot (signed)
#   - plotCorrsRanksAbs: Ranked absolute correlations plot
#   - CorrelatedVars: Per-variable correlation summary plot
#
################################################################################

FindMostCorrelatedVariables2 <- function(Data, Corrmethod = "spearman", TargetVariable = NA, Corrlimit = 0, Flipped = TRUE) {
  # Reference: https://alastairrushworth.github.io/Calculating-and-visualising-correlation-coefficients-with-inspectdf/
  library(inspectdf)
  library(dplyr)

  # Helper function: Calculate 95% confidence intervals as quantiles
  # Used for visualizing uncertainty in correlation estimates
  quantiles_95 <- function(x) {
    r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }

  # Compute all pairwise correlations using inspectdf package
  # inspect_cor calculates correlation matrix with confidence intervals
  Corrs <- Data %>%
    inspect_cor(method = Corrmethod) %>%
    show_plot(col_palette = 4)

  # Create basic correlation plot with light theme
  pCorrs <-
    Corrs +
    theme_light() +
    theme(axis.text.y = element_text(size = 5)) +
    ggtitle("Correlations")

  # Reshape correlation data from pairwise format to long format
  # This allows aggregation of correlations by variable
  xdf1  <- subset(Corrs$data, select = c("col_1", "corr"))
  xdf2  <- subset(Corrs$data, select = c("col_2", "corr"))
  names(xdf1) <- c("variable", "corr")
  names(xdf2) <- c("variable", "corr")
  CorrsAll <- rbind.data.frame(xdf1, xdf2)

  # Calculate median and 95% CI of absolute correlations for each variable
  # This summarizes how strongly each variable correlates with all others
  m1 <- data.frame(aggregate(abs(CorrsAll$corr), list(CorrsAll$variable), median))
  m1 <- cbind.data.frame(m1, aggregate(CorrsAll$corr, list(CorrsAll$variable), function(x) quantile(abs(x), probs = c(0.025)))[2])
  m1 <- cbind.data.frame(m1, aggregate(CorrsAll$corr, list(CorrsAll$variable), function(x) quantile(abs(x), probs = c(0.975)))[2])
  names(m1) <- c("variable", "corr", "lower", "upper")
  m1$Col <- .1

  # Create plot showing per-variable correlation strength
  # Variables ordered by median absolute correlation
  pCorrelatedVars <-
    ggplot(data = m1) +
    geom_crossbar(aes(x = corr, y = reorder(variable,-abs(corr)), xmin = lower, xmax = upper, color = factor(Col), fill = factor(Col)),
                  alpha =.3, width = 0.8, size = 0.3) +
    scale_fill_manual(values = "dodgerblue") +
    scale_color_manual(values = "blue") +
    scale_alpha_manual( values  = 0.1 ) +
    #stat_summary(fun.data = quantiles_95, geom = "boxplot", fill = "dodgerblue", alpha = 0.3, width = 0.5, position = "dodge") +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1 )) +
    labs(title = "Correlations and 95% CI") +
    guides(color = "none", fill = "none")

  # Optional: Flip coordinates for better readability
  if (Flipped == TRUE)
    pCorrelatedVars <- pCorrelatedVars + coord_flip()


  # Prepare data for pairwise correlation plots
  # Add indicators for whether target variable is involved in each pair
  CorResData <- data.frame(Corrs$data)
  CorResData$TargetPresent <- ifelse(CorResData$col_1 %in% TargetVariable | CorResData$col_2 %in% TargetVariable, "Yes", "No")
  CorResData$TargetPresentAlpha <- ifelse(CorResData$col_1 %in% TargetVariable | CorResData$col_2 %in% TargetVariable, 0.4, 0.2)

  # Create absolute correlation version
  # For negative correlations, adjust CI bounds to reflect absolute values
  CorResDataAbs <- CorResData
  CorResDataAbs$diff_to_lower <- abs(CorResDataAbs$corr - CorResDataAbs$lower)
  CorResDataAbs$diff_to_upper <- abs(CorResDataAbs$corr - CorResDataAbs$upper)
  CorResDataAbs$corr <- abs(CorResDataAbs$corr)
  CorResDataAbs$lower[CorResDataAbs$sign == "Negative"] <- CorResDataAbs$corr[CorResDataAbs$sign == "Negative"] - CorResDataAbs$diff_to_upper[CorResDataAbs$sign == "Negative"]
  CorResDataAbs$upper[CorResDataAbs$sign == "Negative"] <- CorResDataAbs$corr[CorResDataAbs$sign == "Negative"] + CorResDataAbs$diff_to_lower[CorResDataAbs$sign == "Negative"]
  # Filter by minimum correlation threshold
  CorResDataAbs <- CorResDataAbs[CorResDataAbs$corr >= Corrlimit, ]

  # Create plot showing all pairwise correlations (signed)
  # Highlights pairs involving the target variable if specified
  pCorrsRanks <- ggplot(data = CorResData) +
    geom_crossbar(aes(x = corr, y = pair, xmin = lower, xmax = upper, color = factor(TargetPresent), fill = factor(TargetPresent)),
                  alpha = CorResData$TargetPresentAlpha, width = 0.8, size = 0.3) +
    theme_light() +
    theme(legend.position = c(.2, .8), axis.text.y = element_text(size = 10)) +
    scale_color_manual(values = c("blue", "darkblue")) +
    scale_fill_manual(values = c("dodgerblue", "darkblue"))

  # Add appropriate labels depending on whether target variable is specified
  if (!missing(TargetVariable)) {
    pCorrsRanks <- pCorrsRanks + labs(title = "Sorted mutual correlations", x = "Correlation coefficient (Spearman's rho)", y = "Variables pairs", fill = paste0(TargetVariable, " in pair")) +
      guides(color = "none")
  } else {
    pCorrsRanks <- pCorrsRanks + labs(title = "Sorted mutual correlations", x = "Correlation coefficient (Spearman's rho)", y = "Variables pairs") +
      guides(color = "none", fill = "none")
  }

  # Create plot showing absolute correlations
  # Color-codes by correlation sign (positive/negative)
  pCorrsRanksAbs <- ggplot(data = CorResDataAbs) +
    geom_crossbar(aes(x = corr, y = pair, xmin = lower, xmax = upper,
                      color = factor(sign), fill = factor(TargetPresent)), alpha = CorResDataAbs$TargetPresentAlpha, width = 0.8, size = 0.3) +
    theme_light() +
    theme(legend.position = c(.2, .8), axis.text.y = element_text(size = 10)) +
    scale_fill_manual(values = c("dodgerblue", "darkblue")) +
    labs(color = "Sign")

  # Handle color coding based on presence of both positive and negative correlations
  if(length(unique(CorResData$sign)) > 1) {
    pCorrsRanksAbs <- pCorrsRanksAbs +
      scale_color_manual(values = c("firebrick", "forestgreen"), labels = c("negative", "positive"))
  }  else {
    pCorrsRanksAbs <- pCorrsRanksAbs +
      scale_color_manual(values = "blue") +
      guides(color = "none")
  }

  # Add appropriate labels for absolute correlation plot
  if (!missing(TargetVariable)) {
    pCorrsRanksAbs <- pCorrsRanksAbs + labs(title = "Sorted mutual correlations", x = "Correlation coefficient (Spearman's rho)", y = "Variables pairs", fill = paste0(TargetVariable, " in pair"))
  } else {
    pCorrsRanksAbs <- pCorrsRanksAbs + labs(title = "Sorted mutual correlations", x = "Correlation coefficient (Spearman's rho)", y = "Variables pairs") +
      guides(fill = "none")
  }

  # Return all correlation results and plots as a list
  return(list(CorrTab = Corrs, plotCoors = pCorrs, plotCorrsRanks = pCorrsRanks, plotCorrsRanksAbs = pCorrsRanksAbs, CorrelatedVars = pCorrelatedVars))
}
