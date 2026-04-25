# Figure S2: Correlation of circulating metabolite levels with quantitative IFNγ in follow-up IGRA (week 14)
#
# Linear regression was used to identify baseline metabolites whose levels
# predicted log2 fold-change of TBAg IFNγ from baseline to week 14 follow-up,
# adjusting for age, sex, BMI, presence of BCG scar, and exposure risk score.
# This analysis uses ALL subjects (not just the strict subset).
#
# Output: Volcano plot with regression coefficient (beta) on x-axis and
#         -log10(p-value) on y-axis.

# load packages -----------------------------------------------------------
library(pacman)
pacman::p_load(dplyr, tidyverse, ggplot2, EnhancedVolcano)

# Load data ---------------------------------------------------------------
source("../00_load_data.R")

# Set output directory
OUTPUT_DIR <- file.path(REPO_DIR, "output/")
FIGURE_DIR <- file.path(REPO_DIR, "figures/")

create_output_dirs()

# Load metabolite name mapping
name <- read.csv(METABOLITE_NAME_MAP)
colnames(name) <- c("metabolites", "name")

# Prepare data ------------------------------------------------------------
# This analysis uses ALL subjects (not just strict subset) with valid Log2FC_TBonly
# Log2FC_TBonly = log2(TBAg at week 14) - log2(TBAg at baseline)
# This matches the exact analysis in analysisHMDB_fix.withBCGscar.v2.R

if (!("Log2FC_TBonly" %in% colnames(klinis))) {
  stop("Log2FC_TBonly not found in clinical data. Ensure INFECT_log2fc_tbag.csv is in data/ and 00_load_data.R has been sourced.")
}

# Use all subjects with valid Log2FC_TBonly
klinis_hc <- klinis %>%
  dplyr::filter(!is.na(Log2FC_TBonly))

# Subset metabolites to subjects with valid Log2FC
met_hc <- metabolites[rownames(metabolites) %in% klinis_hc$IdInfect, , drop = FALSE]

# Align order by IdInfect
klinis_hc <- klinis_hc %>% dplyr::arrange(IdInfect)
met_hc <- met_hc[order(rownames(met_hc)), , drop = FALSE]

# Verify alignment
stopifnot(all(rownames(met_hc) == klinis_hc$IdInfect))

# Get metabolite columns
metabolite_cols <- colnames(met_hc)

# Ensure covariates are properly typed
klinis_hc$agecontact <- as.numeric(klinis_hc$agecontact)
klinis_hc$BMI_contact <- as.numeric(klinis_hc$BMI_contact)
klinis_hc$sex_contact <- as.factor(klinis_hc$sex_contact)
klinis_hc$BCG_scar <- as.factor(klinis_hc$BCG_scar)
klinis_hc$exp_score <- as.numeric(klinis_hc$exp_score)

# Linear regression: Log2FC_TBonly ~ metabolite + covariates ---------------
# Outcome: log2 fold-change of TBAg IFNγ from baseline to week 14
# Predictor: each baseline metabolite
# Covariates: age, sex, BMI, BCG scar, exposure risk score

result_lm <- data.frame(
  metabolites = character(),
  estimate = numeric(),
  std_error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (m in metabolite_cols) {
  # Build data frame for this metabolite
  df <- data.frame(
    outcome = as.numeric(klinis_hc$Log2FC_TBonly),
    metabolite = as.numeric(met_hc[[m]]),
    age = klinis_hc$agecontact,
    sex = klinis_hc$sex_contact,
    bmi = klinis_hc$BMI_contact,
    bcg = klinis_hc$BCG_scar,
    exp = klinis_hc$exp_score
  )
  
  # Remove incomplete cases
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 10) {
    result_lm <- rbind(result_lm, data.frame(
      metabolites = m,
      estimate = NA,
      std_error = NA,
      t_value = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # Fit linear model
  fit <- tryCatch(
    lm(outcome ~ metabolite + age + sex + bmi + bcg + exp, data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    result_lm <- rbind(result_lm, data.frame(
      metabolites = m,
      estimate = NA,
      std_error = NA,
      t_value = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  s <- summary(fit)$coefficients
  
  # Extract metabolite coefficient (row 2)
  result_lm <- rbind(result_lm, data.frame(
    metabolites = m,
    estimate = s[2, 1],
    std_error = s[2, 2],
    t_value = s[2, 3],
    p_value = s[2, 4],
    stringsAsFactors = FALSE
  ))
}

# Adjust p-values (Benjamini-Hochberg)
result_lm$p_adj <- p.adjust(result_lm$p_value, method = "BH")

# Join with metabolite names
result_lm <- dplyr::left_join(result_lm, name, by = "metabolites")

# Determine significance
SIG_P <- 0.05
result_lm$significant <- result_lm$p_value < SIG_P

# Color grouping
result_lm <- result_lm %>%
  dplyr::mutate(
    group = dplyr::case_when(
      is.na(p_value) ~ "Not tested",
      p_value >= SIG_P ~ "Not significant",
      estimate > 0 & p_value < SIG_P ~ "Positive association",
      estimate < 0 & p_value < SIG_P ~ "Negative association",
      TRUE ~ "Not significant"
    )
  )

# Save results
write.csv(result_lm,
          file.path(OUTPUT_DIR, "FigS2_linear_regression_IGRA_quantitative.csv"),
          row.names = FALSE)

# Summary statistics (all metabolites)
n_sig <- sum(result_lm$significant, na.rm = TRUE)
n_total <- sum(!is.na(result_lm$p_value))
cat("Significant associations (p < 0.05):", n_sig, "out of", n_total, "\n")
cat("Proportion:", round(n_sig / n_total * 100, 1), "%\n")

# =============================================================================
# ANALYSIS ON DIFFERENTIALLY ABUNDANT METABOLITES ONLY
# Manuscript states: "Almost all (94%) differentially abundant metabolites
# also significantly correlated with quantitative IGRA results at follow-up"
# =============================================================================

# Load differentially abundant metabolites from strict logistic regression
da_file <- file.path(OUTPUT_DIR, "3. significant_metabolites_subset.csv")

if (file.exists(da_file)) {
  da_metabolites <- read.csv(da_file)
  da_ids <- da_metabolites$metabolites
  
  cat("\n=== Differentially abundant metabolites ===\n")
  cat("Total DA metabolites:", length(da_ids), "\n")
  
  # Subset results to DA metabolites only
  result_lm_da <- result_lm %>% dplyr::filter(metabolites %in% da_ids)
  
  n_sig_da <- sum(result_lm_da$significant, na.rm = TRUE)
  n_total_da <- sum(!is.na(result_lm_da$p_value))
  
  cat("DA metabolites with significant correlation (p < 0.05):", n_sig_da, "out of", n_total_da, "\n")
  cat("Proportion:", round(n_sig_da / n_total_da * 100, 1), "%\n")
  
  # Save DA-only results
  write.csv(result_lm_da,
            file.path(OUTPUT_DIR, "FigS2_linear_regression_DA_only.csv"),
            row.names = FALSE)
} else {
  cat("\nNote: DA metabolite file not found:", da_file, "\n")
  cat("Run 01_logistic_regression_only.R first to generate the strict analysis results.\n")
}

# Volcano plot ------------------------------------------------------------
plot_data <- result_lm %>%
  dplyr::filter(!is.na(estimate), !is.na(p_value)) %>%
  dplyr::mutate(
    neg_log10_p = -log10(p_value)
  )

# Define colors (matching main manuscript volcano plot style)
# Positive correlation (estimate > 0): higher metabolite -> higher IFN-gamma -> RED (converter-like)
# Negative correlation (estimate < 0): higher metabolite -> lower IFN-gamma -> BLUE (persistently negative-like)
keyvals.colour <- ifelse(
  plot_data$p_value >= SIG_P, "black",
  ifelse(plot_data$p_value < SIG_P & plot_data$estimate > 0, "red",
         ifelse(plot_data$p_value < SIG_P & plot_data$estimate < 0, "blue",
                "black"))
)
keyvals.colour[is.na(keyvals.colour)] <- "black"
names(keyvals.colour)[keyvals.colour == "black"] <- "Not significant"
names(keyvals.colour)[keyvals.colour == "red"] <- "Positive correlation with\nTBAg-Nil IFN-gamma at follow-up"
names(keyvals.colour)[keyvals.colour == "blue"] <- "Negative correlation with\nTBAg-Nil IFN-gamma at follow-up"

# Create volcano plot (matching main manuscript styling)
plot_s2 <- EnhancedVolcano(
  plot_data,
  lab = plot_data$name,
  titleLabSize = 1,
  subtitleLabSize = 1,
  subtitle = "",
  x = "estimate",
  y = "p_value",
  legendPosition = "top",
  legendLabSize = 19,
  legendIconSize = 5,
  legendDropLevels = FALSE,
  ylab = bquote(~ italic(p-value)),
  xlim = c(-3, 3),
  ylim = c(0, 4),
  colCustom = keyvals.colour,
  title = "",
  pCutoff = SIG_P,
  FCcutoff = 0,
  pointSize = 5.0,
  labSize = 5.25,
  colAlpha = 0.3,
  gridlines.minor = FALSE,
  xlab = "Linear Regression Coefficient Estimate",
  drawConnectors = FALSE,
  caption = paste0("total = ", nrow(plot_data), " metabolites"),
  widthConnectors = 1
)

# Adjust y-axis to show p-value scale
plot_s2 <- plot_s2 +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 1.30103, 2, 3), labels = c(0, 0.1, 0.05, 0.01, 0.001))

plot_s2

# Save figure
ggsave(
  file.path(FIGURE_DIR, "FigureS2_IGRA_quantitative_correlation.png"),
  plot_s2,
  width = 10,
  height = 9,
  dpi = 600
)

ggsave(
  file.path(FIGURE_DIR, "FigureS2_IGRA_quantitative_correlation.svg"),
  plot_s2,
  width = 10,
  height = 9
)

cat("Figure S2 saved to:", FIGURE_DIR, "\n")
cat("Regression results saved to:", OUTPUT_DIR, "\n")
