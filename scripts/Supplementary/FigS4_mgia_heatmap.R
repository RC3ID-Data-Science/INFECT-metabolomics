# Supplementary Figure 4. Clinical variables are not associated with mycobacterial growth inhibition
#
# Panel A: Heatmap showing Spearman correlations between 12 clinical variables and CFU in the MGIA (n=42).
#          No significant associations were observed after FDR correction (all FDR > 0.05).
# Panel B: Heatmap of Spearman correlations between 1,334 plasma metabolites and 12 clinical variables.
#          Only significant correlations (FDR < 0.05) are shown in color (green: positive, pink: negative);
#          non-significant correlations appear grey.

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# =============================================================================
# PANEL A: Clinical variables vs CFU in MGIA
# =============================================================================

# Check if CFU data exists
if (exists("cfu_mgia") || exists("clinical_mgia")) {
  
  # Clinical variables to correlate with CFU
  clinical_vars_a <- c("Alcohol", "Age", "Hemoglobin", "Sleep", "Probiotic",
                       "Antibiotic", "Fruit_freq", "Vegetable_freq", "Smoking",
                       "BMI", "Stress", "Sex")
  
  # Use available variables
  available_vars_a <- clinical_vars_a[clinical_vars_a %in% names(clinical_mgia)]
  
  if (length(available_vars_a) > 0 && exists("cfu_mgia")) {
    cat("=== Panel A: Clinical variables vs CFU ===\n")
    cat("Clinical variables:", paste(available_vars_a, collapse = ", "), "\n")
    
    # Calculate Spearman correlations between clinical variables and CFU
    n_vars <- length(available_vars_a)
    cor_a <- numeric(n_vars)
    pval_a <- numeric(n_vars)
    
    for (j in 1:n_vars) {
      clin_vals <- as.numeric(clinical_mgia[[available_vars_a[j]]])
      cfu_vals <- as.numeric(cfu_mgia)
      
      complete_idx <- complete.cases(clin_vals, cfu_vals)
      
      if (sum(complete_idx) > 5) {
        test_result <- cor.test(clin_vals[complete_idx], cfu_vals[complete_idx],
                                method = "spearman", exact = FALSE)
        cor_a[j] <- test_result$estimate
        pval_a[j] <- test_result$p.value
      }
    }
    
    # FDR correction
    fdr_a <- p.adjust(pval_a, method = "BH")
    
    cat("Panel A correlations:\n")
    for (j in 1:n_vars) {
      cat(sprintf("  %s: rho=%.3f, p=%.4f, FDR=%.4f\n", 
                  available_vars_a[j], cor_a[j], pval_a[j], fdr_a[j]))
    }
    
    n_sig_a <- sum(fdr_a < 0.05, na.rm = TRUE)
    cat("Significant associations (FDR < 0.05):", n_sig_a, "\n")
    
    # Create single-row heatmap for Panel A
    cor_matrix_a <- matrix(cor_a, nrow = 1)
    colnames(cor_matrix_a) <- available_vars_a
    rownames(cor_matrix_a) <- "CFU"
    
    # Mask non-significant
    cor_matrix_a_masked <- cor_matrix_a
    cor_matrix_a_masked[, fdr_a >= 0.05] <- NA
    
    col_fun_a <- colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                            c("#D01C8B", "#F1B6DA", "#F7F7F7", "#B8E186", "#4DAC26"))
    
    ht_a <- Heatmap(cor_matrix_a_masked,
                    name = "Spearman rho",
                    col = col_fun_a,
                    na_col = "gray90",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    column_names_rot = 45,
                    column_names_side = "bottom",
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(
                      title = "Spearman rho",
                      title_gp = gpar(fontsize = 10),
                      labels_gp = gpar(fontsize = 8),
                      legend_height = unit(3, "cm")
                    ))
  }
}

# =============================================================================
# PANEL B: 1334 metabolites vs clinical variables
# =============================================================================

# Define clinical variables to use
clinical_vars_b <- c("Alcohol", "Age", "Hemoglobin", "Sleep", "Probiotic",
                     "Antibiotic", "Fruit_freq", "Vegetable_freq", "Smoking",
                     "BMI", "Stress", "Sex")

# Check if all clinical variables exist in clinical_mgia
available_vars_b <- clinical_vars_b[clinical_vars_b %in% names(clinical_mgia)]
cat("\n=== Panel B: Metabolites vs clinical variables ===\n")
cat("Using clinical variables:", paste(available_vars_b, collapse = ", "), "\n")

# Get metabolite data
metabolite_matrix <- as.matrix(metabolite_data_mgia)
n_metabolites <- ncol(metabolite_matrix)
n_participants <- nrow(metabolite_matrix)

cat("Number of metabolites:", n_metabolites, "\n")
cat("Number of participants:", n_participants, "\n")

# Get clinical data
clinical_matrix <- as.matrix(clinical_mgia[, available_vars_b])

# Calculate Spearman correlations and p-values
cor_matrix_b <- matrix(NA, nrow = n_metabolites, ncol = length(available_vars_b))
pval_matrix_b <- matrix(NA, nrow = n_metabolites, ncol = length(available_vars_b))

rownames(cor_matrix_b) <- colnames(metabolite_matrix)
colnames(cor_matrix_b) <- available_vars_b
rownames(pval_matrix_b) <- colnames(metabolite_matrix)
colnames(pval_matrix_b) <- available_vars_b

cat("Calculating correlations...\n")
for (i in 1:n_metabolites) {
  for (j in 1:length(available_vars_b)) {
    met_vals <- as.numeric(metabolite_matrix[, i])
    clin_vals <- as.numeric(clinical_matrix[, j])

    # Remove NAs
    complete_idx <- complete.cases(met_vals, clin_vals)

    if (sum(complete_idx) > 5) {
      test_result <- cor.test(met_vals[complete_idx], clin_vals[complete_idx],
                              method = "spearman", exact = FALSE)
      cor_matrix_b[i, j] <- test_result$estimate
      pval_matrix_b[i, j] <- test_result$p.value
    }
  }
}

# Apply FDR correction
cat("Applying FDR correction...\n")
pval_vector_b <- as.vector(pval_matrix_b)
fdr_vector_b <- p.adjust(pval_vector_b, method = "BH")
fdr_matrix_b <- matrix(fdr_vector_b, nrow = n_metabolites, ncol = length(available_vars_b))

# Create masked correlation matrix (only show significant)
cor_matrix_b_masked <- cor_matrix_b
cor_matrix_b_masked[fdr_matrix_b >= 0.05] <- NA

# Count significant correlations
n_significant_b <- sum(fdr_matrix_b < 0.05, na.rm = TRUE)
cat("Number of significant correlations (FDR < 0.05):", n_significant_b, "\n")

# Create color function (green for positive, pink for negative)
col_fun_b <- colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                        c("#D01C8B", "#F1B6DA", "#F7F7F7", "#B8E186", "#4DAC26"))

# Hierarchical clustering for rows (metabolites)
cat("Performing hierarchical clustering...\n")
cor_for_cluster <- cor_matrix_b
cor_for_cluster[is.na(cor_for_cluster)] <- 0
row_dist <- dist(cor_for_cluster)
row_clust <- hclust(row_dist, method = "complete")

ht_b <- Heatmap(cor_matrix_b_masked,
                name = "Spearman rho",
                col = col_fun_b,
                na_col = "gray90",
                cluster_rows = row_clust,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_rot = 45,
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(
                  title = "Spearman rho",
                  title_gp = gpar(fontsize = 10),
                  labels_gp = gpar(fontsize = 8),
                  legend_height = unit(4, "cm")
                ),
                column_title = paste0("Correlation of ", format(n_metabolites, big.mark = ","),
                                      " metabolites with clinical variables (n=", n_participants, ")"),
                column_title_gp = gpar(fontsize = 11, fontface = "bold"))

# =============================================================================
# SAVE COMBINED FIGURE
# =============================================================================

png("Fig S4. MGIA clinical heatmaps.png",
    width = 12, height = 10, units = "in", res = 300)

if (exists("ht_a")) {
  # Draw both panels
  draw(ht_a %v% ht_b,
       column_title = "Supplementary Figure 4",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
} else {
  draw(ht_b,
       column_title = "Supplementary Figure 4",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
}

dev.off()

cat("\nFigure saved as 'Fig S4. MGIA clinical heatmaps.png'\n")
cat("Done!\n")
