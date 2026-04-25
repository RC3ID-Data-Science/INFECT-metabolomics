# =============================================================================
# Supplementary Figure 5. Overlapping significant metabolites in Indonesian and Dutch dataset
# INFECT Metabolomics Study - Dutch vs Indonesian Cohort Comparison
# =============================================================================
#
# Panel A: Venn diagram showing overlap of all detected metabolites between cohorts
# Panel B: Venn diagram showing overlap of significant metabolites within shared set
# Panel C: Table of 7 overlapping metabolites with statistics

library(VennDiagram)
library(grid)
library(gridExtra)
library(dplyr)

# =============================================================================
# Extract metabolite lists
# =============================================================================

dutch_all <- correlation$HMDB_ID1
indonesian_all <- all_INFECT_metabolites$metabolites
dutch_sig <- results_name_sig$HMDB_ID1
indonesian_sig <- metabolite_list_INFECT$metabolites

# Calculate shared metabolites
shared_all <- intersect(dutch_all, indonesian_all)

# Significant within shared
dutch_sig_shared <- intersect(dutch_sig, shared_all)
indonesian_sig_shared <- intersect(indonesian_sig, shared_all)
both_sig_shared <- intersect(dutch_sig_shared, indonesian_sig_shared)

# =============================================================================
# Panel A: All Detected Metabolites
# =============================================================================

venn_A <- draw.pairwise.venn(
  area1 = length(dutch_all),
  area2 = length(indonesian_all),
  cross.area = length(shared_all),
  category = c("", ""),
  fill = c("#3498DB", "#E74C3C"),
  alpha = 0.65,
  col = c("#2980B9", "#C0392B"),
  lwd = 3,
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0,
  print.mode = "raw",
  margin = 0.05
)

# =============================================================================
# Panel B: Significant Metabolites within 125 Shared
# =============================================================================

venn_B <- draw.pairwise.venn(
  area1 = length(dutch_sig_shared),
  area2 = length(indonesian_sig_shared),
  cross.area = length(both_sig_shared),
  category = c("", ""),
  fill = c("#3498DB", "#E74C3C"),
  alpha = 0.65,
  col = c("#2980B9", "#C0392B"),
  lwd = 3,
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0,
  print.mode = "raw",
  margin = 0.05
)

# =============================================================================
# Panel C: Table of 7 overlapping metabolites
# =============================================================================

# Extract statistics for the 7 overlapping metabolites
# These should be available from the analysis results

if (exists("results_name_sig") && exists("logres")) {
  
  # Get Dutch MGIA statistics
  dutch_stats <- results_name_sig %>%
    filter(HMDB_ID1 %in% both_sig_shared) %>%
    select(HMDB_ID1, Metabolite_Name, Estimate, P_value) %>%
    rename(
      HMDB_ID = HMDB_ID1,
      Metabolite = Metabolite_Name,
      Dutch_Estimate = Estimate,
      Dutch_P = P_value
    )
  
  # Get Indonesian IGRA statistics
  indo_stats <- logres %>%
    filter(metabolites %in% both_sig_shared) %>%
    select(metabolites, name, log2OR, Pr...z..) %>%
    rename(
      HMDB_ID = metabolites,
      Metabolite = name,
      Indo_log2OR = log2OR,
      Indo_P = Pr...z..
    )
  
  # Combine
  table_c <- full_join(dutch_stats, indo_stats, by = "HMDB_ID") %>%
    mutate(Metabolite = coalesce(Metabolite.x, Metabolite.y)) %>%
    select(HMDB_ID, Metabolite, Indo_log2OR, Indo_P, Dutch_Estimate, Dutch_P) %>%
    arrange(HMDB_ID)
  
} else {
  # Fallback: create table structure if data objects not available
  cat("Warning: Data objects for Panel C not found. Creating empty table structure.\n")
  cat("Expected: results_name_sig (Dutch MGIA results) and logres (Indonesian logistic regression results)\n")
  
  table_c <- data.frame(
    HMDB_ID = character(),
    Metabolite = character(),
    Indo_log2OR = numeric(),
    Indo_P = numeric(),
    Dutch_Estimate = numeric(),
    Dutch_P = numeric(),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Create Combined Figure (Panels A + B)
# =============================================================================

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(0.08, 0.92), "npc"),
                                           widths = unit(c(0.5, 0.5), "npc"))))

# Panel A title
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()

# Panel B title
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()

# Panel A Venn
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(venn_A)
popViewport()

# Panel B Venn
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(venn_B)
popViewport()

# =============================================================================
# Save Figures
# =============================================================================

# Save combined figure (Panels A + B)
png("Fig S5. venn_combined.png", width = 10, height = 5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(0.08, 0.92), "npc"),
                                           widths = unit(c(0.5, 0.5), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(venn_A)
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(venn_B)
popViewport()
dev.off()

# Save as PDF
pdf("Fig S5. venn_combined.pdf", width = 10, height = 5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(0.08, 0.92), "npc"),
                                           widths = unit(c(0.5, 0.5), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("A", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.text("B", x = 0.08, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(venn_A)
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(venn_B)
popViewport()
dev.off()

# Save Panel C table as CSV
if (nrow(table_c) > 0) {
  write.csv(table_c, "Fig S5. overlapping_metabolites_table.csv", row.names = FALSE)
  cat("\nPanel C table saved to: Fig S5. overlapping_metabolites_table.csv\n")
}

cat("\n=== Figures saved! ===\n")
cat("- Fig S5. venn_combined.png\n")
cat("- Fig S5. venn_combined.pdf\n")
if (nrow(table_c) > 0) {
  cat("- Fig S5. overlapping_metabolites_table.csv\n")
}

# =============================================================================
# Summary
# =============================================================================

cat("\n=== SUMMARY ===\n")
cat("Panel A - All Metabolites:\n")
cat("  Dutch (MGIA):", length(dutch_all), "\n")
cat("  Indonesian (IGRA):", length(indonesian_all), "\n")
cat("  Overlap:", length(shared_all), "\n")

cat("\nPanel B - Significant within", length(shared_all), "shared:\n")
cat("  Dutch significant:", length(dutch_sig_shared), "\n")
cat("  Indonesian significant:", length(indonesian_sig_shared), "\n")
cat("  Both significant:", length(both_sig_shared), "\n")

cat("\nPanel C - Overlapping metabolites table:\n")
cat("  Rows:", nrow(table_c), "\n")
