#!/usr/bin/env Rscript
# 00_extract_data.R
# Extracts all data sheets from the supplementary Excel file into CSV format
# Run this once before executing analysis scripts

cat("=== INFECT Metabolomics Data Extraction ===\n\n")

# Check for required package
if (!requireNamespace("readxl", quietly = TRUE)) {
  cat("Installing readxl package...\n")
  install.packages("readxl", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("tools", quietly = TRUE)) {
  install.packages("tools", repos = "https://cloud.r-project.org/")
}

library(readxl)

# Get the directory where this script is located
script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
repo_dir <- dirname(script_dir)
data_dir <- file.path(repo_dir, "data")
excel_file <- file.path(data_dir, "Metabolomics_INFECT_supplementary_data.xlsx")

# Verify Excel file exists
if (!file.exists(excel_file)) {
  stop("Excel file not found: ", excel_file, 
       "\nPlease ensure 'Metabolomics_INFECT_supplementary_data.xlsx' is in the data/ directory.")
}

cat("Reading Excel file:", basename(excel_file), "\n")
sheets <- excel_sheets(excel_file)
cat("Found", length(sheets), "sheets:\n")
print(sheets)
cat("\n")

# Define sheet to filename mapping
# These filenames match what the analysis scripts expect
sheet_mapping <- list(
  "INFECT characteristics" = "INFECT_characteristics.csv",
  "INFECT metabolites" = "INFECT_metabolites.csv",
  "INFECT HMDB-metabolite name" = "INFECT_HMDB_metabolite_name.csv",
  "INFECT log.regression strict" = "INFECT_log_regression_strict.csv",
  "INFECT KEGG name significant" = "INFECT_KEGG_name_significant.csv",
  "INFECT targeted metabolic gene" = "INFECT_targeted_metabolic_gene.csv",
  "INFECT summary metabolic gene" = "INFECT_summary_metabolic_gene.csv",
  "INFECT PRODH - Proline" = "INFECT_PRODH_proline.csv",
  "INFECT AFMID - Tryptophan" = "INFECT_AFMID_tryptophan.csv",
  "300BCG metabolites" = "300BCG_metabolites.csv",
  "300BCG metabolite names" = "300BCG_metabolite_names.csv",
  "300BCG MGIA" = "300BCG_MGIA.csv",
  "300BCG lin.regression" = "300BCG_lin_regression.csv"
)

# Extract each sheet
extracted_files <- c()
errors <- c()

for (sheet_name in names(sheet_mapping)) {
  if (!(sheet_name %in% sheets)) {
    msg <- paste0("WARNING: Sheet '", sheet_name, "' not found in Excel file. Skipping.")
    cat(msg, "\n")
    errors <- c(errors, msg)
    next
  }
  
  output_file <- file.path(data_dir, sheet_mapping[[sheet_name]])
  
  cat("Extracting:", sheet_name, "->", basename(output_file), "\n")
  
  tryCatch({
    df <- read_excel(excel_file, sheet = sheet_name)
    write.csv(df, output_file, row.names = FALSE)
    extracted_files <- c(extracted_files, output_file)
    cat("  ✓ Saved", nrow(df), "rows x", ncol(df), "columns\n")
  }, error = function(e) {
    msg <- paste0("ERROR extracting '", sheet_name, "': ", e$message)
    cat("  ", msg, "\n")
    errors <- c(errors, msg)
  })
}

# Also extract data dictionary (optional, for reference)
if ("Data dictionary" %in% sheets) {
  output_file <- file.path(data_dir, "Data_dictionary.csv")
  cat("\nExtracting: Data dictionary ->", basename(output_file), "\n")
  tryCatch({
    df <- read_excel(excel_file, sheet = "Data dictionary")
    write.csv(df, output_file, row.names = FALSE)
    extracted_files <- c(extracted_files, output_file)
    cat("  ✓ Saved", nrow(df), "rows x", ncol(df), "columns\n")
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
}

cat("\n=== Extraction Complete ===\n")
cat("Files extracted to:", data_dir, "\n")
cat("Total files:", length(extracted_files), "\n")

if (length(errors) > 0) {
  cat("\nWARNINGS/ERRORS (", length(errors), "):\n")
  for (err in errors) {
    cat("  -", err, "\n")
  }
}

cat("\nNext steps:\n")
cat("  1. Source scripts/config.R to load data paths\n")
cat("  2. Run analysis scripts in order\n")
