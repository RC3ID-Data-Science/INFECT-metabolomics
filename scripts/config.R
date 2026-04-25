# =============================================================================
# Configuration File for INFECT Metabolomics Analysis
# =============================================================================
# 
# PURPOSE: Centralize all file paths and parameters used across the analysis
#          pipeline. Users should modify paths below to match their local
#          environment.
#
# USAGE:   Source this file at the beginning of each analysis script:
#          source("config.R")
#
# AUTHOR:  Todia P. Setiabudiawan
# DATE:    2026-04-25
# =============================================================================

# =============================================================================
# 1. BASE DIRECTORIES
# =============================================================================

# Main data directory - where all input data files are stored
# Default: Points to the INFECT metabolomics analysis folder
# UPDATE THIS to your local path:
# DATA_DIR <- "C:/Users/todia/Google Drive/1. PhD/TP0001_INFECT/Metabolomics/1. todia analisis/"

# Repository directory (where this code lives)
REPO_DIR <- "D:/Users/todia/OneDrive - Radboudumc/MK TS RvC Projects - INFECT metabolomics/"

# Extracted data directory (CSV files from supplementary Excel)
DATA_DIR <- file.path(REPO_DIR, "data/")

# Raw data directory (original files from Google Drive)
RAW_DATA_DIR <- "D:/Users/todia/Google Drive/1. PhD/TP0001_INFECT/Metabolomics/1. todia analisis/"

# Output directory for results, figures, and tables
OUTPUT_DIR <- file.path(REPO_DIR, "output/")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Supplementary Excel file (primary data source, from Figshare)
SUPPLEMENTARY_EXCEL <- file.path(DATA_DIR, "Metabolomics_INFECT_supplementary_data.xlsx")

# =============================================================================
# 2. METABOLOMICS DATA FILES
# =============================================================================

# Raw metabolomics results from 4 LC-MS methods (in raw data directory)
HILIC_POS_RAW   <- paste0(RAW_DATA_DIR, "1. 20_0401_INFECT_HILIC-pos_metabolomics_results.xlsx")
HILIC_NEG_RAW   <- paste0(RAW_DATA_DIR, "1. 21_0305_INFECT_HILIC-neg_metabolomics_results.xlsx")
C8_POS_RAW      <- paste0(RAW_DATA_DIR, "1. 20_0401_INFECT_C8-pos_metabolomics_results.xlsx")
C18_NEG_RAW     <- paste0(RAW_DATA_DIR, "1. 20_0401_INFECT_C18-neg_metabolomics_results.xlsx")

# Normalized metabolomics data (output from preprocessing)
# Extracted from supplementary Excel: Sheet "INFECT metabolites"
NORMALIZED_METABOLOMICS <- file.path(DATA_DIR, "INFECT_metabolites.csv")

# Metabolite name mapping
# Extracted from supplementary Excel: Sheet "INFECT HMDB-metabolite name"
METABOLITE_NAME_MAP <- file.path(DATA_DIR, "INFECT_HMDB_metabolite_name.csv")

# MetaboAnalyst covariate analysis result
METABOANALYST_RESULT <- "C:/Users/todia/Downloads/covariate_result_full_metabo.csv"

# =============================================================================
# 3. CLINICAL DATA FILES
# =============================================================================

# Main clinical dataset with all covariates
# Extracted from supplementary Excel: Sheet "INFECT characteristics"
CLINICAL_DATA <- file.path(DATA_DIR, "INFECT_characteristics.csv")

# Outcome assignments (converter vs persistently uninfected)
# NOTE: Now included in CLINICAL_DATA (strict_0.15_conversion column)
# OUTCOME_FILE <- paste0(DATA_DIR, "2. assigned outcome.xlsx")

# Exposure risk score
# NOTE: Now included in CLINICAL_DATA (exp_score column)
# EXPOSURE_SCORE <- paste0(DATA_DIR, "2. exposure score.xlsx")

# IGRA results
# NOTE: Now included in CLINICAL_DATA (NIL1, TBAg.NIL1, NIL1_14, TBAg.NIL1_14 columns)
# IGRA_BASELINE <- paste0(DATA_DIR, "2. IGRA Baseline.xlsx")
# IGRA_FOLLOWUP <- paste0(DATA_DIR, "2. IGRA 14 weeks.xlsx")

# Additional clinical variables (in raw data directory)
CONTACT_CRF     <- paste0(RAW_DATA_DIR, "2. F9 contact baseline CRF.xlsx")
INDEX_CRF       <- paste0(RAW_DATA_DIR, "2. F3 CRF index case.xlsx")
INDEX_XRAY      <- paste0(RAW_DATA_DIR, "2. F4 index xray findings.xlsx")
MICROBIOLOGY    <- paste0(RAW_DATA_DIR, "2. microbiology results bandung.xlsx")
VITAMIN_D       <- paste0(RAW_DATA_DIR, "2. Vitamin D.xlsx")
IRON            <- paste0(RAW_DATA_DIR, "2. serum iron and inflam_2.xlsx")
MICRONUTRIENTS  <- paste0(RAW_DATA_DIR, "2. micronutrients.xlsx")
WBSA            <- paste0(RAW_DATA_DIR, "2. wbsa wholeblood baseline.xlsx")
WGS             <- paste0(RAW_DATA_DIR, "2. wgs lineages index and exclusions.xlsx")
LINKAGE         <- paste0(RAW_DATA_DIR, "2. linkage index contact.xlsx")
HEMATOLOGY      <- paste0(RAW_DATA_DIR, "2. haematology.xlsx")
SOCIOENV        <- paste0(RAW_DATA_DIR, "2. socioenvironment.xlsx")

# Eligible cases and contacts (full cohort)
ELIGIBLE_COHORT <- paste0(RAW_DATA_DIR, "2. Eligible cases and contacts n1347 - repaired.xlsx")

# =============================================================================
# 4. GENETIC DATA FILES
# =============================================================================

# Directory containing GWAS genotype data (imputed, dosage format)
GENOTYPE_DIR <- "D:/Download/General/genotype/"

# Chromosome files (dosage format from GEMMA/IMPUTE2 pipeline)
CHR17_DOSAGE <- paste0(GENOTYPE_DIR, "chr17_TBM_INFECT_imputed_EAS_IND.finalQC.dosage.rds")
CHR22_DOSAGE <- paste0(GENOTYPE_DIR, "chr22_TBM_INFECT_imputed_EAS_IND.finalQC.dosage.rds")

# =============================================================================
# 5. CYTOKINE / FUNCTIONAL DATA
# =============================================================================

# Ex vivo whole blood stimulation assay (WBSA) - batch-corrected cytokines
WBSA_CORRECTED    <- "D:/Users/todia/Google Drive/1. PhD/TP0001_INFECT/Antibody/INFECT_antibody/INFECt_ex_vivo_cytokines_corrected.rds"
WBSA_IFN_CORRECTED <- "D:/Users/todia/Google Drive/1. PhD/TP0001_INFECT/Antibody/INFECT_antibody/INFECt_ex_vivo_cytokines_corrected_with_IFN.rds"

# =============================================================================
# 6. ANALYSIS PARAMETERS
# =============================================================================

# IGRA cut-off definitions (strict cut-offs used in manuscript)
IGRA_BASELINE_CUTOFF <- 0.15    # IU/mL - baseline negative threshold
IGRA_FOLLOWUP_CUTOFF <- 0.7     # IU/mL - follow-up conversion threshold
IGRA_MANUFACTURER_CUTOFF <- 0.35 # IU/mL - standard QFT-GIT cut-off

# Missing value threshold for metabolite filtering
MISSING_VALUE_THRESHOLD <- 0.15  # Remove metabolites with >15% missing values

# Coefficient of variation threshold
CV_THRESHOLD <- 30               # Remove metabolites with CV > 30%

# Statistical significance thresholds
P_VALUE_THRESHOLD <- 0.05        # Uncorrected p-value for exploratory analysis
FDR_THRESHOLD <- 0.05            # Benjamini-Hochberg FDR threshold

# Mixed model parameters
MIXED_MODEL_OPTIMIZER <- "bobyqa"  # Optimization algorithm for glmer

# =============================================================================
# 7. FIGURE OUTPUT PATHS
# =============================================================================

# Directory for manuscript figures
# Default: Local figures directory (creates if doesn't exist)
FIGURE_DIR <- file.path(REPO_DIR, "figures/")

# Alternative: Original manuscript figures folder
# FIGURE_DIR <- "C:/Users/todia/Google Drive/1. PhD/TP0001_INFECT/Metabolomics/Paper/"

# Directory for supplementary/revision figures
REVISION_DIR <- paste0(FIGURE_DIR, "revision/")

# =============================================================================
# 8. METABOLITE HMDB IDs (Key metabolites from the manuscript)
# =============================================================================

# The 7 metabolites featured in Figure 5D
KEY_METABOLITES <- c(
  "HMDB0000162" = "Proline",
  "HMDB0001085" = "Leukotriene B4",
  "HMDB0000929" = "Tryptophan",
  "HMDB0000123" = "Glycine",
  "HMDB0004667" = "13-HODE",
  "HMDB0001139" = "Prostaglandin F2alpha",
  "HMDB0000197" = "Indoleacetic acid"
)

# =============================================================================
# 9. GENETIC VARIANTS OF INTEREST
# =============================================================================

# rs13056032 - PRODH (Arginine/Proline metabolism)
PRODH_SNP <- "rs13056032"
PRODH_CHR <- "chr22"
PRODH_POS <- "22:18938193:T:C"

# rs3642 - AFMID (Tryptophan metabolism)
AFMID_SNP <- "rs3642"
AFMID_CHR <- "chr17"
AFMID_POS <- "17:76162178:G:A"

# =============================================================================
# 10. UTILITY FUNCTIONS
# =============================================================================

# Function to create output directories if they don't exist
create_output_dirs <- function() {
  dirs <- c(OUTPUT_DIR, FIGURE_DIR, REVISION_DIR)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# Function to check if all required files exist
check_required_files <- function() {
  required <- c(
    NORMALIZED_METABOLOMICS,
    CLINICAL_DATA,
    METABOLITE_NAME_MAP
  )
  
  missing <- required[!file.exists(required)]
  
  if (length(missing) > 0) {
    warning("Missing required files:\n", paste(missing, collapse = "\n"))
    return(FALSE)
  }
  
  message("All required files found.")
  return(TRUE)
}

# Helper function to read data directly from the supplementary Excel file
# This avoids extracting CSVs if you prefer working with the Excel directly
read_supplementary <- function(sheet_name) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl", repos = "https://cloud.r-project.org/")
  }
  library(readxl)
  
  if (!file.exists(SUPPLEMENTARY_EXCEL)) {
    stop("Supplementary Excel file not found: ", SUPPLEMENTARY_EXCEL)
  }
  
  readxl::read_excel(SUPPLEMENTARY_EXCEL, sheet = sheet_name)
}

# Convenience wrappers for common sheets
read_clinical_data <- function() read_supplementary("INFECT characteristics")
read_metabolites <- function() read_supplementary("INFECT metabolites")
read_metabolite_names <- function() read_supplementary("INFECT HMDB-metabolite name")
read_regression_results <- function() read_supplementary("INFECT log.regression strict")
read_genetic_results <- function() read_supplementary("INFECT summary metabolic gene")

# =============================================================================
# 11. SESSION INFO
# =============================================================================

cat("\n=== INFECT Metabolomics Configuration ===\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Figure directory:", FIGURE_DIR, "\n")
cat("=====================================\n\n")

# =============================================================================
# NOTES FOR REPRODUCIBILITY
# =============================================================================
#
# 1. All paths use forward slashes (/) for cross-platform compatibility.
#    R automatically handles this on Windows.
#
# 2. DATA SOURCE: The primary data source is the supplementary Excel file
#    (Metabolomics_INFECT_supplementary_data.xlsx) which contains all
#    processed data from the manuscript. This file is available on Figshare
#    with DOI: 10.6084/m9.figshare.31769128.
#
# 3. SETUP INSTRUCTIONS:
#    a) Place the supplementary Excel file in data/ directory
#    b) Run scripts/00_extract_data.R to extract CSV files
#    c) Source this config.R file to load all paths
#    d) Run analysis scripts
#
#    Alternatively, use read_supplementary("sheet_name") to read directly
#    from the Excel file without extracting CSVs.
#
# 4. The raw metabolomics files (HILIC_pos/neg, C8_pos, C18_neg) are only
#    needed for the preprocessing scripts (01_data_preprocessing/).
#    The main analysis scripts only require NORMALIZED_METABOLOMICS.
#
# 5. Genetic data (CHR17_DOSAGE, CHR22_DOSAGE) and cytokine data (WBSA_*) 
#    are only needed for the genotype-metabolite and cytokine analyses.
#    These are NOT included in the supplementary Excel due to size limits.
#
# 6. If you don't have certain data files, the corresponding analysis
#    scripts will fail at the file loading step. Skip those scripts or
#    comment out the relevant sections.
#
# =============================================================================
