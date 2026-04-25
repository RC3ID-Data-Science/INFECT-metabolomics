# 00_load_data.R
# Loads extracted data from supplementary Excel and creates common objects
# Source this at the beginning of analysis scripts

# Ensure dplyr is available for pipe operations
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

cat("=== Loading INFECT Metabolomics Data ===\n")

# Determine script location for relative sourcing
script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
repo_dir <- dirname(script_dir)

# Source config
config_path <- file.path(repo_dir, "scripts", "config.R")
if (file.exists(config_path)) {
  source(config_path)
} else {
  # Fallback: assume config.R is in the same directory
  if (file.exists("config.R")) {
    source("config.R")
  } else if (file.exists("../config.R")) {
    source("../config.R")
  } else {
    stop("config.R not found. Please ensure scripts/config.R exists.")
  }
}

# Load normalized metabolomics data
if (!file.exists(NORMALIZED_METABOLOMICS)) {
  stop("Metabolomics data not found: ", NORMALIZED_METABOLOMICS, 
       "\nRun scripts/00_extract_data.R first.")
}
metabolites <- read.csv(NORMALIZED_METABOLOMICS, check.names = FALSE)
cat("Loaded metabolites:", nrow(metabolites), "x", ncol(metabolites), "\n")

if (!file.exists(CLINICAL_DATA)) {
  stop("Clinical data not found:", CLINICAL_DATA,
       "\nRun scripts/00_extract_data.R first.")
}
klinis <- read.csv(CLINICAL_DATA, check.names = FALSE)
cat("Loaded clinical data:", nrow(klinis), "x", ncol(klinis), "\n")

# Load exact ages if available (manuscript used exact ages, not range midpoints)
# The supplementary Excel contains age ranges for privacy; exact ages are
# provided separately to ensure regression results match the manuscript.
EXACT_AGES_FILE <- file.path(DATA_DIR, "INFECT_exact_ages.csv")
if (file.exists(EXACT_AGES_FILE)) {
  exact_ages <- read.csv(EXACT_AGES_FILE, check.names = FALSE)
  klinis <- dplyr::left_join(klinis, exact_ages, by = "IdInfect")
  if ("agecontact_exact" %in% colnames(klinis)) {
    klinis$agecontact <- klinis$agecontact_exact
    klinis$agecontact_exact <- NULL
    cat("Loaded exact ages from INFECT_exact_ages.csv\n")
  }
} else if ("agecontact" %in% colnames(klinis)) {
  # Fallback: convert agecontact ranges to approximate midpoints
  age_map <- c(
    "5-9" = 7, "10-14" = 12, "15-19" = 17, "20-24" = 22,
    "25-29" = 27, "30-34" = 32, "35-39" = 37, "40-44" = 42,
    "45-49" = 47, "50-54" = 52, "55-59" = 57, "60-64" = 62,
    "65-69" = 67, "70-74" = 72, "75-79" = 77, "80-84" = 82
  )
  klinis$agecontact <- as.numeric(age_map[as.character(klinis$agecontact)])
  cat("WARNING: Using age range midpoints. Results may differ from manuscript.\n")
  cat("        Place INFECT_exact_ages.csv in data/ for exact reproduction.\n")
}

# Load Log2FC_TBonly values if available (needed for linear regression with IGRA quantitative outcomes)
# The supplementary Excel does not contain this pre-computed variable; it is provided separately.
LOG2FC_FILE <- file.path(DATA_DIR, "INFECT_log2fc_tbag.csv")
if (file.exists(LOG2FC_FILE)) {
  log2fc_data <- read.csv(LOG2FC_FILE, check.names = FALSE)
  klinis <- dplyr::left_join(klinis, log2fc_data, by = "IdInfect")
  cat("Loaded Log2FC_TBonly from INFECT_log2fc_tbag.csv\n")
} else {
  cat("WARNING: INFECT_log2fc_tbag.csv not found. Linear regression with Log2FC_TBonly will not be available.\n")
}

# Remove QC columns if present
qc_cols <- grep("^(internal|NA\\.|redundant)", colnames(metabolites), value = TRUE)
if (length(qc_cols) > 0) {
  metabolites <- metabolites %>% dplyr::select(-all_of(qc_cols))
}

# Remove clinical columns from metabolites (they're in klinis)
# The extracted CSV includes strict_0.15_conversion and status for convenience,
# but we want the pure metabolite matrix
clinical_cols_in_met <- intersect(c("strict_0.15_conversion", "status"), colnames(metabolites))
for (col in clinical_cols_in_met) {
  metabolites[[col]] <- NULL
}

# Store subject IDs and filter clinical data
metabolites_subject <- metabolites$IdInfect
klinis <- klinis %>% dplyr::filter(IdInfect %in% metabolites_subject)

# Create merged dataset
merge <- dplyr::left_join(metabolites, klinis, by = "IdInfect")
merge <- merge %>% dplyr::arrange(desc(-IdInfect))
rownames(merge) <- merge[,1]

# Separate back
metabolites <- merge %>% dplyr::select(all_of(colnames(metabolites)))
klinis <- merge %>% dplyr::select(all_of(colnames(klinis)))

# Set rownames and remove IdInfect from metabolites
rownames(metabolites) <- metabolites[,1]
metabolites$IdInfect <- NULL

# Create metabolite name vector
nama.metabolites <- colnames(metabolites)

# Create data_plot for figure scripts
data_plot <- merge %>%
  dplyr::mutate(status = factor(status,
                                levels = c("Converter", "Persistently_uninfected"),
                                labels = c("IGRA-converters", "Persistently\nIGRA-negatives")))

# Ensure key columns are properly typed
klinis$IdInfect <- as.character(klinis$IdInfect)

# Load metabolite name mapping if available
if (file.exists(METABOLITE_NAME_MAP)) {
  metabolite_names <- read.csv(METABOLITE_NAME_MAP, check.names = FALSE)
  cat("Loaded metabolite name mapping:", nrow(metabolite_names), "entries\n")
} else {
  metabolite_names <- NULL
}

cat("=== Data loading complete ===\n")
cat("Objects created: merge (", nrow(merge), "x", ncol(merge), 
    "), metabolites (", nrow(metabolites), "x", ncol(metabolites),
    "), klinis (", nrow(klinis), "x", ncol(klinis), ")\n\n")
