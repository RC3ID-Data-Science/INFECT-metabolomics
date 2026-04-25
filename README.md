# INFECT Metabolomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19773171.svg)](https://doi.org/10.5281/zenodo.19773171)
[![GWAS Catalog](https://img.shields.io/badge/GWAS%20Catalog-GCST90824109-green)](https://www.ebi.ac.uk/gwas/studies/GCST90824109)

R scripts and analysis code for the manuscript: **"Circulating metabolites associated with protection against *Mycobacterium tuberculosis* infection and inhibition of mycobacterial growth"**

---

## Overview

This repository contains the R analysis scripts for a study investigating circulating metabolites and genetic variants associated with protection against *Mycobacterium tuberculosis* (Mtb) infection. The study combines:

1. **Untargeted metabolomics** of serum from Indonesian tuberculosis household contacts (INFECT cohort)
2. **Genetic association analysis** of metabolic pathway genes
3. **Validation** using ex vivo mycobacterial growth inhibition assay (MGIA) in healthy Dutch adults (300BCG cohort)

---

## Study Design

### Indonesian Cohort (INFECT)
- **Setting**: Bandung, Indonesia (2014-2018)
- **Participants**: 199 IGRA-negative household contacts of sputum smear-positive tuberculosis patients
- **Follow-up**: 14 weeks
- **Metabolomics**: 528 annotated serum metabolites measured by LC-MS/MS (after QC from 537 initially measured)
- **Genetics**: 4,407,262 genetic variants after QC (HumanOmniExpressExome-8 v1.0)

### Dutch Validation Cohort (300BCG)
- **Setting**: The Netherlands (low TB incidence)
- **Participants**: 42 healthy adults (before BCG vaccination)
- **Assay**: Mycobacterial Growth Inhibition Assay (MGIA) using PBMCs
- **Metabolomics**: 1,334 annotated plasma metabolites by flow-injection TOF-MS

---

## Repository Structure

```
INFECT-metabolomics/
├── scripts/
│   ├── 00_extract_data.R                    # Extract CSV files from supplementary Excel
│   ├── 00_load_data.R                       # Load data and create common objects for analysis
│   ├── config.R                             # Central configuration (paths, parameters)
│   ├── 01_data_preprocessing/
│   │   ├── 01_load_metabolomics_data.R      # Load raw LC-MS data, create compound IDs
│   │   └── 02_normalize_metabolomics.R      # MetaboAnalystR normalization, QC, imputation
│   ├── 02_clinical_data/
│   │   └── 01_prepare_clinical_covariates.R # Merge clinical variables, exposure scoring
│   ├── 03_differential_analysis/
│   │   ├── 01_logistic_regression_only.R    # Main analysis: strict logistic regression (BMI + covariates)
│   │   ├── 01_differential_metabolomics.R   # Main analysis: logistic regression (all covariates)
│   │   ├── 02_differential_with_BCGscar.R   # Sensitivity analysis: including BCG scar as covariate
│   │   ├── 03_differential_no_BCGscar.R     # Sensitivity analysis: excluding BCG scar covariate
│   │   └── Fig5_violin_plots.R              # Figure 5: MGIA violin plots + genotype associations
│   ├── 04_pathway_analysis/
│   │   ├── arachidonic_acid.R               # Arachidonic acid metabolism pathway visualization
│   │   ├── glutathione.R                    # Glutathione metabolism pathway visualization
│   │   ├── glycerophospholipid.R            # Glycerophospholipid metabolism pathway visualization
│   │   ├── glycine.R                        # Glycine/serine metabolism pathway visualization
│   │   ├── proline.R                        # Arginine/proline metabolism pathway visualization
│   │   └── tryptophan.R                     # Tryptophan metabolism pathway visualization
│   ├── 05_figures/
│   │   ├── Fig2_IGRA_status.R               # Figure 2: Correlation heatmap + volcano plots
│   │   └── Fig6_cytokine_correlation.R      # Figure 6: Cytokine correlation plots
│   ├── 06_300BCG_MGIA/
│   │   ├── 300BCG_correlation.R             # 300BCG metabolite-MGIA correlation analysis
│   │   └── 300BCG_meta_MGIA.R               # 300BCG meta-analysis
│   ├── 07_genetic_analysis/
│   │   ├── GWAS_cumulative_plots.R          # Cumulative GWAS distribution plots
│   │   └── GWAS_distribution_plots.R        # GWAS result visualizations (Fig 4A, 4C, 4D)
│   └── Supplementary/
│       ├── FigS2_IGRA_quantitative_correlation.R  # Figure S2: Linear regression vs quantitative IGRA (68 sig. metabolites)
│       ├── FigS4_mgia_heatmap.R                   # Figure S4: MGIA supplementary heatmaps
│       ├── FigS5_venn_diagram.R                   # Figure S5: Venn diagram (Dutch vs Indonesian metabolites)
│       └── FigS8_cfu_clinical.R                   # Figure S8: CFU clinical correlation
├── data/
│   ├── Metabolomics_INFECT_supplementary_data.xlsx  # Primary data source (all processed data)
│   ├── INFECT_characteristics.csv                   # Clinical data and covariates
│   ├── INFECT_metabolites.csv                       # Normalized metabolomics data
│   ├── INFECT_HMDB_metabolite_name.csv              # HMDB ID to metabolite name mapping
│   ├── INFECT_exact_ages.csv                        # Exact ages (manuscript uses these, not range midpoints)
│   ├── INFECT_log2fc_tbag.csv                       # Log2FC of TBAg IFNγ (for quantitative IGRA analysis)
│   ├── INFECT_log_regression_strict.csv             # Logistic regression results
│   ├── INFECT_KEGG_name_significant.csv             # KEGG pathway mappings
│   ├── INFECT_targeted_metabolic_gene.csv           # Targeted metabolic gene list
│   ├── INFECT_summary_metabolic_gene.csv            # GWAS summary statistics
│   ├── INFECT_PRODH_proline.csv                     # PRODH genotype-phenotype data
│   ├── INFECT_AFMID_tryptophan.csv                  # AFMID genotype-phenotype data
│   ├── 300BCG_metabolites.csv                       # Dutch validation cohort metabolomics
│   ├── 300BCG_metabolite_names.csv                  # Dutch cohort metabolite names
│   ├── 300BCG_MGIA.csv                              # MGIA assay results
│   ├── 300BCG_lin_regression.csv                    # Dutch cohort regression results
│   └── Data_dictionary.csv                          # Data dictionary for all variables
├── output/                                    # Analysis results (CSV tables)
├── figures/                                   # Generated figures (PNG, SVG)
└── README.md
```

**Notes**:
- Scripts in `03_differential_analysis/` labeled "with_BCGscar" and "no_BCGscar" are sensitivity analyses. The main manuscript presents results with all covariates (age, sex, BMI, BCG scar, exposure score) included in a single logistic regression model.
- The `INFECT_exact_ages.csv` and `INFECT_log2fc_tbag.csv` files are required for exact reproduction of manuscript results. The supplementary Excel contains age ranges for privacy and does not include the pre-computed Log2FC_TBonly variable.
- The repository also contains additional scripts (`06_300BCG_MGIA/`, some supplementary figures) that were used for exploratory analyses but are **not presented in the manuscript**.

---

## Key Findings

### Differential Metabolomics
- **64 metabolites** significantly associated with protection (strict logistic regression, p < 0.05)
- **68 metabolites** significantly correlated with quantitative IGRA IFNγ at follow-up (linear regression with Log2FC_TBonly, p < 0.05)
- Almost all (65.6%) differentially abundant metabolites also significantly correlated with quantitative IGRA results

### Metabolic Pathways Associated with Protection
1. **Arachidonic acid metabolism** (FDR = 0.005)
2. **Glutathione metabolism** (FDR = 0.017)
3. **Arginine and proline metabolism** (FDR = 0.017)
4. **Tryptophan metabolism** (FDR = 0.087)

### Key Metabolites (6 shared between cohorts)
| Metabolite | Pathway | Association with Protection |
|------------|---------|----------------------------|
| Prostaglandin F2α | Arachidonic acid metabolism | Lower levels protective |
| Leukotriene B4 | Arachidonic acid metabolism | Lower levels protective |
| 13-HODE | Linoleic acid metabolism | Lower levels protective |
| Tryptophan | Tryptophan metabolism | Lower levels protective |
| Glycine | Glycine/serine metabolism | Lower levels protective |
| Proline | Arginine/proline metabolism | Lower levels protective |

### Genetic Variants
- **rs3642** (near *AFMID*, tryptophan metabolism): P = 1.8 × 10⁻⁵
- **rs13056032** (near *PRODH*, arginine/proline metabolism): P = 9.7 × 10⁻⁵

---

## Requirements

### Software
- **R** >= 4.4.1 (tested with R 4.4.3)
- **RStudio** (recommended)

### R Packages

#### Data manipulation and visualization
```r
install.packages(c("tidyverse", "ggplot2", "ggpubr", "cowplot", 
                   "ComplexHeatmap", "circlize", "RColorBrewer",
                   "EnhancedVolcano", "tableone", "VennDiagram", "grid"))
```

#### Statistical analysis
```r
install.packages(c("stats", "MASS"))  # Base R packages
```

#### Pathway analysis
- **Cytoscape** >= 3.10 with Metscape 3.10 app
- **MetaboAnalyst** 5.0 (web-based)

#### Genetic analysis
```r
# GEMMA (command-line tool)
# Download from: https://github.com/genetics-statistics/GEMMA

# R packages for genetic analysis
install.packages(c("LDlinkR", "otargen"))
```

---

## Quick Start

### Step 1: Download Supplementary Data

Download the supplementary Excel file from [Figshare](https://doi.org/10.6084/m9.figshare.31769128) and place it in `data/`.

### Step 2: Extract CSV Files

```r
source("scripts/00_extract_data.R")
```

This creates the following CSV files in `data/`:
- `INFECT_characteristics.csv` - Clinical data and covariates
- `INFECT_metabolites.csv` - Normalized metabolomics data
- `INFECT_HMDB_metabolite_name.csv` - HMDB ID to metabolite name mapping
- `INFECT_log_regression_strict.csv` - Logistic regression results
- `INFECT_KEGG_name_significant.csv` - KEGG pathway mappings
- `INFECT_targeted_metabolic_gene.csv` - Targeted metabolic gene list
- `INFECT_summary_metabolic_gene.csv` - GWAS summary statistics
- `INFECT_PRODH_proline.csv` - PRODH genotype-phenotype data
- `INFECT_AFMID_tryptophan.csv` - AFMID genotype-phenotype data
- `300BCG_metabolites.csv` - Dutch validation cohort metabolomics
- `300BCG_metabolite_names.csv` - Dutch cohort metabolite names
- `300BCG_MGIA.csv` - MGIA assay results
- `300BCG_lin_regression.csv` - Dutch cohort regression results

### Step 3: Add Required Auxiliary Files

For **exact reproduction** of manuscript results, you need two additional files that are not in the supplementary Excel:

1. **`INFECT_exact_ages.csv`** - Exact ages for all participants (the Excel contains age ranges for privacy)
2. **`INFECT_log2fc_tbag.csv`** - Pre-computed Log2FC_TBonly values (log2 fold-change of TBAg IFNγ from baseline to week 14)

Place these files in `data/` before running analyses.

### Step 4: Run Analyses

```r
source("scripts/config.R")
source("scripts/00_load_data.R")

# Now run any analysis script
source("scripts/03_differential_analysis/01_logistic_regression_only.R")
source("scripts/Supplementary/FigS2_IGRA_quantitative_correlation.R")
```

---

## Data Availability

- **Metabolomics data (INFECT Indonesian cohort)**: Metabolomics Workbench (study ID pending)
- **Metabolomics data (300BCG Dutch cohort)**: [GitLab - bcg300](https://gitlab.com/xavier-lab-computation/public/bcg300)
- **GWAS summary statistics**: [NHGRI-EBI GWAS Catalog - GCST90824109](https://www.ebi.ac.uk/gwas/studies/GCST90824109)
- **Supplementary data**: [Figshare](https://doi.org/10.6084/m9.figshare.31769128) (DOI: [10.6084/m9.figshare.31769128](https://doi.org/10.6084/m9.figshare.31769128))

---

## Analysis Workflow

### 1. Data Preprocessing
- Log2-transformation of LC-MS peak intensities
- Imputation of values below limit of detection (half of minimum detected value)
- Removal of metabolites with >25% missing values or CV >30%

### 2. Differential Metabolomics
- **Logistic regression** with `glm()` (binomial family)
- **Covariates**: age, sex, BMI, BCG scar presence, exposure risk score
- **Uncorrected P < 0.05** for exploratory analysis
- **FDR correction** (Benjamini-Hochberg) for multiple testing
- **Sensitivity analyses**: with and without BCG scar as covariate

### 3. Quantitative IGRA Correlation (Figure S2)
- **Linear regression** with `lm()` 
- **Outcome**: Log2FC_TBonly (log2 fold-change of TBAg IFNγ, baseline → week 14)
- **Predictor**: each baseline metabolite
- **Covariates**: age, sex, BMI, BCG scar, exposure risk score
- **Subjects**: all 199 participants (not just strict subset)
- **Result**: 68 metabolites significantly correlated (p < 0.05)

### 4. Pathway Analysis
- **Metscape 3.10** (Cytoscape app) for network visualization
- **MetaboAnalyst 5.0** for KEGG overrepresentation analysis
- Background: all 537 measurable annotated metabolites

### 5. Genetic Association
- **GEMMA** for mixed-model GWAS (adjusts for relatedness)
- LD pruning with LDlinkR (r² < 0.1, East Asian reference)
- **Open Targets Platform** for eQTL lookup via `otargen`
- Bonferroni correction for filtered variants

### 6. MGIA Validation
- Linear regression between metabolite levels and log CFU
- No covariate adjustment (N=42, limited sample size)
- Uncorrected P < 0.1 for pathway analysis power
- Cross-cohort metabolite intersection

---

## Reproducibility Notes

### Age Data
The supplementary Excel contains age ranges (e.g., "10-14") for participant privacy. The manuscript analysis used **exact ages**. To reproduce exact results:
- Place `INFECT_exact_ages.csv` in `data/`
- `00_load_data.R` will automatically use exact ages when available
- Without this file, age range midpoints will be used (results may differ slightly)

### Log2FC_TBonly
The supplementary Excel does not contain the pre-computed Log2FC_TBonly variable (log2 fold-change of TBAg IFNγ). This variable is required for Figure S2 (quantitative IGRA correlation). To reproduce:
- Place `INFECT_log2fc_tbag.csv` in `data/`
- `00_load_data.R` will automatically load it
- Without this file, the quantitative IGRA analysis will not be available

### Key Reproduced Results
| Analysis | Expected Result |
|----------|----------------|
| Strict logistic regression (64 metabolites) | 64 significant metabolites (p < 0.05) |
| Linear regression vs Log2FC_TBonly (Fig S2) | 68 significant metabolites (p < 0.05) |

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{pediatamasetiabudiawanCirculatingMetabolitesAssociated2026,
  title = {Circulating Metabolites Associated with Protection against {{Mycobacterium}} Tuberculosis Infection and Inhibition of Mycobacterial Growth [in Press]},
  author = {Setiabudiawan, Todia Pediatama and Apriani, Lika and {Avila-Pacheco}, Julian and Ardiansyah, Edwin and Kumar, Vinod and Alisjahbana, Bachti and Verrall, Ayesha and DiNardo, Andrew and Indrati, Agnes and Mourits, Vera and {de Bree}, Charlotte and Moorlag, Simone and Clish, Clary and Joosten, Leo and {van Laarhoven}, Arjan and Meijgaarden, Krista and Netea, Mihai and Joosten, Simone and Koeken, Valerie and Hill, Philip and {van Crevel}, Reinout},
  year = 2026,
  journal = {Communications Biology},
  doi = {10.21203/rs.3.rs-6882340/v1},
  abstract = {A significant proportion of individuals who are heavily exposed to infectious tuberculosis patients do not acquire Mycobacterium tuberculosis (Mtb) infection, as detected by an interferon gamma release assay (IGRA). We examined circulating metabolite profiles and metabolic genotypes in 199 heavily exposed IGRA-negative tuberculosis household contacts in Indonesia. Based on differentially abundant metabolites, activity of several pathways including arachidonic acid, arginine and proline, glutathione, and tryptophan metabolism, correlated with a negative IGRA at three months. SNPs near PRODH (involved in arginine and proline metabolism) were associated with circulating proline concentrations and persistently negative IGRA results, while SNPs near AFMID (involved in tryptophan metabolism) were associated with IGRA conversion. For further validation, plasma metabolomic profiles were correlated with mycobacterial growth inhibition by peripheral blood mononuclear cells from individuals in a low-incidence setting. Lower circulating concentrations of six metabolites, including leukotriene B4, proline, glycine, and tryptophan correlated with better growth control in non-exposed healthy individuals, as well as with stronger protection against IGRA-conversion in the Indonesian tuberculosis household contacts. Collectively, these data support the notion that circulating metabolites may impact innate host defense against Mtb infection, and that metabolic interventions may prevent tuberculosis infection and disease.},
  keywords = {accepted,IGRA,innate immunity,metabolites,metabolomics,tuberculosis}
}
```

Related publications:
- Verrall AJ, et al. Early Clearance of Mycobacterium tuberculosis: The INFECT 
  Case Contact Cohort Study in Indonesia. J Infect Dis. 2020;221:1351-1360.
- Setiabudiawan TP, et al. Immune correlates of early clearance of Mycobacterium 
  tuberculosis among tuberculosis household contacts in Indonesia. Nat Commun. 
  2025;16:309.
- van Meijgaarden KE, et al. BCG vaccination-induced acquired control of 
  mycobacterial growth differs from growth control preexisting to BCG vaccination. 
  Nat Commun. 2024;15:114.

---

## Ethical Approval

- **Indonesian study**: Health Research Ethics Committee Universitas Padjadjaran 
  (14/UN6.C2.1.2/KEPK/PN/2014) and Southern Health and Disability Ethics 
  Committee New Zealand (13/STH/132)
- **Dutch study (300BCG)**: Arnhem-Nijmegen Medical Ethical Committee 
  (NL58553.091.16)

---

## Contact

For questions about the code or analysis, please open an issue on GitHub or contact:

- **Todia P. Setiabudiawan** - todia.setiabudiawan@radboudumc.nl

**Affiliation**:
- Radboud University Medical Center, Nijmegen, The Netherlands
- Research Center for Care and Control of Infectious Disease (RC3ID), Universitas Padjadjaran, Bandung, Indonesia 

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

The data and code are provided for academic research purposes. Please cite the 
original publication if you use any part of this work.

---

## Acknowledgements

This work was supported by:
- University of Otago and Mercy Hospital, Dunedin, New Zealand
- European Union's Seventh Framework Programme (FP7/2007-2013), grant 305279
- Royal Netherlands Academy of Arts and Sciences (09-PD-14)
- The Netherlands Organization for Scientific Research (VIDI grant 017.106.310)
- ERC Advanced Grant (833247) and Spinoza Grant
- Danish National Research Foundation (DNRF108)

We thank all study participants and field teams in Indonesia and the Netherlands.
