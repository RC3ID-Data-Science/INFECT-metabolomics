
# load package ------------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggpubr)
library(data.table)
library(conflicted)
library(gtsummary)
conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::first)
library(easystats)
library(factoextra)
library(ggprism)
library(FactoMineR)
library(broom)
`%out%` <- negate(`%in%`)
library(janitor)
library(purrr)
library(data.table)
library(ComplexHeatmap)
library(rstatix)
library(extrafont)

theme_set(theme_pubr(base_family = "Helvetica Neue"))

# load data ---------------------------------------------------------------
all_300BCG = read_excel("190321 300BCG Dataset wide.xlsx")
meta_cyto <- load("17.02.2021 Metabolomics 300BCG dataset.RData")
#meta_cyto <- as.data.frame(read.csv("300BCG_17.02.2021_Metabolomics_300BCG_dataset.csv"))
all_INFECT_metabolites <- read.csv("3. result logres metabolites_subset.csv") %>% select(metabolites,name)
metabolite_list_INFECT <- read.csv("3. significant_metabolites_subset.csv") %>% select(metabolites,name)
metabolite_list_300BCG <- read_excel("300BCG_08.08.2023_Metabolites_HMDBIDS.xlsx")
MGIA_300BCG <- read_excel("300BCG_MGIA.xlsx") %>%
  rename(PatientID = ID)

MGIA_300BCG |> 
  select(age, gender) |> 
  gtsummary::tbl_summary()
?tbl_summary

metabolite_list_300BCG <- metabolite_list_300BCG %>%
  select(ionIdx,
         score,
         `label (bona fide)`,
         HMDB_ID1)

subset_metabolite_list_300BCG <- metabolite_list_300BCG %>%
  filter(score == 100)

metabolite_select <- colnames(mydata[53:1425])
#metabolite_select <- c(as.character(subset_metabolite_list_300BCG$ionIdx))
duplicated(subset_metabolite_list_300BCG$ionIdx)

klinis_meta = left_join(MGIA_300BCG,mydata,by="PatientID")
join_meta_MGIA <- left_join(MGIA_300BCG,mydata,by="PatientID")
join_meta_MGIA <- join_meta_MGIA %>%
  select(PatientID,
         `pre vacc`,
         metabolite_select) %>%
  column_to_rownames("PatientID")



# distribution ------------------------------------------------------------

dist = MGIA_300BCG |> 
  select(`pre vacc`) |> 
  rename(MGIA = `pre vacc`)

dist$Classification = if_else(dist$MGIA > 1.26, "No control","Controlling")

ggpubr::gghistogram(dist, x = "MGIA", bins = 20,
                    fill = "Classification", palette = c("blue","red"), alpha=0.3, legend ="none") +
  geom_vline(xintercept = 2.38, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = 1.28, linetype = "solid", color = "black", size = 1) +
  labs(x = "log CFU", y = "Count") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )
?labs
#right dashed line inoculum
#left straight line 50% of inoculum

?gghistogram
# correlation -------------------------------------------------------------

library(Hmisc)
library(corrplot)


# Calculate the correlation matrix with p-values
# Combine the matrices
combined_matrix <- as.matrix(join_meta_MGIA)

# Calculate the correlation matrix with p-values
cor_matrix <- rcorr(combined_matrix, type="spearman")

# Extract the correlations between cytokines and metabolites
n_cytokines <- length(1)
n_metabolites <- length(metabolite_select)
cor_subset <- cor_matrix$r[1:n_cytokines, (n_cytokines+1):(n_cytokines+n_metabolites)]
p_subset <- cor_matrix$P[1:n_cytokines, (n_cytokines+1):(n_cytokines+n_metabolites)]

cor_subset_t <- as.data.frame(cor_subset)
p_subset_t <- as.data.frame(p_subset)

correlation <- add_column(cor_subset_t,p_subset_t)
correlation <- correlation %>%
  rownames_to_column("ionIdx")
correlation <- left_join(correlation,metabolite_list_300BCG, by="ionIdx")
correlation$p.adj.BH <- p.adjust(correlation$p_subset, "BH")
correlation$p.adj.Bonferroni <- p.adjust(correlation$p_subset, "bonferroni")

correlation <- correlation %>%
  rename(label = `label (bona fide)`)

correlation <- correlation %>%
  arrange(-desc(cor_subset))

correlation <- correlation %>% filter(
  score == 100 
)


correlation.sig <- correlation %>% filter(
  score == 100 
) %>%
  filter(p_subset<0.1)



write.csv(correlation,"300BCG_correlationMGIA_sig.csv")

neg_corr <- correlation.sig %>% filter(cor_subset < 0)
pos_corr <- correlation.sig %>% filter(cor_subset >0)

top10percent_neg <- neg_corr %>%
  arrange(-desc(cor_subset)) %>%
  top_frac(0.1, -cor_subset)

top10percent_pos <- pos_corr %>%
  arrange(-desc(cor_subset)) %>%
  top_frac(0.1, cor_subset)

top_cor_meta <- add_case(top10percent_pos,top10percent_neg)

top15_neg <- neg_corr %>%
  arrange(-desc(cor_subset)) %>%
  slice_min(n=15,cor_subset)

top15_pos <- pos_corr %>%
  arrange(-desc(cor_subset)) %>%
  slice_max(n=15,cor_subset)

top15_cor_meta <- add_case(top15_pos,top15_neg)
write.csv(top15_cor_meta,"300BCG_top15_cor_MGIA.csv")
write.csv(top_cor_meta,"300BCG_top_cor_MGIA.csv")
write.csv(correlation,"300BCG_all_cor_MGIA.csv")
write.csv(correlation.sig,"300BCG_all_cor.sig_MGIA.csv")

ggpubr::ggbarplot(correlation.sig,x="label",y="cor_subset",
                  sort.val = "asc",
                  orientation  = "reverse"
                  )

join_meta_MGIA <- join_meta_MGIA %>%
  rename(MGIA = `pre vacc`)



# correlation plot --------------------------------------------------------



correlation.filter <- correlation %>%
  filter(p.adj.BH < 0.1)
# Cluster the columns (cytokines)
hc <- hclust(dist(t(cor_subset_t)), method = "complete")
cytokine_order <- hc$order

# Create a correlation plot with p-values
corrplot(cor_subset_t, method="color", type="full", 
         p.mat = p_subset_t, sig.level = 0.05, insig = "blank",
         addCoef.col = "black", # Add correlation coefficients
         tl.col="black", tl.srt=45, # Text label color and rotation
         cl.align.text = "l", # Align correlation coefficients to the left
         cl.offset = 0.3, # Offset correlation coefficients
         number.cex = 0.75, # Size of correlation coefficients
         diag = FALSE) # Hide correlation coefficients on the diagonal


# Create the scatter plot
ggplot(join_meta_MGIA, aes(x = `1541`, y = `MGIA`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_minimal() +
  labs(#x = "9", 
       y = "Log CFU",
       title = "Beta-glucan Level vs IL-8 E.coli")

z <- rcorr(join_meta_MGIA$`1528`, join_meta_MGIA$MGIA,"spearman")

ggp(join_meta_MGIA, "9", "MGIA")

str(join_meta_MGIA)


# volcano plot ------------------------------------------------------------
str(correlation)

keyvals.colour <- ifelse(
  correlation$p_subset > 0.05, 'black',
  ifelse(correlation$p_subset< 0.05 & correlation$cor_subset < 0, 'blue',
         ifelse(correlation$p_subset < 0.05 & correlation$cor_subset > 0,'red',
                'black')))
keyvals.colour[is.na(keyvals.colour)] <- 'lightgrey'
names(keyvals.colour)[keyvals.colour == 'lightgrey'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Positively correlated with CFU'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Negatively correlated with CFU'

library(EnhancedVolcano)
# First, filter the data for the specific metabolites
specific_metabolites <- c("Proline", "Prostaglandin F2a", "Indole-3-acetate", "LysoPE(0:0/20:2(11Z,14Z))", "Ternatin")
filtered_correlation <- correlation[correlation$label %in% specific_metabolites, ]

# Separate positive and negative correlations
positive_corr <- filtered_correlation[filtered_correlation$cor_subset > 0, ]
negative_corr <- filtered_correlation[filtered_correlation$cor_subset < 0, ]

# Get the top 2 with lowest p-value for each
top_positive <- top_cor_meta[order(positive_corr$p_subset), ][1:2, ]
top_negative <- negative_corr[order(negative_corr$p_subset), ][1:2, ]

# Combine the top metabolites
top_metabolites <- rbind(top_positive, top_negative)

# Create a vector of labels, where only the top metabolites are labeled
label_vector <- ifelse(correlation$label %in% specific_metabolites, correlation$label, "")

# Now, modify the EnhancedVolcano function call
plot <- EnhancedVolcano(correlation,
                        lab = correlation$label,
                        titleLabSize = 1,
                        subtitleLabSize = 1,
                        subtitle = "",
                        x = "cor_subset",
                        y = 'p_subset',
                        legendPosition = "none",
                        legendLabSize = 15,
                        legendIconSize = 5,
                        legendDropLevels = FALSE,
                        #ylab = paste("-Log10",bquote(~ italic(P))),
                        xlim = c(-2,2),
                        ylim = c(0, 5),
                        colCustom = keyvals.colour,
                        title = '',
                        pCutoff = 0.05,
                        FCcutoff = 0,
                        pointSize = 5.0,
                        labSize = 5.25,
                        colAlpha = 0.15,
                        gridlines.minor = FALSE,
                        xlab = "Spearman's rho",
                        drawConnectors = TRUE,
                        widthConnectors = 1)

#png(file="C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\Paper\\figure5_subset_pval.v5.png",
#   width=1000, height=800)

plot

#plot2 <- plot +
#  scale_y_continuous(limits = c(0,5), breaks = c(0,1,1.30103,2,3,4,5),labels = c(0,0.1,0.05,0.01,0.001,0.0001,0.00001))
#
#plot2




# linear model ------------------------------------------------------------


MGIA_300BCG$PatientID = as.character(MGIA_300BCG$PatientID)


klinis_meta = klinis_meta |> 
  mutate(gender = factor(gender,
                         levels=c("F","M"),
                         labels=c(1,2)),
         age = floor(age),
         classification = if_else(classification=="already controlling",
                                    "control",
                                    "no control"),
         classification = factor(classification,
                                 levels=c("no control",
                                          "control")))
str(klinis_meta)                        
head(klinis_meta)

F_klinis_meta = klinis_meta |> filter(gender == "1")
M_klinis_meta = klinis_meta |> filter(gender == "2")

lm = lm(`pre vacc` ~ `40`, data = F_klinis_meta)
summary(lm)


# check for sex and metabolites -------------------------------------------

Hemoglobin_300BCG = all_300BCG |> 
  select(PatientID, `HGB(mmol/L).WB.1`) 

Hemoglobin_300BCG$`HGB(mmol/L).WB.1` = as.numeric(Hemoglobin_300BCG$`HGB(mmol/L).WB.1`)

klinis_meta = left_join(
  klinis_meta,
  Hemoglobin_300BCG,
  "PatientID"
)

# Prepare data - merge metabolites with clinical data for MGIA participants
mgia_metabolites <- klinis_meta %>%
  select(PatientID, age, gender, BMI, `HGB(mmol/L).WB.1`,
         all_of(metabolite_select))  # All metabolite columns

# Function to test association between one variable and all metabolites
test_confounder_metabolite_association <- function(confounder_var, confounder_name, data) {
  
  results_list <- list()
  
  for (metabolite in metabolite_select) {
    
    # For continuous confounders (age, BMI, hemoglobin)
    if (is.numeric(data[[confounder_var]])) {
      cor_test <- cor.test(data[[confounder_var]], data[[metabolite]],
                           method = "spearman")
      results_list[[metabolite]] <- data.frame(
        metabolite = metabolite,
        confounder = confounder_name,
        test = "Spearman",
        statistic = cor_test$estimate,
        p.value = cor_test$p.value
      )
    }
    # For categorical confounders (sex/gender)
    else if (is.factor(data[[confounder_var]])) {
      wilcox_test <- wilcox.test(data[[metabolite]] ~ data[[confounder_var]])
      results_list[[metabolite]] <- data.frame(
        metabolite = metabolite,
        confounder = confounder_name,
        test = "Wilcoxon",
        statistic = wilcox_test$statistic,
        p.value = wilcox_test$p.value
      )
    }
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Adjust for multiple testing
  results_df$p.adjusted <- p.adjust(results_df$p.value, method = "fdr")
  
  return(results_df)
}

# Test each confounder
sex_metabolites <- test_confounder_metabolite_association("gender", "Sex", mgia_metabolites)
age_metabolites <- test_confounder_metabolite_association("age", "Age", mgia_metabolites)
bmi_metabolites <- test_confounder_metabolite_association("BMI", "BMI", mgia_metabolites)
hgb_metabolites <- test_confounder_metabolite_association("HGB(mmol/L).WB.1", "HGB(mmol/L).WB.1", mgia_metabolites)


# Summarize findings
cat("=== CONFOUNDER-METABOLITE ASSOCIATIONS IN MGIA SUBSET (n=42) ===\n\n")

cat("SEX and Metabolites:\n")
cat("  Significant (p < 0.05):", sum(sex_metabolites$p.value < 0.05),
    "out of", nrow(sex_metabolites), "metabolites\n")
cat("  FDR-significant (q < 0.05):", sum(sex_metabolites$p.adjusted < 0.05), "\n\n")

cat("AGE and Metabolites:\n")
cat("  Significant (p < 0.05):", sum(age_metabolites$p.value < 0.05),
    "out of", nrow(age_metabolites), "metabolites\n")
cat("  FDR-significant (q < 0.05):", sum(age_metabolites$p.adjusted < 0.05), "\n\n")

cat("BMI and Metabolites:\n")
cat("  Significant (p < 0.05):", sum(bmi_metabolites$p.value < 0.05),
    "out of", nrow(bmi_metabolites), "metabolites\n")
cat("  FDR-significant (q < 0.05):", sum(bmi_metabolites$p.adjusted < 0.05), "\n\n")

cat("HEMOGLOBIN and Metabolites:\n")
cat("  Significant (p < 0.05):", sum(hgb_metabolites$p.value < 0.05),
    "out of", nrow(hgb_metabolites), "metabolites\n")
cat("  FDR-significant (q < 0.05):", sum(hgb_metabolites$p.adjusted < 0.05), "\n\n")


# heatmap of correlation --------------------------------------------------

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

################################################################################
# PREPARE DATA FOR HEATMAP
################################################################################
metabolite_names <- metabolite_list_300BCG %>%
  filter(ionIdx %in% metabolite_ids, score == 100) %>%
  select(ionIdx, `label (bona fide)`, score) |> 
  rename(
    name = `label (bona fide)`
  )

# Prepare data for ggplot
plot_data <- all_confounders %>%
  mutate(
    sig_label = case_when(
      p.adjusted < 0.05 ~ "FDR < 0.05",
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05 ~ "p < 0.05",
      TRUE ~ "NS"
    ),
    sig_label = factor(sig_label,
                       levels = c("FDR < 0.05", "p < 0.001", "p < 0.01",
                                  "p < 0.05", "NS"))
  )

# Add metabolite names
if(exists("metabolite_names")) {
  plot_data <- plot_data %>%
    left_join(
      metabolite_names %>%
        filter(score == 100) %>%
        select(ionIdx, name),
      by = c("metabolite" = "ionIdx")
    )
}

# Filter to only significant associations for cleaner visualization
plot_data_sig <- plot_data %>%
  filter(p.value < 0.05)

# Create ggplot heatmap
p_heatmap <- ggplot(plot_data_sig,
                    aes(x = confounder, y = factor(metabolite))) +
  geom_tile(aes(fill = statistic), color = "white", size = 0.5) +
  geom_text(aes(label = ifelse(sig_label == "FDR < 0.05", "✓", "")),
            size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "Association\nStrength",
    guide = guide_colorbar(barwidth = 1, barheight = 10)
  ) +
  facet_grid(sig_label ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Confounder-Metabolite Associations in MGIA Subset (n=42)",
    subtitle = paste0("Showing ", nrow(plot_data_sig),
                      " significant associations (p < 0.05) out of ",
                      nrow(plot_data), " total tests"),
    x = NULL,
    y = "Metabolite"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey60"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right"
  )

ggsave("MGIA_confounder_metabolite_heatmap_ggplot.pdf",
       p_heatmap, width = 10, height = 12)

cat("ggplot version saved as: MGIA_confounder_metabolite_heatmap_ggplot.pdf\n")

# linear model ------------------------------------------------------------




create_model <- function(predictor) {
  formula <- as.formula(paste0("`pre vacc` ~ `", predictor, "`"
                               )
                               )
  lm(formula, data = klinis_meta)
}


models <- map(metabolite_select, create_model)


results <- map_dfr(models, tidy, .id = "ionIdx") %>%
  filter(term != "(Intercept)" #& term != "age" 
         & term != "gender2") %>%
  mutate(ionIdx = metabolite_select[as.integer(ionIdx)])


results <- results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

results.name = left_join(results,metabolite_list_300BCG,"ionIdx")
results.name = results.name |> 
  filter(score == 100)


results.name.sig = results.name |> 
  filter(score == 100) |> 
  filter(p.value < 0.1)

write.csv(results.name.sig,"300BCG_lm_gender.metaMGIA.v2.csv")
write_rds(results.name.sig,"300BCG_lm_gender.metaMGIA.v2.rds")


z <- read.csv("3. significant_metabolites_subset.csv")

both.sig <- intersect(z$metabolites,results.name.sig$HMDB_ID1)

all.igra <- read.csv("3. result logres metabolites_subset.csv")
all.igra.met <- all.igra$metabolites

both.available <- intersect(all.igra.met,results.name$HMDB_ID1)

sig.igra.fromall <- intersect(z$metabolites,both.available)

both.sig2 <- intersect(sig.igra.fromall,results.name.sig$HMDB_ID1)

negate(intersect(sig.igra.fromall,both.available))

outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

intersect(results.name.sig$HMDB_ID1,both.available)
outersect(sig.igra.fromall,both.available)
library(Barnard)

# Create the 2x2 contingency table
data_matrix <- matrix(c(7, 10,    # First row
                        30, 78),   # Second row
                      nrow = 2, 
                      byrow = TRUE)

# Add row and column names for clarity
rownames(data_matrix) <- c("Associated with IGRA", "Not associated with IGRA")
colnames(data_matrix) <- c("Associated with MGIA", "Not associated with MGIA")

# View the table
print(data_matrix)

# Perform Barnard's exact test
#barnard.test(data_matrix)

# Or more explicitly:
barnard.test(7,    # successes in group 1 (Associated with IGRA)
             17,    # total in group 1 (7 + 10)
             30,    # successes in group 2 (Not associated with IGRA)
             108)   # total in group 2 (30 + 78)

# volcano linear model ----------------------------------------------------
str(results.name)

keyvals.colour <- ifelse(
  results.name$p.value > 0.1, 'black',
  ifelse(results.name$p.value< 0.1 & results.name$estimate < 0, 'blue',
         ifelse(results.name$p.value < 0.1 & results.name$estimate > 0,'red',
                'black')))
keyvals.colour[is.na(keyvals.colour)] <- 'lightgrey'
names(keyvals.colour)[keyvals.colour == 'lightgrey'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Positively correlated with CFU'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Negatively correlated with CFU'

library(EnhancedVolcano)


results.name <- results.name %>%
  mutate(`label (bona fide)` = case_when(
    `label (bona fide)` == "Indole-3-acetate" ~ "Indoleacetic acid",
    `label (bona fide)` == "Indoleacitic acid" ~ "Indoleacetic acid",
    `label (bona fide)` == "13S-hydroxyoctadecadienoic acid" ~ "13-HODE",
    `label (bona fide)` == "4-Guanidinobutanoate" ~ "4-Guanidinobutanoic acid",
    `label (bona fide)` == "Prostaglandin F2a" ~ "Prostaglandin F2 alpha",
    TRUE ~ `label (bona fide)`
  ))



specific_metabolites = results.name %>% filter(HMDB_ID1 %in% both.sig2)
specific_metabolites = specific_metabolites$`label (bona fide)`



# First, filter the data for the specific metabolites
#specific_metabolites <- c("Proline", "Prostaglandin F2a", "Indole-3-acetate", "LysoPE(0:0/20:2(11Z,14Z))", "Ternatin")
filtered_correlation <- results.name[results.name$`label (bona fide)` %in% specific_metabolites, ]

# Separate positive and negative correlations
positive_corr <- filtered_correlation[filtered_correlation$estimate > 0, ]
negative_corr <- filtered_correlation[filtered_correlation$estimate < 0, ]

# Get the top 2 with lowest p-value for each
top_positive <- positive_corr[order(positive_corr$p.value), ][1:2, ]
top_negative <- negative_corr[order(negative_corr$p.value), ][1:2, ]

# Combine the top metabolites
top_metabolites <- rbind(top_positive, top_negative)

# Create a vector of labels, where only the top metabolites are labeled
label_vector <- ifelse(results.name$`label (bona fide)` %in% specific_metabolites, results.name$`label (bona fide)`, "")

plot <- EnhancedVolcano(results.name,
                        lab = label_vector,
                        titleLabSize = 1,
                        subtitleLabSize = 1,
                        subtitle = "",
                        x = "estimate",
                        y = 'p.value',
                        legendPosition = "none",
                        legendLabSize = 15,
                        legendIconSize = 5,
                        legendDropLevels = FALSE,
                        ylim = c(0, 4.5),
                        colCustom = keyvals.colour,
                        title = '',
                        pCutoff = 0.1,
                        FCcutoff = 0,
                        pointSize = 5.0,
                        labSize = 5.25,
                        colAlpha = 0.15,
                        gridlines.minor = FALSE,
                        xlab = "Linear Regression Coefficient Estimate",
                        ylab = bquote(~italic("P")-value),
                        drawConnectors = TRUE,
                        widthConnectors = 1,
                        max.overlaps = Inf,
                        boxedLabels = TRUE,
                        min.segment.length = 0.5
) +
  scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c("0.1", "0.01", "0.001", "0.0001"))

plot

ggsave("fig5a.png",
       width = 20,
       height = 15,
       units = "cm")

# stratified by sex -------------------------------------------------------

sex_strat = merged_data |> 
  select(
    PatientID,
    `pre vacc`,
    gender,
    significant.both$ionIdx
  ) |> 
  mutate(
    gender = factor(
      gender,
      levels = c("F","M")
    )
  )



# linear model F and M ----------------------------------------------------


create_model <- function(predictor) {
  formula <- as.formula(paste0("`pre vacc` ~ `", predictor, "`"# + gender"
  )
  )
  lm(formula, data = F_klinis_meta)
}


models <- map(metabolite_select, create_model)


results <- map_dfr(models, tidy, .id = "ionIdx") %>%
  filter(term != "(Intercept)" #& term != "age" 
         & term != "gender2") %>%
  mutate(ionIdx = metabolite_select[as.integer(ionIdx)])


results <- results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

results.name_F = left_join(results,metabolite_list_300BCG,"ionIdx")
results.name_F = results.name_F |> 
  filter(score == 100)


results.name_F.sig = results.name_F |> 
  filter(score == 100) |> 
  filter(p.value < 0.1)

write.csv(results.name_F.sig,"300BCG_lm_female.metaMGIA.v2.csv")
write_rds(results.name_F.sig,"300BCG_lm_female.metaMGIA.v2.rds")




create_model <- function(predictor) {
  formula <- as.formula(paste0("`pre vacc` ~ `", predictor, "`"# + gender"
  )
  )
  lm(formula, data = M_klinis_meta)
}


models <- map(metabolite_select, create_model)


results <- map_dfr(models, tidy, .id = "ionIdx") %>%
  filter(term != "(Intercept)" #& term != "age" 
         & term != "gender2") %>%
  mutate(ionIdx = metabolite_select[as.integer(ionIdx)])


results <- results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

results.name_M = left_join(results,metabolite_list_300BCG,"ionIdx")
results.name_M = results.name_M |> 
  filter(score == 100)


results.name_M.sig = results.name_M |> 
  filter(score == 100) |> 
  filter(p.value < 0.1)



z <- read.csv("3. significant_metabolites_subset.csv")

both.sig <- intersect(z$metabolites,results.name_F.sig$HMDB_ID1)

all.igra <- read.csv("3. result logres metabolites_subset.csv")
all.igra.met <- all.igra$metabolites

both.available <- intersect(all.igra.met,results.name$HMDB_ID1)

sig.igra.fromall <- intersect(z$metabolites,both.available)

both.sig2 <- intersect(sig.igra.fromall,results.name_F.sig$HMDB_ID1)

negate(intersect(sig.igra.fromall,both.available))

outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

intersect(results.name.sig$HMDB_ID1,both.available)

# check normality ---------------------------------------------------------

# 4. Combined visual and statistical approach
library(car)
qqPlot(join_meta_MGIA$MGIA)
qqPlot(join_meta_MGIA$`858`)
# 5. Check both variables at once
# Create a function for comprehensive normality check
check_normality <- function(data, variable) {
  # Shapiro test
  shapiro_test <- shapiro.test(data[[variable]])
  
  # Create plots
  par(mfrow = c(1,2))
  hist(data[[variable]], main=paste("Histogram of", variable))
  qqnorm(data[[variable]], main=paste("Q-Q Plot of", variable))
  qqline(data[[variable]])
  
  # Print results
  cat("\nNormality test results for", variable, ":\n")
  cat("Shapiro-Wilk statistic:", shapiro_test$statistic, "\n")
  cat("p-value:", shapiro_test$p.value, "\n")
  cat("Interpretation:", ifelse(shapiro_test$p.value > 0.05, 
                                "Data appears normal (p > 0.05)", 
                                "Data does not appear normal (p < 0.05)"), "\n\n")
}

# Use the function
check_normality(join_meta_MGIA, "MGIA")
check_normality(join_meta_MGIA, "858")  # or whatever your column name is

check_normality(join_meta_MGIA, "MGIA")
check_normality(join_meta_MGIA, "40")  # or whatever your column name is


# correlation spearman plot -----------------------------------------------

significant.both <- results.name.sig |> 
  filter(HMDB_ID1 %in% both.sig)

list.sig.both = significant.both$ionIdx

significant.both.corplot <- join_meta_MGIA |> 
  select(MGIA,list.sig.both)


# test --------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

# Function to create a single scatter plot
create_metabolite_plot <- function(data, metabolite_col, metabolite_name) {
  # Prepare data
  plot_data <- data %>%
    select(`MGIA`, all_of(metabolite_col))
  
  # Create formula using the actual column name
  formula_str <- paste("`MGIA` ~", paste("`", metabolite_col, "`", sep=""))
  
  # Calculate statistics
  model <- lm(as.formula(formula_str), data = plot_data)
  model_stats <- summary(model)
  
  # Extract statistics
  stats <- list(
    beta = round(coef(model)[2], 3),
    p_val = round(coef(model_stats)[2,4], 3),
    x_pos = max(plot_data[[metabolite_col]]),
    y_pos = max(plot_data$`MGIA`)
  )
  
  # Create plot with proper handling of column names
  ggplot(plot_data, aes(x = .data[[metabolite_col]], y = `MGIA`)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue", fill = "darkgrey", alpha = 0.5) +
    labs(
      title = metabolite_name,
      x = "Log (peak ion intensity)",
      y = "Log CFU"
    ) +
    theme_pubr() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 20
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    annotate(
      "text",
      x = (stats$x_pos-0.075),
      y = 0.5,
      label = paste("β =", stats$beta, "\n",
                    "P =", stats$p_val),
      hjust = 0,
      vjust = 1
    )
}

# Define metabolites to plot
metabolites <- list(
  list(col = "8", name = "Glycine"),
  list(col = "40", name = "Proline"),
  list(col = "304", name = "Tryptophan"),
  list(col = "193", name = "Indoleacetic acid"),
  list(col = "786", name = "Leukotriene B4"),
  list(col = "645", name = "13-HODE"),
  list(col = "858", name = "Prostaglandin F2α")
)

# Create all plots
plots <- map(metabolites, ~create_metabolite_plot(
  join_meta_MGIA, 
  .x$col, 
  .x$name
))

export = 
  join_meta_MGIA |> 
  select(MGIA,
         `304`,
         `193`,
         `786`,
         `645`,
         `858`,
         `40`,
         `8`)

write.csv(export,"Fig 5C_INFECT metabolomics.csv")

# Combine plots in one rows
combined_plot <- wrap_plots(plots, nrow = 4) +
  plot_layout(guides = "collect") & 
  theme(plot.margin = margin(5, 5, 5, 5))

# Display the combined plot
print(combined_plot)

ggsave("combined corplot.png",
       width = 20,
       height = 40,
       units = "cm")

# Optionally, save the plot with appropriate dimensions
# ggsave("metabolite_plots.png", combinhttp://127.0.0.1:13511/graphics/plot_zoom_png?width=1168&height=392ed_plot, width = 20, height = 5, dpi = 300)