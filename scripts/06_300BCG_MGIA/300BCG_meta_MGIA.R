# metabolites and MGIA 300BCG #


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

# load data ---------------------------------------------------------------

meta_cyto <- load("17.02.2021 Metabolomics 300BCG dataset.RData")
#meta_cyto <- as.data.frame(read.csv("300BCG_17.02.2021_Metabolomics_300BCG_dataset.csv"))
all_INFECT_metabolites <- read.csv("3. result logres metabolites_subset.csv") %>% select(metabolites,name)
metabolite_list_INFECT <- read.csv("3. significant_metabolites_subset.csv") %>% select(metabolites,name)
metabolite_list_300BCG <- read_excel("300BCG_08.08.2023_Metabolites_HMDBIDS.xlsx")
MGIA_300BCG <- read_excel("300BCG_MGIA.xlsx") %>%
  rename(PatientID = ID)

not_sig_INFECT <- setdiff(all_INFECT_metabolites$metabolites,metabolite_list_INFECT$metabolites)

metabolite_list_300BCG$ionIdx <- as.numeric(metabolite_list_300BCG$ionIdx)
str(metabolite_list_300BCG)


subset_metabolite_list_300BCG <- metabolite_list_300BCG %>%
  filter(score == 100) %>%
  filter(if_any(starts_with("HMDB_ID"), ~ . %in% metabolite_list_INFECT$metabolites))


subset_metabolite_list_300BCG_v2 <- metabolite_list_300BCG %>%
  filter(score == 100) %>%
  filter(HMDB_ID1 %in% metabolite_list_INFECT$metabolites)

subset_notsig_metabolite_list_300BCG_v2 <- metabolite_list_300BCG %>%
  filter(score == 100) %>%
  filter(HMDB_ID1 %in% not_sig_INFECT)

all_metabolite <- c(names(mydata[53:1425]))
#subset_metabolite_list_300BCG_v3 <- metabolite_list_300BCG %>%
#  filter(score == 100) %>%
#  filter(ionIdx %in% all_metabolite)

subset_metabolite_list_300BCG_v3 <- metabolite_list_300BCG %>%
  filter(score == 100) %>%
  filter(HMDB_ID1 %in% all_INFECT_metabolites$metabolites)



subset_metabolite_list_300BCG_v2 <- left_join(subset_metabolite_list_300BCG_v2,metabolite_list_INFECT, by = c("HMDB_ID1" ="metabolites"))
subset_notsig_metabolite_list_300BCG_v2 <- left_join(subset_notsig_metabolite_list_300BCG_v2,all_INFECT_metabolites, by = c("HMDB_ID1" ="metabolites"))
subset_metabolite_list_300BCG_v3 <- left_join(subset_metabolite_list_300BCG_v3,all_INFECT_metabolites, by = c("HMDB_ID1" ="metabolites"))



#metabolite_select <- c(as.character(subset_metabolite_list_300BCG$ionIdx))

metabolite_select <- c(as.character(subset_metabolite_list_300BCG_v2$ionIdx))
metabolite_select_notsig <- c(as.character(subset_notsig_metabolite_list_300BCG_v2$ionIdx))

metabolite_select_all <- c(as.character(subset_metabolite_list_300BCG_v3$ionIdx))

meta_300BCG.filtered <- mydata %>% select(PatientID, Gender, Age, BMI, metabolite_select) %>%
  filter(PatientID %in% MGIA_300BCG$PatientID)

meta_300BCG.filtered.notsig <- mydata %>% select(PatientID, Gender, Age, BMI, metabolite_select_notsig) %>%
  filter(PatientID %in% MGIA_300BCG$PatientID)

meta_300BCG.all <- mydata %>% select(PatientID, Gender, Age, BMI, metabolite_select_all) %>%
  filter(PatientID %in% MGIA_300BCG$PatientID)

join_meta_MGIA <- left_join(meta_300BCG.filtered,MGIA_300BCG,by="PatientID")
name <- c(subset_metabolite_list_300BCG_v2$name)
#name <- subset_metabolite_list_300BCG_v3 %>% filter(ionIdx %in% x)
#name <- c(subset_metabolite_list_300BCG_v3$`label (bona fide)`)
join_meta_MGIA <- join_meta_MGIA %>%
  rename_with(~ name, matches("^\\d+$"))

join_meta_MGIA <- join_meta_MGIA %>%
  mutate(classification2 = if_else(classification=="already controlling","control","no control"))


join_meta_MGIA$classification <- factor(join_meta_MGIA$classification,
                                        levels = c("already controlling","acquired control","no control"),
                                        labels = c("already controlling","acquired control","no control"))




join_meta_MGIA.notsig <- left_join(meta_300BCG.filtered.notsig,MGIA_300BCG,by="PatientID")
name.notsig <- c(subset_notsig_metabolite_list_300BCG_v2$name)
#name <- subset_metabolite_list_300BCG_v3 %>% filter(ionIdx %in% x)
#name <- c(subset_metabolite_list_300BCG_v3$`label (bona fide)`)
join_meta_MGIA.notsig <- join_meta_MGIA.notsig %>%
  rename_with(~ name.notsig, matches("^\\d+$"))

join_meta_MGIA.notsig <- join_meta_MGIA.notsig %>%
  mutate(classification2 = if_else(classification=="already controlling","control","no control"))


join_meta_MGIA.notsig$classification <- factor(join_meta_MGIA.notsig$classification,
                                        levels = c("already controlling","acquired control","no control"),
                                        labels = c("already controlling","acquired control","no control"))



join_meta_MGIA.all <- left_join(meta_300BCG.all,MGIA_300BCG,by="PatientID")
join_meta_MGIA.all <- join_meta_MGIA.all %>%
  mutate(classification2 = if_else(classification=="already controlling","control","no control"))
join_meta_MGIA.notsig$classification2 <- factor(join_meta_MGIA.notsig$classification2,
                                               levels = c("no control","control"),
                                               labels = c("no control","control")
)

name.all <- c(subset_metabolite_list_300BCG_v3$name)
#name <- subset_metabolite_list_300BCG_v3 %>% filter(ionIdx %in% x)
#name <- c(subset_metabolite_list_300BCG_v3$`label (bona fide)`)
join_meta_MGIA.all <- join_meta_MGIA.all %>%
  rename_with(~ name.all, matches("^\\d+$"))

# reshape long ------------------------------------------------------------

## Reshape the data to long format
#join_meta_MGIA_long <- join_meta_MGIA %>%
#  pivot_longer(cols = c(name), 
#               names_to = "metabolite", 
#               values_to = "value")
#
## Reshape the data to long format
#join_meta_MGIA.notsig_long <- join_meta_MGIA.notsig %>%
#  pivot_longer(cols = c(name.notsig), 
#               names_to = "metabolite", 
#               values_to = "value")

join_meta_MGIA.all_long <- join_meta_MGIA.all %>%
  pivot_longer(cols = c(name.all), 
               names_to = "metabolite", 
               values_to = "value")


# outliers ----------------------------------------------------------------


outliers <- join_meta_MGIA.all_long %>%
  group_by(metabolite) %>%
  identify_outliers(value)

#outliers.notsig <- join_meta_MGIA.notsig_long %>%
#  group_by(metabolite) %>%
#  identify_outliers(value)

# Summarize the number of outliers per metabolite
outlier_summary <- outliers %>%
  group_by(metabolite) %>%
  summarise(
    total_outliers = n(),
    extreme_outliers = sum(is.extreme == TRUE)
  ) %>%
  arrange(desc(total_outliers))

# Summarize the number of outliers per metabolite
#outlier_summary.notsig <- outliers.notsig %>%
#  group_by(metabolite) %>%
#  summarise(
#    total_outliers = n(),
#    extreme_outliers = sum(is.extreme == TRUE)
#  ) %>%
#  arrange(desc(total_outliers))


ggplot(outlier_summary, aes(x = reorder(metabolite, -total_outliers), y = total_outliers)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_bar(aes(y = extreme_outliers), stat = "identity", fill = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Number of Outliers per Metabolite",
       x = "Metabolite",
       y = "Number of Outliers") +
  scale_y_continuous(expand = c(0, 0)) +
  annotate("text", x = Inf, y = Inf, label = "Dark blue: Extreme outliers", 
           hjust = 1.1, vjust = 2, size = 3)

ggsave("outliers_per_metabolite.png", width = 15, height = 8)

multi_outlier_samples <- outliers %>%
  #outliers.notsig %>%
  group_by(PatientID) %>%
  summarise(
    outlier_count = n(),
    metabolites = paste(metabolite, collapse = ", ")
  ) %>%
  arrange(desc(outlier_count))


# Create a dataframe of outliers to remove
outliers_to_remove <- outliers %>%
  select(PatientID, metabolite)

join_meta_MGIA.all_long_cleaned <- join_meta_MGIA.all_long %>%
  anti_join(outliers_to_remove, by = c("PatientID", "metabolite"))

# Reshape the data back to wide format
join_meta_MGIA.all_cleaned <- join_meta_MGIA.all_long_cleaned %>%
  pivot_wider(names_from = metabolite, values_from = value)  # Remove the temporary row_id column

# Print summary of removed outliers
outlier_summary <- outliers %>%
  group_by(metabolite) %>%
  summarise(outliers_removed = n())

# Calculate the percentage of data points removed
total_datapoints <- nrow(join_meta_MGIA.all) * (ncol(join_meta_MGIA.all) - 6)  # Subtract non-metabolite columns
removed_datapoints <- sum(outlier_summary$outliers_removed)
percent_removed <- (removed_datapoints / total_datapoints) * 100

cat(sprintf("Total data points removed: %d (%.2f%%)\n", 
            removed_datapoints, percent_removed))

# Compare the dimensions of the original and cleaned datasets
cat("Original dataset dimensions:", dim(join_meta_MGIA.all), "\n")
cat("Cleaned dataset dimensions:", dim(join_meta_MGIA.all_cleaned), "\n")


# cleaned df --------------------------------------------------------------

join_meta_MGIA <- join_meta_MGIA.all_cleaned %>% select(colnames(join_meta_MGIA.all_cleaned[1:11]), name)

join_meta_MGIA.notsig <- join_meta_MGIA.all_cleaned %>% select(colnames(join_meta_MGIA.all_cleaned[1:11]), name.notsig)

join_meta_MGIA_long <- join_meta_MGIA.all_long_cleaned %>% filter(metabolite %in% name)
join_meta_MGIA.notsig_long <- join_meta_MGIA.all_long_cleaned %>% filter(metabolite %in% name.notsig)

# normality ---------------------------------------------------------------


# Density plot
density <- ggdensity(join_meta_MGIA_long, x="value",facet.by = "metabolite", fill = "lightgray", scales = "free")
# QQ plot
qqplot <- ggqqplot(join_meta_MGIA_long, x="value",facet.by = "metabolite", scales = "free")

density_MGIA <- ggdensity(join_meta_MGIA_long, x="pre vacc",facet.by = "metabolite", fill = "lightgray", scales = "free")

join_meta_MGIA_long$mgia <- join_meta_MGIA_long$`pre vacc`

normality_results_MGIA <- join_meta_MGIA_long %>%
  filter(metabolite=="Histamine") %>%
  shapiro_test(mgia) %>%
  mutate(
    normal = ifelse(p > 0.05, "Yes", "No")
  )

# Density plot
density.notsig <- ggdensity(join_meta_MGIA.notsig_long, x="value",facet.by = "metabolite", fill = "lightgray", scales = "free")
# QQ plot
qqplot.notsig <- ggqqplot(join_meta_MGIA.notsig_long, x="value",facet.by = "metabolite", scales = "free")

normality_results <- join_meta_MGIA_long %>%
  group_by(metabolite, classification2) %>%
  shapiro_test(value) %>%
  ungroup() %>%
  mutate(
    normal = ifelse(p > 0.05, "Yes", "No")
  )

normality_results.notsig <- join_meta_MGIA.notsig_long %>%
  group_by(metabolite, classification2) %>%
  shapiro_test(value) %>%
  ungroup() %>%
  mutate(
    normal = ifelse(p > 0.05, "Yes", "No")
  )



# -------------------------------------------------------------------------

gtsummary::tbl_summary(join_meta_MGIA %>% select(Gender, Age, BMI
                                                 #,metabolite_select
                                                 ,name
                                                 ,classification2), by = "classification2") %>%
  add_p()



gtsummary::tbl_summary(join_meta_MGIA.notsig %>% select(Gender, Age, BMI
                                                 #,metabolite_select
                                                 ,name.notsig
                                                 ,classification2), by = "classification2") %>%
  add_p()




#join_meta_MGIA <- column_to_rownames(join_meta_MGIA,"PatientID")



# comparison between group ------------------------------------------------


# Perform t-test for each metabolite
results_sig <- join_meta_MGIA_long %>%
  group_by(metabolite) %>%
  t_test(value ~ classification2) %>%
  #wilcox_effsize(value ~ classification2) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

results_sig.effectsize <- join_meta_MGIA_long %>%
  group_by(metabolite) %>%
  wilcox_effsize(value ~ classification2)

results_sig <- left_join(results_sig, results_sig.effectsize, by = "metabolite")


results_sig.metabolites <- as.character(unlist(c(results_sig %>% filter(p < 0.05) %>% select(metabolite))))

# Perform t-test for each metabolite
results_not.sig <- join_meta_MGIA.notsig_long %>%
  group_by(metabolite) %>%
  wilcox_test(value ~ classification2) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

results_not.sig.effectsize <- join_meta_MGIA.notsig_long %>%
  group_by(metabolite) %>%
  wilcox_effsize(value ~ classification2)

results_not.sig <- left_join(results_not.sig, results_not.sig.effectsize, by = "metabolite")

results_not.sig.metabolites <- unlist(c(results_not.sig %>% filter(p < 0.05) %>% select(metabolite)))


# boxplot -----------------------------------------------------------------


#ggpubr::ggboxplot(join_meta_MGIA %>% 
#                    pivot_longer(name),"classification","value",facet.by = "name", scales="free_y",
#                  add = "dotplot",
#                  short.panel.labs = TRUE) + stat_compare_means(paired=FALSE, method = "wilcox.test",ref.group = "no control",
#                                                        comparisons = list(c("already controlling","acquired control"),
#                                                                           c("already controlling","no control")),
#                                                        aes(label = paste0("", after_stat(p.format)))) + theme_pubr()
#


ggpubr::ggboxplot(join_meta_MGIA_long %>% filter(metabolite %in% results_sig.metabolites), color = "classification2",palette = c("red","blue"),"classification2","value",facet.by = "metabolite", scales="free",
                  add = "dotplot") + stat_compare_means(method = "wilcox.test") + theme_pubr()


ggpubr::ggboxplot(join_meta_MGIA.notsig_long %>% filter(metabolite %in% results_not.sig.metabolites), color = "classification2",palette = c("red","blue"),"classification2","value",facet.by = "metabolite", scales="free",
                  add = "dotplot") + stat_compare_means(method = "wilcox.test") + theme_pubr()



# PCA ---------------------------------------------------------------------


# Perform PCA on the metabolite variables
pca_result <- PCA(join_meta_MGIA[, name],
                  #join_meta_MGIA.notsig[, name.notsig],
                  #join_meta_MGIA[, 4:1296], 
                  scale = TRUE, graph = FALSE
                  )

fviz_pca_ind(pca_result,
             axes = c(1, 2),
                habillage = as.factor(join_meta_MGIA$classification2),
                palette = c("#209ebf","red"),
                addEllipses = TRUE,
                ellipse.level = 0.95,
                legend.title = "Classification",
                title = "PCA stratified by MGIA result",
             ggtheme = theme_classic2(),
             labelsize = 4,
             xlab = paste0("PC1 (", round(pca_result$eig[, 2][1], 2), "%)"),
             ylab = paste0("PC2 (", round(pca_result$eig[, 2][2], 2), "%)"),
             repel = TRUE)

fviz_pca_biplot(pca_result,
             axes = c(1, 2),
             habillage = as.factor(join_meta_MGIA$classification2),
             palette = c("#209ebf","red"),
             #palette = "startrek",
             geom = "point", pointsize = 2,
             col.var = "black",
             alpha.ind = 1,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             legend.title = "Classification",
             title = "PCA stratified by MGIA result",
             ggtheme = theme_prism(),
             labelsize = 4,
             xlab = paste0("PC1 (", round(pca_result$eig[, 2][1], 2), "%)"),
             ylab = paste0("PC2 (", round(pca_result$eig[, 2][2], 2), "%)"),
             repel = TRUE)+
  scale_shape_manual(values=c(19,19,19))



# correlation signif-------------------------------------------------------------



PreVacc_stat <- data.frame(correlation(data = data.frame(join_meta_MGIA[name]
                                                         #join_meta_MGIA[4:1296]
                                                         ), data2 = data.frame(join_meta_MGIA$`pre vacc`),
                         method = "pearson",
                         p_adjust = "none",
                         easystats.standardize_names = TRUE ))
PreVacc_stat$name <- name
PreVacc_stat.sig <- as.character(unlist(PreVacc_stat %>% filter(p<0.05) %>% select(name)))
PreVacc_stat.sig
PreVacc_cor <- join_meta_MGIA %>%
  pivot_longer(cols = c(PreVacc_stat.sig), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "pre vacc",
            color = "black",
            alpha = 0.3,
            facet.by = "column", ncol = 3, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "pearson",
           #cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..rr.label.., ", p = ", ..p.label..)),
            label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 4
           ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "log MFI", y = "log CFU", title = "Correlation of pre-vaccination MGIA logCFU with metabolites associated with IGRA at follow-up in INFECT") +
  theme_pubr()+
  theme(panel.spacing = unit(2, "lines"))
#1150x400
PreVacc_cor



# correlation notsignif ---------------------------------------------------


PreVacc_stat.notsig <- data.frame(correlation(data = data.frame(join_meta_MGIA.notsig[name.notsig]
                                                         #join_meta_MGIA.notsig[4:1296]
), data2 = data.frame(join_meta_MGIA.notsig$`pre vacc`),
method = "pearson",
p_adjust = "none",
easystats.standardize_names = TRUE ))

PreVacc_stat.notsig$name <- name.notsig
PreVacc_stat.notsig.meta <- as.character(unlist(PreVacc_stat.notsig %>% filter(p<0.05) %>% select(name)))
PreVacc_stat.notsig.meta
PreVacc_cor.notsig <- join_meta_MGIA.notsig %>%
  pivot_longer(cols = c(PreVacc_stat.notsig.meta), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "pre vacc",
            facet.by = "column", ncol = 3, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "pearson",
           #cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..rr.label.., ", p = ", ..p.label..)),
           label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 3
  ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "logMFI", y = "logCFU", title = "Correlation Plots: pre-vaccination MGIA logCFU with metabolites") +
  theme_pubr()

PreVacc_cor.notsig



# linear regression -------------------------------------------------------
library(tidyverse)
library(broom)



# test --------------------------------------------------------------------

# Function to safely extract vector from potentially list column
safe_extract <- function(df, col_name) {
  col <- df[[col_name]]
  if (is.list(col)) {
    return(unlist(col))
  } else {
    return(col)
  }
}

# Function to safely run linear model
safe_lm <- function(metabolite, data) {
  tryCatch({
    # Safely extract vectors
    y <- safe_extract(data, "pre vacc")
    x <- safe_extract(data, metabolite)
    age <- safe_extract(data, "Age")
    bmi <- safe_extract(data, "BMI")
    gender <- safe_extract(data, "Gender")
    
    # Create a new dataframe for the model
    model_data <- data.frame(y = y, x = x, Age = age, BMI = bmi, Gender = gender)
    
    # Run the model
    model <- lm(y ~ x 
                #+ Age
                #+ BMI 
                #* Gender
                # + Gender
                , data = model_data)
    tidy_result <- tidy(model)
    
    # Extract the row for the metabolite (usually the second row)
    metabolite_row <- tidy_result[2, ]
    c(metabolite = metabolite, 
      estimate = metabolite_row$estimate, 
      std.error = metabolite_row$std.error, 
      statistic = metabolite_row$statistic, 
      p.value = metabolite_row$p.value)
  }, error = function(e) {
    c(metabolite = metabolite, 
      estimate = NA, 
      std.error = NA, 
      statistic = NA, 
      p.value = NA)
  })
}



# Initialize result.lin.mod as a list
result.lin.mod <- list()

# Run models for significant  metabolites
for (i in name) {
  result.lin.mod[[i]] <- safe_lm(i, join_meta_MGIA)
}


# #Run models for not significant metabolites
#for (i in name.notsig) {
#  result.lin.mod[[i]] <- safe_lm(i, join_meta_MGIA.notsig)
#}

# Convert list to dataframe
result.lin.mod <- do.call(rbind, result.lin.mod) %>% as.data.frame()

# Convert to a tibble for easier handling
result.lin.mod <- as_tibble(result.lin.mod)

# Adjust p-values for multiple comparisons
result.lin.mod <- result.lin.mod %>%
  mutate(p.adjusted = p.adjust(p.value, method = "BH"))

#
## Identify significant metabolites (you can adjust the threshold as needed)
#sig.metabolites <- result.lin.mod %>%
#  filter(p.adjusted < 0.05) %>%
#  arrange(p.adjusted)


######################### volcano fullset ########################


sig.metabolites <- data.frame(result.lin.mod)
#sig.metabolites$p.adj = p.adjust(sig.metabolites$Pr...t.., method = "BH")

#sig.metabolites <- sig.metabolites %>% filter(Pr...t.. < 0.05)

#sig.metabolites <- left_join(sig.metabolites,metabolite_list_300BCG %>% filter(score == 100) %>% select(ionIdx,`label (bona fide)`,HMDB_ID1) %>% rename(V1 = ionIdx) %>% mutate(V1 = as.character(V1)))

#sig.metabolites <- sig.metabolites %>%
#  dplyr::rename(metabolites = V1)

#name <- read.csv("C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\1. todia analisis\\boston_metabolites_name.csv")
#name <- name[2:3]
##colnames(name) <- c("metabolites","name")#

#sig.metabolites <- left_join(sig.metabolites,name)

#gt(sig.metabolites)

str(sig.metabolites)
sig.metabolites$p.value <- as.numeric(sig.metabolites$p.value)
sig.metabolites$estimate <- as.numeric(sig.metabolites$estimate)
sig.metabolites$statistic <- as.numeric(sig.metabolites$statistic)
sig.metabolites$std.error <- as.numeric(sig.metabolites$std.error)

keyvals.colour <- ifelse(
  sig.metabolites$`p.value` > 0.1, 'black',
  ifelse(sig.metabolites$`p.value` < 0.1 & sig.metabolites$estimate < 0, 'blue',
         ifelse(sig.metabolites$`p.value` < 0.1 & sig.metabolites$estimate > 0,'red',
                'black')))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Increase logCFU'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Decrease logCFU'


library(EnhancedVolcano)

plot <- EnhancedVolcano(sig.metabolites,
                        lab = sig.metabolites$metabolite,
                        titleLabSize = 1,
                        subtitleLabSize = 1,
                        subtitle = "",
                        #x='logFC',
                        x = 'estimate',
                        #y = 'p.adj',
                        #y='P.Value',
                        y =  'p.value',
                        legendPosition = "top",
                        legendLabSize = 15,
                        legendIconSize = 5,
                        legendDropLevels = FALSE,
                        ylab = bquote(~italic("P value")),
                        #xlim = c(-3,3),
                        #ylim = c(0, 2),
                        colCustom = keyvals.colour,
                        title = '',
                        pCutoff = 0.1,
                        FCcutoff = 0,
                        pointSize = 5.0,
                        labSize = 5.25,
                        colAlpha = 0.3,
                        gridlines.minor = FALSE,
                        xlab = "Linear Regression Coefficient Estimate",
                        drawConnectors = FALSE,
                        caption =NULL,
                        #caption = paste0("total = ", nrow(sig.metabolites), " metabolites"),
                        widthConnectors = 1)

#png(file="C:\\Users\\todia\\Google Drive\\1. PhD\\TP0001_INFECT\\Metabolomics\\Paper\\figure5_subset_pval.v5.png",
#   width=1000, height=800)

plot + theme_prism()

plot2 <- plot +
  scale_y_continuous(limits = c(0,2), breaks = c(0,1,1.30103,2),labels = c(0,0.1,0.05,0.01))

plot2 + theme_prism()



# chi-squared -------------------------------------------------------------

# Create the 2x2 contingency table
table_data <- matrix(c(3, 14, 7, 101), nrow = 2, byrow = TRUE)

# Set row and column names
rownames(table_data) <- c("Associated with EC", "Not associated with EC")
colnames(table_data) <- c("Significant correlation with MGIA", "Non-significant correlation with MGIA")

# Create the table
contingency_table <- as.table(table_data)

library("graphics")
mosaicplot(contingency_table, shade = TRUE, las=2
           #,main = "housetasks"
           )

# install.packages("vcd")
library("vcd")
# plot just a subset of the table
assoc(head(dt, 5), shade = TRUE, las=3)

# Data
ec_associated <- c(3, 14)  # Significant, Non-significant
not_ec_associated <- c(7, 101)  # Significant, Non-significant

# Calculate percentages
ec_percentage <- ec_associated[1] / sum(ec_associated) * 100
not_ec_percentage <- not_ec_associated[1] / sum(not_ec_associated) * 100

# Create a data frame
data <- data.frame(
  Group = c("Associated with\nIGRA at follow-up", "Not associated with\nIGRA at follow-up"),
  Percentage = c(ec_percentage, not_ec_percentage)
)

# Create the bar plot
barplot <- ggplot(data, aes(x = Group, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 1), width = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            vjust = -0.5, size = 5) +
  ylim(0, max(data$Percentage) * 1.1) +  # Adjust y-axis limit to fit labels
  labs(title = "Percentage of metabolites significantly correlated with MGIA",
       y = "Percentage", 
       x = "") +
  theme_prism() +
  theme(legend.position = "none",  # Remove legend as it's redundant
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Rotate x-axis labels
        panel.spacing = unit(0.1, "lines"))   # Remove legend as it's redundant

library(viridis)
barplot + scale_fill_brewer(palette = "Set2")


#barplot + scale_fill_viridis(discrete = TRUE)
# Print the plot
print(barplot)
#
#chisq <- janitor::chisq.test(contingency_table, correct = FALSE)
#chisq
#
#fisher_test <- janitor::fisher.test(contingency_table)
#print(fisher_test)

library(DescTools)
DescTools::BarnardTest(contingency_table)

# sex interaction ---------------------------------------------------------

PreVacc_stat.sig.interaction <- as.character(unlist(sig.metabolites %>% filter(p.value<0.1) %>% select(metabolite)))

PreVacc_stat.sig.interaction

PreVacc_cor.interaction <- join_meta_MGIA %>%
  pivot_longer(cols = c(PreVacc_stat.sig.interaction), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "pre vacc",
            color = "gender",
            facet.by = "column", ncol = 3, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "gender", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "pearson",
           #cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..rr.label.., ", p = ", ..p.label..)),
           label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 3
  ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "logMFI", y = "logCFU", title = "Correlation Plots: pre-vaccination MGIA logCFU with metabolites") +
  theme_pubr()

PreVacc_cor.interaction


# cytokines ---------------------------------------------------------------

wide_dataset_cyto <- read_excel("190321 300BCG Dataset wide.xlsx")
meta <- mydata %>% select(PatientID ,metabolite_select)

v1 <- c("IL1b_T3_log.1",
        "IL1b_T4_log.1",
         "IL6_T3_log.1",
         "IL6_T4_log.1",
         "TNF_T3_log.1",
        "IFNg_W3_log.1",
        "IFNg_W4_log.1",
        "IL17_W4_log.1"
)

v3 <- c("IL1b_T3_log.3",
        "IL1b_T4_log.3",
         "IL6_T3_log.3",
         "IL6_T4_log.3",
         "TNF_T3_log.3",
        "IFNg_W3_log.3",
        "IFNg_W4_log.3",
        "IL17_W4_log.3"
)
#olink <- wide_dataset_cyto[c(1,444:516)]
olink <- wide_dataset_cyto[c(1,522:594)]

cyto <- wide_dataset_cyto %>% select(PatientID,v1)
cytokines <- left_join(meta,cyto)
cytokines <- cytokines %>%
  rename_with(~ name, matches("^\\d+$"))
MGIA_cyto <- left_join(MGIA_300BCG,cytokines)
MGIA_olink <- left_join(MGIA_300BCG,olink)
cytokines <- column_to_rownames(cytokines,"PatientID")
cytokines[18:25] <- sapply(cytokines[18:25], as.numeric)
str(cytokines)


x <- data.frame(correlation(data = cytokines[1:17], data2 = cytokines[18:25], 
                                       method = "spearman",
                                       p_adjust = "none"))


PostVacc_cor <- cytokines %>%
  pivot_longer(cols = c("13-HODE","PGF2alpha","Indoleacetic acid",
                        "Proline"), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "IL1b_T3_log.1",
            facet.by = "column", ncol = 4, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..r.., ", p = ", ..p.format..)),
           label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 3
  ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "logMFI", y = "IL-1b", title = "Correlation Plots: IL-1b upon S.aureus stimulation with metabolites") +
  theme_pubr()
PostVacc_cor


PostVacc_cor <- cytokines %>%
  pivot_longer(cols = c("13-HODE","PGF2alpha","Indoleacetic acid",
                        "Proline"), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "IL6_T4_log.1",
            facet.by = "column", ncol = 4, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..r.., ", p = ", ..p.format..)),
           label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 3
  ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "logMFI", y = "IL-6", title = "Correlation Plots: IL-6 upon Mtb stimulation with metabolites") +
  theme_pubr()
PostVacc_cor


PostVacc_cor <- cytokines %>%
  pivot_longer(cols = c("13-HODE","PGF2alpha","Indoleacetic acid",
                        "Proline"), names_to = "column", values_to = "value") %>%
  ggscatter(x = "value", y = "TNF_T3_log.1",
            facet.by = "column", ncol = 4, scales = "free_x",
            add = "reg.line",
            add.params = list(color = "blue", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(method = "spearman",
           cor.coef.name = "rho",
           #aes(label = paste0("rho = ", ..r.., ", p = ", ..p.format..)),
           label.y.npc = "bottom", 
           label.x.npc = "middle",
           size = 3
  ) +
  geom_hline(yintercept = 2.38, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = 1.19, linetype = "dotted", color = "red") +
  labs(x = "logMFI", y = "TNF", title = "Correlation Plots: IL-6 upon Mtb stimulation with metabolites") +
  theme_pubr()
PostVacc_cor





# heatmap cyto-meta -------------------------------------------------------


heatmap_corr2_mito <- x %>%
  filter(p < 0.05) %>%
  select(Parameter1, Parameter2, rho) %>%
  pivot_wider(names_from = Parameter1, values_from = rho)

metabolites_list2 <- intersect(metabolites_list, colnames(heatmap_corr2_mito))

heatmap_corr2_mito <- heatmap_corr2_mito%>%
  select(Parameter2, all_of(metabolites_list2)) %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(1) %>%
  rownames_to_column("Parameter1") %>%
  left_join(significant_metabolites_name, by = "Parameter1") %>%
  column_to_rownames("name") %>%
  select(-Parameter1)

matrix_mito <- heatmap_corr2_mito %>%
  map_dfc(as.numeric) %>%
  as.matrix()

matrix_mito[is.na(matrix_mito)] <- 0

rownames(matrix_mito) <- rownames(heatmap_corr2_mito)

rg <- max(abs(matrix_mito))


z <- data.frame(matrix_mito)
z <- rownames_to_column(z,"name")

pathway <- read.csv("INFECT_metabolic pathway edited.v3.csv")
pathway <- pathway %>% select(name,pathway_name)
pathway <- left_join(z,pathway)
pathway <- pathway %>% select(name,pathway_name)
pathway$pathway_name <- if_else(is.na(pathway$pathway_name),"",pathway$pathway_name)


mito_heatmap <- matrix_mito %>%
  ComplexHeatmap::pheatmap(scale = "none",
                           name = "Spearman's rho",
                           #labels_col = labels,
                           labels_row = pathway$pathway_name,
                           clustering_distance_cols = "correlation",
                           clustering_distance_rows = "correlation",
                           cluster_rows = T,
                           cluster_cols = T,
                           #breaks = seq(0.5, -0.5, length.out = 7),
                           breaks = seq(1, -1, length.out = 7),
                           #color = viridis(7),
                           color = colorRampPalette(c("firebrick3", "white", "navy"))(7),
                           #annotation_col = after_underscore,
                           #annotation_colors = anno_colors,
                           fontsize_row = 14,
                           fontsize = 14,
                           #display_numbers = stat_z,
                           #fontsize = 40,
                           #legend_breaks = c(-1,-0.5,0,0.5,1),
                           show_colnames = TRUE,
                           #fontsize_col = 16,
                           #cellwidth = 15,
                           cellwidth = 25,
                           border_color = "grey",
                           legend = TRUE,
                           number_color = "black",
                           fontsize_number = 9,
                           angle_col = "45"
                           #annotation_legend = FALSE,
                           #annotation = FALSE,
                           #column_split=after_underscore$Stimuli
                           #legend_labels = c(-1,-0.5,0,0.5,1),
                           #cutree_rows = 3,
                           #clustering_distance_rows = "euclid",
                           #cutree_cols = 3,
                           ,left_annotation = rowAnnotation(foo = anno_text(pathway$name,
                                                                            #z1$Source,
                                                                            location = 0.99,just = "right",
                                                                            gp = gpar(fontsize = 14)))
  )

mito_heatmap



# MGIA cyto ---------------------------------------------------------------


str(MGIA_cyto)
MGIA_cyto[25:32] <-  sapply(MGIA_cyto[25:32], as.numeric)
cor_MGIA_cyto <-data.frame(correlation(data = MGIA_cyto, select = colnames(MGIA_cyto[25:32]), select2 = "pre vacc", 
                                             method = "spearman",
                                             p_adjust = "none"))

?correlation

str(MGIA_olink)
MGIA_olink[8:80] <- sapply(MGIA_olink[8:80], as.numeric)
cor_MGIA_olink <-data.frame(correlation(data = MGIA_olink, select = colnames(MGIA_olink[8:80]), select2 = "pre vacc", 
                                       method = "spearman",
                                       p_adjust = "none"))



MGIA_meta <- left_join(join_meta_MGIA,olink, by = "PatientID")
str(MGIA_meta)
MGIA_meta[28:100] <- sapply(MGIA_meta[28:100], as.numeric)
cor_MGIA_meta <-data.frame(correlation(data = MGIA_meta, select = colnames(MGIA_meta[28:100]), select2 = colnames(MGIA_meta[5:21]), 
                                        method = "spearman",
                                        p_adjust = "fdr"))


heatmap_corr2_mito <- cor_MGIA_meta %>%
  filter(p < 0.05) %>%
  select(Parameter1, Parameter2, rho) %>%
  pivot_wider(names_from = Parameter1, values_from = rho)

metabolites_list2 <- metabolite_select

heatmap_corr2_mito <- heatmap_corr2_mito%>%
  select(Parameter2, all_of(metabolites_list2)) %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(1) %>%
  rownames_to_column("Parameter1") %>%
  left_join(significant_metabolites_name, by = "Parameter1") %>%
  column_to_rownames("name") %>%
  select(-Parameter1)

heatmap_corr2_mito <- data.frame(heatmap_corr2_mito %>% column_to_rownames("Parameter2"))

matrix_mito <- heatmap_corr2_mito %>%
  map_dfc(as.numeric) %>%
  as.matrix()

matrix_mito[is.na(matrix_mito)] <- 0

rownames(matrix_mito) <- rownames(heatmap_corr2_mito)

rg <- max(abs(matrix_mito))


z <- data.frame(matrix_mito)
z <- rownames_to_column(z,"name")

#pathway <- read.csv("INFECT_metabolic pathway edited.v3.csv")
#pathway <- pathway %>% select(name,pathway_name)
#pathway <- left_join(z,pathway)
#pathway <- pathway %>% select(name,pathway_name)
#pathway$pathway_name <- if_else(is.na(pathway$pathway_name),"",pathway$pathway_name)


mito_heatmap <- matrix_mito %>%
  ComplexHeatmap::pheatmap(scale = "none",
                           name = "Spearman's rho",
                           #labels_col = labels,
                           #labels_row = pathway$pathway_name,
                           clustering_distance_cols = "correlation",
                           clustering_distance_rows = "correlation",
                           cluster_rows = T,
                           cluster_cols = T,
                           #breaks = seq(0.5, -0.5, length.out = 7),
                           breaks = seq(1, -1, length.out = 7),
                           #color = viridis(7),
                           color = colorRampPalette(c("firebrick3", "white", "navy"))(7),
                           #annotation_col = after_underscore,
                           #annotation_colors = anno_colors,
                           fontsize_row = 14,
                           fontsize = 14,
                           #display_numbers = stat_z,
                           #fontsize = 40,
                           #legend_breaks = c(-1,-0.5,0,0.5,1),
                           show_colnames = TRUE,
                           #fontsize_col = 16,
                           #cellwidth = 15,
                           cellwidth = 15,
                           border_color = NA,
                           legend = TRUE,
                           number_color = "black",
                           fontsize_number = 9,
                           angle_col = "45"
                           #annotation_legend = FALSE,
                           #annotation = FALSE,
                           #column_split=after_underscore$Stimuli
                           #legend_labels = c(-1,-0.5,0,0.5,1),
                           #cutree_rows = 3,
                           #clustering_distance_rows = "euclid",
                           #cutree_cols = 3,
                           #,left_annotation = rowAnnotation(foo = anno_text(pathway$name,
                           #                                                 #z1$Source,
                           #                                                 location = 0.99,just = "right",
                           #                                                 gp = gpar(fontsize = 14)))
  )

mito_heatmap


