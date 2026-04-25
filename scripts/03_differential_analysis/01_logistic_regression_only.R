# load packages -----------------------------------------------------------
library(pacman)
pacman::p_load(sjPlot,janitor,stringi,epiDisplay,dplyr,openxlsx,data.table,tableone,tidyverse,ggplot2,viridis,readxl,EnhancedVolcano,pheatmap,RColorBrewer,gplots)
library(sjPlot)
library(janitor)
library(stringi)
library(epiDisplay)
library(dplyr)
library(openxlsx)
library(data.table)
library(tidyverse)
library(tableone)
library(ggplot2)
library(viridis)

# Load data from extracted CSVs -------------------------------------------
source("../00_load_data.R")

# Set output directory for results
dir <- OUTPUT_DIR

# Create output directories if needed
create_output_dirs()

# Load metabolite name mapping
name <- read.csv(METABOLITE_NAME_MAP)
colnames(name) <- c("metabolites","name")

# Additional derived variables
klinis$Mtb.strain2 <- ifelse(klinis$Mtb.strain=="East-Asian (Beijing)","1","0")

# change data types -----------------------------------------------------
test<-klinis$status =="Converter"
klinis$status[test] <- "0"
test<-klinis$status=="Persistently_uninfected"
klinis$status[test] <- "1"
test<-klinis$DM =="2"
klinis$DM[test] <- "1"
test<-klinis$DM =="3"
klinis$DM[test] <- "2"

smoking_map <- c("Never smoked" = "1", "Current smoker" = "2", "Quit >6 mo ago" = "3")
klinis$smoking <- smoking_map[as.character(klinis$smoking)]
klinis$smoking <- sapply(klinis$smoking,as.ordered)
klinis$smoking <- ordered(klinis$smoking,levels=c("1","2","3"))

klinis$hours_with_caseindex <- as.numeric(klinis$hours_with_caseindex)
klinis$status <- as.numeric(klinis$status)
klinis$extent <- as.numeric(klinis$extent)
klinis$extent <- klinis$extent/100
klinis$IdInfect <- sapply(klinis$IdInfect,as.character)
klinis$DM <- sapply(klinis$DM,as.factor)

test<-klinis$strict_0.15_conversion =="Converter"
klinis$strict_0.15_conversion[test] <- "0"
test<-klinis$strict_0.15_conversion=="Persistently_uninfected"
klinis$strict_0.15_conversion[test] <- "1"
test<-klinis$strict_0.15_conversion=="Uncertain"
klinis$strict_0.15_conversion[test] <- "2"

klinis$cavity <- as.factor(klinis$cavity)
klinis$IdInfect <- sapply(klinis$IdInfect,as.character)

factor_cols <- c("sex_contact","BCG_scar","DM","smoking","Mtb.strain","cavity")
factor_cols <- factor_cols[factor_cols %in% colnames(klinis)]
if (length(factor_cols) > 0) klinis[,factor_cols] <- lapply(klinis[,factor_cols], as.factor)

ordered_cols <- c("smoking")
ordered_cols <- ordered_cols[ordered_cols %in% colnames(klinis)]
if (length(ordered_cols) > 0) klinis[,ordered_cols] <- lapply(klinis[,ordered_cols], ordered)

numeric_cols <- c("NIL1","TBAg.NIL1","Mitogen.Nil1","NIL1_14","TBAg.NIL1_14","Mitogen.Nil1_14")
numeric_cols <- numeric_cols[numeric_cols %in% colnames(klinis)]
if (length(numeric_cols) > 0) klinis[,numeric_cols] <- lapply(klinis[,numeric_cols], as.numeric)

# Remove variables with >50% missing
klinis_numeric <- as.data.frame(lapply(klinis, as.numeric))
merge_numeric <- cbind(klinis_numeric,metabolites)

missing <- klinis %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.))))
missing <- as.data.frame(t(missing))
missing$percent <- missing$V1/nrow(klinis)
missing <- missing %>% filter(percent >0.5)
missing_var <- c(rownames(missing))

klinis_numeric <- klinis_numeric %>% select(-missing_var)
klinis <- klinis %>% select(-missing_var)

nama.variabel <- colnames(klinis)
nama.variabel <- nama.variabel [! nama.variabel %in% c("IdInfect","IdPatient")]

klinis_2 <- klinis %>% filter(strict_0.15_conversion %in% c("0","1"))

klinis_2$strict_0.15_conversion <- as.character(klinis_2$strict_0.15_conversion)
klinis_2$strict_0.15_conversion <- as.factor(klinis_2$strict_0.15_conversion)

klinis_strain <- klinis %>% filter(!is.na(Mtb.strain2))

# Table 1 output
TableOne <- CreateTableOne(vars = nama.variabel,strata=c("status"),data=klinis)
Stat_metabolomic<- print(TableOne, nonnormal = TRUE)
write.csv(Stat_metabolomic,file.path(dir,"stat_metabolomic.csv"))

TableOne2 <- CreateTableOne(vars = nama.variabel,strata="strict_0.15_conversion",data=klinis_2)
Stat_metabolomic<- print(TableOne2, nonnormal = TRUE)
write.csv(Stat_metabolomic,file.path(dir,"stat_metabolomic_without_missing.csv"))

TableOne2 <- CreateTableOne(vars = nama.variabel,strata="status",data=klinis_strain)
Stat_metabolomic<- print(TableOne2, nonnormal = TRUE)
write.csv(Stat_metabolomic,file.path(dir,"stat_metabolomic_strain.csv"))

nama.metabolites <- colnames(metabolites)

########### merge for strict #########################
merge_strict1 <- merge %>% filter(strict_0.15_conversion == "Converter")
merge_strict2 <- merge %>% filter(strict_0.15_conversion == "Persistently_uninfected")
merge_strict <- rbind(merge_strict1,merge_strict2)
rownames(merge_strict) <- merge_strict$IdInfect
merge_strict <- arrange(merge_strict,IdInfect)
rownames(merge_strict) <- merge_strict$IdInfect
met_strict <- merge_strict %>% select(nama.metabolites)
klinis_strict <- merge_strict %>% select(!colnames(metabolites))

test<-klinis_strict$status =="Converter"
klinis_strict$status[test] <- "0"
test<-klinis_strict$status =="Persistently_uninfected"
klinis_strict$status[test] <- "1"

test<-klinis_strict$DM =="2"
klinis_strict$DM[test] <- "1"
test<-klinis_strict$DM =="3"
klinis_strict$DM[test] <- "2"

klinis_strict$smoking <- smoking_map[as.character(klinis_strict$smoking)]
klinis_strict$smoking <- sapply(klinis_strict$smoking,as.ordered)
klinis_strict$smoking <- ordered(klinis_strict$smoking,levels=c("1","2","3"))

klinis_strict$hours_with_caseindex <- as.numeric(klinis_strict$hours_with_caseindex)
klinis_strict$status <- as.numeric(klinis_strict$status)
klinis_strict$extent <- as.numeric(klinis_strict$extent)
klinis_strict$extent <- klinis_strict$extent/100
klinis_strict$IdInfect <- sapply(klinis_strict$IdInfect,as.character)
klinis_strict$DM <- sapply(klinis_strict$DM,as.factor)

if ("strict_0.15_conversion" %in% colnames(klinis_strict)) {
  test <- klinis_strict$strict_0.15_conversion == "Converter"
  klinis_strict$strict_0.15_conversion[test] <- "0"
  test <- klinis_strict$strict_0.15_conversion == "Persistently_uninfected"
  klinis_strict$strict_0.15_conversion[test] <- "1"
  klinis_strict$strict_0.15_conversion <- as.numeric(klinis_strict$strict_0.15_conversion)
}

klinis_strict$smear_result_3 <- as.numeric(klinis_strict$smear_result_3)
klinis_strict$extent <- klinis_strict$extent/100
str(klinis_strict)
klinis_strict$IdInfect <- sapply(klinis_strict$IdInfect,as.character)
klinis_strict$DM <- sapply(klinis_strict$DM,as.factor)

if (is.character(klinis_strict$strict_0.15_conversion)) {
  test<-klinis_strict$strict_0.15_conversion=="Converter"
  klinis_strict$strict_0.15_conversion[test] <- "0"
  test<-klinis_strict$strict_0.15_conversion=="Persistently_uninfected"
  klinis_strict$strict_0.15_conversion[test] <- "1"
}

klinis_strict$cavity <- as.factor(klinis_strict$cavity)
klinis_strict$IdInfect <- sapply(klinis_strict$IdInfect,as.character)

factor_cols_strict <- c("status","sex_contact","BCG_scar","DM","smoking","Mtb.strain","cavity")
factor_cols_strict <- factor_cols_strict[factor_cols_strict %in% colnames(klinis_strict)]
if (length(factor_cols_strict) > 0) klinis_strict[,factor_cols_strict] <- lapply(klinis_strict[,factor_cols_strict], as.factor)

ordered_cols_strict <- c("smoking")
ordered_cols_strict <- ordered_cols_strict[ordered_cols_strict %in% colnames(klinis_strict)]
if (length(ordered_cols_strict) > 0) klinis_strict[,ordered_cols_strict] <- lapply(klinis_strict[,ordered_cols_strict], ordered)

numeric_cols_strict <- c("NIL1","TBAg.NIL1","Mitogen.Nil1","NIL1_14","TBAg.NIL1_14","Mitogen.Nil1_14")
numeric_cols_strict <- numeric_cols_strict[numeric_cols_strict %in% colnames(klinis_strict)]
if (length(numeric_cols_strict) > 0) klinis_strict[,numeric_cols_strict] <- lapply(klinis_strict[,numeric_cols_strict], as.numeric)

common_cols <- intersect(colnames(klinis), colnames(klinis_strict))
klinis_strict <- klinis_strict %>% select(all_of(common_cols))
z <- colnames(klinis_strict)
vars <- z[-(1:2)]

x <- CreateTableOne(vars = vars, data = klinis_strict, strata = "status")
x <- print(x, nonnormal = TRUE)
write.csv(x, file.path(OUTPUT_DIR, "karakteristik_strict.csv"))

write.csv(klinis_strict, file.path(OUTPUT_DIR, "2.klinis_strict_metabolomics.csv"))

################## logistic regression all samples ####################
klinis <- arrange(klinis,desc(IdInfect))
metabolites <- rownames_to_column(metabolites,var="IdInfect")
metabolites <- arrange(metabolites,desc(IdInfect))
metabolites <- column_to_rownames(metabolites,var="IdInfect")

result.lin.mod<- NULL
za <- NULL

for (x in nama.metabolites) {
  lin.mod <- glm(klinis$status ~ metabolites[,x]+klinis$agecontact+klinis$sex_contact+klinis$exp_score+klinis$BCG_scar, family = "binomial",na.action = na.exclude)
  summary.lin.mod<-summary(lin.mod)
  OR <- epiDisplay::logistic.display(lin.mod)
  OR <- OR$table
  OR <- OR[1,1:3]
  p.val.lin.mod<-summary.lin.mod$coefficients[2,1:4]
  result.lin.mod <- rbind(result.lin.mod,c(x,p.val.lin.mod,OR))
}

summary.lin.mod

######################### volcano fullset ########################

sig.metabolites <- data.frame(rbind(result.lin.mod))
sig.metabolites$p.adj = p.adjust(sig.metabolites$Pr...z.., method = "BH")

sig.metabolites <- sig.metabolites %>%
  dplyr::rename(metabolites = V1)

sig.metabolites <- left_join(sig.metabolites,name)

str(sig.metabolites)
sig.metabolites$Pr...z.. <- as.numeric(sig.metabolites$Pr...z..)
sig.metabolites$Estimate <- as.numeric(sig.metabolites$Estimate)
sig.metabolites$z.value <- as.numeric(sig.metabolites$z.value)
sig.metabolites$OR <- as.numeric(sig.metabolites$OR)
sig.metabolites$log2OR <- log(sig.metabolites$OR,2)
converter_metabolites <- sig.metabolites %>%  filter(sig.metabolites$Estimate < 0)
converter_metabolites <- converter_metabolites %>% arrange (converter_metabolites$Pr...z..)
converter_metabolites <- converter_metabolites %>% filter (converter_metabolites$Pr...z.. < 0.05)
write.csv(converter_metabolites,file.path(OUTPUT_DIR, "3. converter metabolites_fullset.csv"))

PU_metabolites <- sig.metabolites %>%  filter(sig.metabolites$Estimate > 0)
PU_metabolites <- PU_metabolites %>% arrange (PU_metabolites$Pr...z..)
PU_metabolites <- PU_metabolites %>% filter (PU_metabolites$Pr...z.. < 0.05)
PU_metabolites_HMDB <- c(PU_metabolites$metabolites)
PU_metabolites_name <- c(PU_metabolites$name)
write.csv(PU_metabolites,file.path(OUTPUT_DIR, "3. PU_metabolites_fullset.csv"))

log.res.metabolites <- sig.metabolites  %>% filter (sig.metabolites$Pr...z.. < 0.05)
fullset_sig <- c(log.res.metabolites$metabolites)
write.csv(log.res.metabolites,file.path(dir,"3. significant metabolites_fullset.csv"))
str(sig.metabolites)

keyvals.colour <- ifelse(
  sig.metabolites$`Pr...z..` > 0.05, 'black',
  ifelse(sig.metabolites$`Pr...z..` < 0.05 & sig.metabolites$Estimate > 0, 'blue',
         ifelse(sig.metabolites$`Pr...z..` < 0.05 & sig.metabolites$Estimate < 0,'red',
                'black')))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Higher in IGRA converter'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Higher in Persistently IGRA-negatives'

plot <- EnhancedVolcano(sig.metabolites,
                        lab = sig.metabolites$name,
                        titleLabSize = 1,
                        subtitleLabSize = 1,
                        subtitle = "",
                        x = 'Estimate',
                        y =  'Pr...z..',
                        legendPosition = "top",
                        legendLabSize = 15,
                        legendIconSize = 5,
                        legendDropLevels = FALSE,
                        ylab = bquote(~italic("P value")),
                        xlim = c(-3,3),
                        ylim = c(0, 3.5),
                        colCustom = keyvals.colour,
                        title = '',
                        pCutoff = 0.05,
                        FCcutoff = 0,
                        pointSize = 5.0,
                        labSize = 5.25,
                        colAlpha = 0.3,
                        gridlines.minor = FALSE,
                        xlab = "Logistic Regression Coefficient Estimate",
                        drawConnectors = FALSE,
                        caption =NULL,
                        widthConnectors = 1)

plot2 <- plot +
  scale_y_continuous(limits = c(0,3), breaks = c(0,1,1.30103,2,3),labels = c(0,0.1,0.05,0.01,0.001))

plot2
write.csv(sig.metabolites,file.path(REVISION_DIR, "S1A.csv"))
ggsave(file.path(FIGURE_DIR, "figure4_full.v4.png"), plot2, width = 11, height = 10, dpi = 600)

################## logistic regression strict selecction##################################################
klinis_strict <- arrange(klinis_strict,desc(IdInfect))
met_strict <- rownames_to_column(met_strict,var="IdInfect")
met_strict <- arrange(met_strict,desc(IdInfect))
met_strict <- column_to_rownames(met_strict,var="IdInfect")

result.lin.mod<- NULL
za <- NULL

for (x in nama.metabolites) {
  lin.mod <- glm(klinis_strict$strict_0.15_conversion ~ met_strict[,x]+klinis_strict$agecontact+klinis_strict$sex_contact + klinis_strict$BMI_contact + klinis_strict$BCG_scar + klinis_strict$exp_score, family = "binomial",na.action = na.exclude)
  summary.lin.mod<-summary(lin.mod)
  OR <- epiDisplay::logistic.display(lin.mod)
  OR <- OR$table
  OR <- OR[1,1:3]
  p.val.lin.mod<-summary.lin.mod$coefficients[2,1:4]
  result.lin.mod <- rbind(result.lin.mod,c(x,p.val.lin.mod,OR))
}

summary.lin.mod

##============ volcano logistic regrestion strict dataset===============================

sig.metabolites_sub <- data.frame(rbind(result.lin.mod))
sig.metabolites_sub$p.adj = p.adjust(sig.metabolites_sub$Pr...z.., method = "BH")

sig.metabolites_sub <- sig.metabolites_sub %>%
  dplyr::rename(metabolites = V1)

sig.metabolites_sub <- left_join(sig.metabolites_sub,name)

str(sig.metabolites_sub)
sig.metabolites_sub$Pr...z.. <- as.numeric(sig.metabolites_sub$Pr...z..)
sig.metabolites_sub$Estimate <- as.numeric(sig.metabolites_sub$Estimate)
sig.metabolites_sub$z.value <- as.numeric(sig.metabolites_sub$z.value)
sig.metabolites_sub$OR <- as.numeric(sig.metabolites_sub$OR)
sig.metabolites_sub$log2OR <- log(sig.metabolites_sub$OR,2)

write.csv(sig.metabolites_sub,file.path(OUTPUT_DIR, "3. result logres metabolites_subset.csv"))

converter_metabolites <- sig.metabolites_sub %>%  filter(sig.metabolites_sub$Estimate < 0)
converter_metabolites <- converter_metabolites %>% arrange (converter_metabolites$Pr...z..)
converter_metabolites <- converter_metabolites %>% filter (converter_metabolites$Pr...z.. < 0.05)
write.csv(converter_metabolites,file.path(OUTPUT_DIR, "3. converter metabolites_subset.csv"))

PU_metabolites <- sig.metabolites_sub %>%  filter(sig.metabolites_sub$Estimate > 0)
PU_metabolites <- PU_metabolites %>% arrange (PU_metabolites$Pr...z..)
PU_metabolites <- PU_metabolites %>% filter (PU_metabolites$Pr...z.. < 0.05)
PU_metabolites_HMDB <- c(PU_metabolites$metabolites)
PU_metabolites_name <- c(PU_metabolites$name)
write.csv(PU_metabolites,file.path(OUTPUT_DIR, "3. PU_metabolites_subset.csv"))

log.res.metabolites_sub <- sig.metabolites_sub  %>% filter (sig.metabolites_sub$Pr...z.. < 0.05)
strictset_sig <- c(log.res.metabolites_sub$metabolites)
write.csv(log.res.metabolites_sub,file.path(dir,"3. significant_metabolites_subset.csv"))
str(sig.metabolites_sub)

z <- log.res.metabolites_sub

keyvals.colour <- ifelse(
  sig.metabolites_sub$`Pr...z..` > 0.05, 'black',
  ifelse(sig.metabolites_sub$`Pr...z..` < 0.05 & sig.metabolites_sub$Estimate > 0, 'blue',
         ifelse(sig.metabolites_sub$`Pr...z..` < 0.05 & sig.metabolites_sub$Estimate < 0,'red',
                'black')))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Higher in\nIGRA converters'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Higher in\nPersistently IGRA-negatives'

plot <- EnhancedVolcano(sig.metabolites_sub,
                lab = sig.metabolites_sub$name,
                titleLabSize = 1,
                subtitleLabSize = 1,
                subtitle = "",
                x = 'Estimate',
                y =  'Pr...z..',
                legendPosition = "top",
                legendLabSize = 19,
                legendIconSize = 5,
                legendDropLevels = FALSE,
                ylab = bquote(~ italic("P value")),
                xlim = c(-3,3),
                ylim = c(0, 4),
                colCustom = keyvals.colour,
                title = '',
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 5.0,
                labSize = 5.25,
                colAlpha = 0.3,
                gridlines.minor = FALSE,
                xlab = "Logistic Regression Coefficient Estimate",
                drawConnectors = FALSE,
                caption = NULL,
                widthConnectors = 1)

plot2 <- plot +
  scale_y_continuous(limits = c(0,3), breaks = c(0,1,1.30103,2,3),labels = c(0,0.1,0.05,0.01,0.001))

plot2
write.csv(sig.metabolites_sub,file.path(REVISION_DIR, "2A.csv"))
ggsave(file.path(FIGURE_DIR, "figure5_subset_pval.v10.png"), plot2, width = 10, height = 10, dpi = 600)

############# heatmap of significant metabolites from the two sets ##################
library(gplots)
combined <- c(fullset_sig,strictset_sig)
combined
intersection <- intersect(fullset_sig,strictset_sig)

combined_x <- setdiff(combined, intersection)

strictset_sig_select_x <- sig.metabolites_sub %>% filter(metabolites %in% strictset_sig)
strictset_sig <- c(strictset_sig_select_x$metabolites)
strictset_sig_select <- sig.metabolites_sub %>% filter(metabolites %in% combined)

fullset_sig_select <- sig.metabolites %>% filter(metabolites %in% combined)

strictset_sig_select2 <- strictset_sig_select %>% select(metabolites,name,Estimate)
fullset_sig_select2 <- fullset_sig_select %>% select(metabolites,name,Estimate)

combined_sig_select <- left_join(strictset_sig_select2,fullset_sig_select2, by="metabolites")
combined_sig_select_x <- combined_sig_select %>% filter(metabolites %in% combined_x)
combined_sig_select_intersect <- combined_sig_select %>% filter(metabolites %in% intersection)

combined_sig_select <- combined_sig_select %>% select(metabolites,name.x,Estimate.x,Estimate.y)
combined_sig_select_x <- combined_sig_select_x %>% select(metabolites,name.x,Estimate.x,Estimate.y)
combined_sig_select_intersect <- combined_sig_select_intersect %>% select(metabolites,name.x,Estimate.x,Estimate.y)

combined_sig_select <- column_to_rownames(combined_sig_select, var="name.x")
combined_sig_select_x <- column_to_rownames(combined_sig_select_x, var="name.x")
combined_sig_select_intersect <- column_to_rownames(combined_sig_select_intersect, var="name.x")

colnames(combined_sig_select) <- c("HMDB","subset_selection","fullset")
colnames(combined_sig_select_x) <- c("HMDB","subset_selection","fullset")
colnames(combined_sig_select_intersect) <- c("HMDB","subset_selection","fullset")

combined_sig_select$sig.fullset <- ifelse(combined_sig_select$HMDB %in% fullset_sig, "YES","NO")
combined_sig_select$sig.subset <- ifelse(combined_sig_select$HMDB %in% strictset_sig, "YES","NO")
combined_sig_select <- arrange(combined_sig_select,sig.fullset,sig.subset,desc(fullset))

combined_sig_select_matrix <- as.matrix(combined_sig_select %>% select(subset_selection,fullset))
combined_sig_select_x_matrix  <- as.matrix(combined_sig_select_x %>% select(subset_selection,fullset))
combined_sig_select_intersect_matrix  <- as.matrix(combined_sig_select_intersect %>% select(subset_selection,fullset))
str(combined_sig_select_intersect)

selected_metabolites <- as.data.frame(combined_sig_select)
selected_metabolites$selected <- ifelse(!between(selected_metabolites$subset_selection,-0.3,0.3)&!between(selected_metabolites$fullset,-0.3,0.3),"YES","NO")
selected_metabolites <- selected_metabolites %>% filter(selected_metabolites$selected == "YES")
write.csv(selected_metabolites,file.path(dir,"selected_metabolites.csv"))

list_selected_metabolites <- rownames(selected_metabolites)
HMDB_selected_metabolites <- selected_metabolites$HMDB

coul <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100)
coul2 <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(256))

col_annotation <- as.data.frame(combined_sig_select %>% select(sig.fullset,sig.subset))
colnames(col_annotation) <- c("Significant in full set only","Significant in subset only")
test <- col_annotation
test$test <- if_else(test$`Significant in full set only`=="YES" & test$`Significant in subset only`=="NO","Significant in full set only",if_else(test$`Significant in full set only`=="NO" & test$`Significant in subset only`=="YES","Significant in subset only","Significant in both datasets"))

Fullset        <- c("#FDE725FF","#35B779FF","#31688EFF")
names(Fullset) <- c("Significant in full set only","Significant in subset only","Significant in both datasets")

combined_sig_select_matrix_2 <- as.data.frame(combined_sig_select_matrix)
combined_sig_select_matrix_2 <- combined_sig_select_matrix_2 %>% select(subset_selection)
combined_sig_select_matrix_2 <- rownames_to_column(combined_sig_select_matrix_2,var="meta")
combined_sig_select_matrix_2 <- combined_sig_select_matrix_2 %>% filter(meta %in% strictset_sig)
combined_sig_select_matrix_2 <- column_to_rownames(combined_sig_select_matrix_2, var="meta")

anno_colors2 <- list(`test` = Fullset)

test <- test %>% select(test)
anno_colors <- list(test = Fullset)
rg <- max(abs(combined_sig_select_matrix), na.rm = TRUE);
if(!is.finite(rg) || rg == 0) rg <- 1;

heatmap_output <- combined_sig_select_matrix  %>% 
  pheatmap(scale = "none",
           clustering_distance_cols = "correlation",
           cluster_rows = F,
           cluster_cols = F,
           breaks = seq(-rg, rg, length.out = 100),
           color = redblue(100),
           annotation_row = test,
           annotation_colors = anno_colors,
           angle_col = 0,
           cellwidth = 100,
           cellheight = 26.5,
           fontsize_row = 22,
           fontsize_col = 22,
           border_color = "white",
           labels_col = c("Subset","Full set"),
           anno_width = 70,
           annotation_names_legend = FALSE,
           annotation_names_row = FALSE
  )

heatmap_output

write.csv(combined_sig_select,file.path(REVISION_DIR, "2B.CSV"))
ggsave(file.path(FIGURE_DIR, "figure8_only subset metabolites.v6.png"), heatmap_output, width = 20, height = 30, dpi = 600)

cat("Logistic regression analysis complete.\n")
cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Figures saved to:", FIGURE_DIR, "\n")
