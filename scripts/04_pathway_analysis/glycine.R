library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pathview)
# dataset -----------------------------------------------------------------

merge = read_rds("INFECT_metabolomics_dataset.rds")

# pathway -----------------------------------------------------------------

gly_pathway <- keggGet("hsa00260")

gly_compounds <- names(gly_pathway[[1]]$COMPOUND)

gly_hmdb = read.csv("name_map_gly.csv")

matching_cols <- colnames(merge)[colnames(merge) %in% gly_hmdb$HMDB]

merge_gly = merge |> 
  select(all_of(matching_cols), strict_0.15_conversion, IdInfect)


# Rename columns in your dataset
merge_gly_renamed <- merge_gly

metabolite_cols <- gly_hmdb$Match[match(names(merge_gly_renamed)[names(merge_gly_renamed) %in% gly_hmdb$HMDB], 
                                        gly_hmdb$HMDB)]

names(merge_gly_renamed)[names(merge_gly_renamed) %in% gly_hmdb$HMDB] <- 
  gly_hmdb$Match[match(names(merge_gly_renamed)[names(merge_gly_renamed) %in% gly_hmdb$HMDB], 
                       gly_hmdb$HMDB)]


# comparison --------------------------------------------------------------

# 2. Load libraries
library(pathview)
library(dplyr)
library(tidyr)


# 4. Prepare data: Calculate Log2 Fold Change (L2FC)
# We will compare 'Persistently_uninfected' vs 'Uncertain'
df_filtered <- merge_gly_renamed %>%
  filter(strict_0.15_conversion %in% c("Persistently_uninfected", "Converter"))

# Calculate mean levels for each group and handle potential 0s with a small pseudocount
pseudocount <- 0.01
mean_levels <- df_filtered %>%
  group_by(strict_0.15_conversion) %>%
  summarise(across(all_of(metabolite_cols), ~ mean(.x + pseudocount, na.rm = TRUE))) %>%
  pivot_longer(cols = -strict_0.15_conversion, names_to = "Metabolite", values_to = "Mean_Level") %>%
  pivot_wider(names_from = strict_0.15_conversion, values_from = Mean_Level)


l2fc_data <- mean_levels %>%
  mutate(
    L2FC = (`Converter` - `Persistently_uninfected`)
  )

kegg_id_map <- gly_hmdb |> 
  filter(gly_hmdb$Match %in% l2fc_data$Metabolite)

kegg_id_map = kegg_id_map |> 
  rename(Metabolite = Match)

l2fc_data = left_join(l2fc_data,kegg_id_map,"Metabolite")
logres = read.csv("3. all logres metabolites_subset.csv")

logres = logres |> 
  rename(HMDB = metabolites)

l2fc_data = left_join(l2fc_data,logres, "HMDB")

pathview_data <- l2fc_data$Estimate
names(pathview_data) <- l2fc_data$Query

# Remove any metabolites that didn't map (if any)
pathview_data <- pathview_data[!is.na(names(pathview_data))]

print("Data vector prepared for pathview:")
print(pathview_data)

# 7. Run pathview
# This saves a diagram file (PNG and PDF) to your current R working directory
pv_out <- pathview(
  cpd.data = pathview_data,
  low = "red",
  high = "blue",
  pathway.id = "hsa00260",     # KEGG ID for "Arachidonic acid metabolism"
  species = "hsa",           # 'hsa' for human. Change if your data is from another species (e.g., 'mmu' for mouse)
  out.suffix = "L2FC_gly_pathway",
  cpd.legend.name = "Log2 Fold Change\n(IGRA Converter / Persistently-IGRA negatives)",
  kegg.native = TRUE         # Use native KEGG diagram style
)

print(paste("Pathway visualization saved as:", pv_out$plot.data.cpd))
