library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(dplyr)
library(conflicted)
library(patchwork)
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
# install.packages("extrafont")
library(extrafont)
# Load the fonts
loadfonts(device = "win")  # for Windows

theme_set(theme_pubr(base_family = "Helvetica"))


# Set default theme once
my_theme <- theme_pubr() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, family = "Helvetica", size = 20),
    axis.title.y = element_text(family = "Helvetica", margin = margin(r = 10)),
    axis.text = element_text(color = "black", family = "Helvetica"),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 0, vjust = 0.5)  # Ensure proper alignment for two-line text
  )

theme_set(my_theme)


# Function to create standardized metabolite plot
create_metabolite_plot <- function(data, metabolite_id, title) {
  
  #metabolite_id = "HMDB0000929"
  
  
  # Create formula using the actual column name
  formula_str <- paste("`status` ~", paste("`", metabolite_id, "` + agecontact + sex_contact + exp_score + BMI_contact + BCG_scar", sep=""))
  
  # Calculate statistics
  model <- glm(as.formula(formula_str), data = data, family = "binomial")
  model_stats <- summary(model)
  # Calculate y position with more space for annotation
  y_range <- range(data[[metabolite_id]], na.rm = TRUE)
  y_span <- diff(y_range)
  
  stats <- list(
    beta = round(coef(model)[2], 3),
    OR = round(exp(coef(model)[2]), 3),
    p_val = round(coef(model_stats)[2,4], 3),
    y_pos = y_range[2] + (y_span * 0.2) 
  )

data = data |> 
    mutate(status = factor(status,
                           levels = c("Persistently\nIGRA-negatives", "IGRA-converters")
    )
    )
  
  
  ggboxplot(data,
            x = "status",
            y = metabolite_id,
            width = 0.25,
            color = "status",
            palette = c("black", "black"),
            add = "jitter",
            add.params = list(color = "black", alpha = 0.25)
  ) +
    geom_violin(aes(fill = status), alpha = 0.2) +
    scale_fill_manual(values = c("#3366CC", "#CC3333")) +
    my_theme +
    scale_y_continuous(limits = c(y_range[1], y_range[2] + (y_span * 0.3))) +  # Add 30% of range to top
    labs(
      title = title,
      x = "",
      y = "Log2 (peak ion Intensity)"
    ) +
    annotate(
      "text",
      x = 1.5,  # Center between groups
      y = stats$y_pos,
      label = paste("OR of persistently IGRA-negatives =", stats$OR, "\n",
                    "P =", stats$p_val),
      hjust = 0.5,
      vjust = 0
    )
}

clin <- read_rds('INFECT_clinical data_1347.rds')

phenotype = clin |> filter(status %in% c("Persistently IGRA-negatives", "IGRA converters"))

# Data preparation with two-line label
data_plot <- phenotype %>%
#   mutate(status = factor(status,
#                          levels = c("Persistently_uninfected", "Converter"),
#                          labels = c("Persistently\nIGRA-negatives", "IGRA-converters")  # Added line break
#   ))
  mutate(status = factor(status,
                         levels = c("IGRA converters", "Persistently IGRA-negatives")
  )) #%>%
  #filter(strict_0.15_conversion != "Uncertain")

data_plot$sex = factor(data_plot$sex,
                               levels = c("Female","Male"))


dir_gen = "D:\\Download\\General\\genotype\\"
chr13 = read_rds(paste0(dir_gen,"chr13_TBM_INFECT_imputed_EAS_IND.finalQC.dosage.rds"))
chr3 = read_rds(paste0(dir_gen,"chr3_TBM_INFECT_imputed_EAS_IND.finalQC.dosage.rds"))
chr22 = read_rds(paste0(dir_gen,"chr22_TBM_INFECT_imputed_EAS_IND.finalQC.dosage.rds"))


data_plot = data_plot |> 
  dplyr::select(
    IdInfect,
    status
  )

data_plot$IdInfect = as.character(data_plot$IdInfect)
str(data_plot)


chr13_CLYBL = chr13[rownames(chr13) == "13:100270229:C:T", ] 
rownames(chr13_CLYBL) = "rs2793763"
str(chr13_CLYBL)

chr13_CLYBL = as.data.frame(t(chr13_CLYBL)) |> 
  rownames_to_column("IdInfect") |> 
  mutate(IdInfect = gsub("X","",IdInfect)) |> 
  mutate(
    `rs2793763` = if_else(
      `rs2793763` < 0.5,
      "CC", if_else(
        `rs2793763` > 1.5,
        "TT",
        "CT"
      )
    )
    ) |> 
  mutate(
    `rs2793763` = factor(
      `rs2793763`,
      levels = c("CC","CT","TT")
    )
  )

data_plot$IdInfect = as.character(data_plot$IdInfect)
chr13_CLYBL = left_join(data_plot,
                        chr13_CLYBL,
                        "IdInfect")

library(tidyplots)
chr13_CLYBL = chr13_CLYBL |> filter(!is.na(rs2793763))
chr13_CLYBL |> 
  tidyplot(x=rs2793763, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("maroon", "skyblue")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

gtsummary::tbl_summary(chr13_CLYBL |> dplyr::select(status,rs2793763), rs2793763)

ggsave("Fig 1E. INFECT genetic_ECvsCon_proportion.png")






chr22_PRODH = chr22_PRODH |> 
  filter(!is.na(`rs13056032`))

ggpubr::ggboxplot(
  chr22_PRODH,
  title = "Proline",
  x = "rs13056032",
  width = 0.6,
  y = "HMDB0000162",
  ylab = "Log2 (peak ion intensity)",
  add = c("jitter"),
  add.params = list(color = "status", alpha = 0.3)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c(1,2), c(1,3), c(2,3)),  # Pairwise comparisons
    method = "wilcox.test",
    label = "p.format"
  )

chr22_PRODH |> 
  tidyplot(x=rs13056032, y=HMDB0000162, dodge_width = 0) |>
  add_boxplot() |>
  add_data_points_beeswarm() |> 
  adjust_legend_position("top")





chr22_PRODH_v2 = chr22_PRODH |> 
  mutate(
    `rs13056032` = factor(
      `rs13056032`,
      levels = c("TT","TC","CC"),
      labels = c("TT","TC+CC","TC+CC")
    )
  )

ggpubr::ggboxplot(
  chr22_PRODH_v2,
  title = "Proline",
  x = "rs13056032",
  width = 0.6,
  y = "HMDB0000162",
  ylab = "Log2 (peak ion intensity)",
  add = c("jitter"),
  add.params = list(alpha = 0.3)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c(1,2)),
    method = "wilcox.test",
    label = "p.format"
  )


library(tidyplots)

chr22_PRODH_v2 |> 
  tidyplot(x=rs13056032, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" , "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

ggsave("PRODH proportion_dom.png")

write.csv(chr22_PRODH,"chr22_PRODH.csv")


# proline second rsid -----------------------------------------------------


proline = data_plot |> 
  select(
    IdInfect,
    status,
    HMDB0000929
  )


chr22_PRODH = chr22[rownames(chr22) == "22:18938190:G:A", ] 
rownames(chr22_PRODH) = "rs13054729"
str(chr22_PRODH)

#chr22_PRODH = chr22_PRODH |> 
#  filter(!is.na(`rs13054729`))


chr22_PRODH = as.data.frame(t(chr22_PRODH)) |> 
  rownames_to_column("IdInfect") |> 
  mutate(IdInfect = gsub("X","",IdInfect)) |> 
  mutate(
    `rs13054729` = if_else(
      `rs13054729` < 0.5,
      "GG", if_else(
        `rs13054729` > 1.5,
        "AA",
        "GA"
      )
    )
  ) |> 
  mutate(
    `rs13054729` = factor(
      `rs13054729`,
      levels = c("GG","GA","AA")
    )
  )

#chr22_PRODH_v2 = chr22_PRODH |> 
#  mutate(
#    `rs13054729` = factor(
#      `rs13054729`,
#      levels = c("TT","TC","CC"),
#      labels = c("TT","TC+CC","TC+CC")
#    )
#  )
#

ggpubr::ggboxplot(
  chr22_PRODH,
  title = "Proline",
  x = "rs13054729",
  width = 0.6,
  y = "HMDB0000162",
  ylab = "Log2 (peak ion intensity)",
  add = c("jitter"),
  add.params = list(alpha = 0.3)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c(1,2)),
    method = "wilcox.test",
    label = "p.format"
  )

library(tidyplots)

chr22_PRODH_v2 |> 
  tidyplot(x=rs13056032, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" , "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

gtsummary::tbl_summary(chr22_PRODH_v2 |> select(
  status, rs13056032
), by = "rs13056032") |> 
  add_p()

?add_p()

ggsave("PRODH proportion_dom.png")

write.csv(chr22_PRODH,"chr22_PRODH.csv")

# tryptophan --------------------------------------------------------------

tryptophan = data_plot |> 
  select(
    IdInfect,
    status,
    HMDB0000929
  )

tryptophan$IdInfect = as.character(tryptophan$IdInfect)
str(tryptophan)

chr3_CA = chr3[rownames(chr3) == "3:53557547:T:C", ] 
rownames(chr3_CA) = "rs77915569"
str(chr3_CA)

chr3_CA = as.data.frame(t(chr3_CA)) |> 
  rownames_to_column("IdInfect") |> 
  mutate(IdInfect = gsub("X","",IdInfect)) |> 
  mutate(
    `rs77915569` = if_else(
      `rs77915569` < 0.5,
      "TT", if_else(
        `rs77915569` > 1.5,
        "CC",
        "TC"
      )
    )
  ) |> 
  mutate(
    `rs77915569` = factor(
      `rs77915569`,
      levels = c("TT","TC","CC")
    )
  )

data_plot$IdInfect = as.character(data_plot$IdInfect)
chr3_CA = left_join(data_plot,
                    chr3_CA,
                    "IdInfect")

chr3_CA |> 
  tidyplot(x=rs77915569, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" ,"purple", "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)


chr17_AFMID = chr17_AFMID |> 
  filter(!is.na(`rs3642`))

ggpubr::ggboxplot(
  chr17_AFMID,
  title = "tryptophan",
  x = "rs3642",
  y = "HMDB0000929",
  ylab = "MFI"
) +
  ggpubr::stat_compare_means(  # Pairwise comparisons
    method = "wilcox.test",
    label = "p.format"
  )

chr17_AFMID_v2 = chr17_AFMID |> 
  mutate(
    `rs3642` = factor(
      `rs3642`,
      levels = c("GG","GA","AA"),
      labels = c("GG","TC+CC","TC+CC")
    )
  )

ggpubr::ggboxplot(
  chr22_PRODH_v2,
  title = "Proline",
  x = "rs13056032",
  y = "HMDB0000162",
  ylab = "MFI"
) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )



# tryptophan v2 -----------------------------------------------------------

tryptophan = data_plot |> 
  select(
    IdInfect,
    status,
    HMDB0000929
  )
tryptophan$IdInfect = as.character(tryptophan$IdInfect)
str(tryptophan)

chr17_AFMID = chr17[rownames(chr17) == "17:76162178:G:A", ] 
rownames(chr17_AFMID) = "rs3642"
str(chr17_AFMID)

chr17_AFMID = as.data.frame(t(chr17_AFMID)) |> 
  rownames_to_column("IdInfect") |> 
  mutate(IdInfect = gsub("X","",IdInfect)) |> 
  mutate(
    `rs3642` = if_else(
      `rs3642` < 0.5,
      "GG", if_else(
        `rs3642` > 1.5,
        "AA",
        "GA"
      )
    )
  ) |> 
  mutate(
    `rs3642` = factor(
      `rs3642`,
      levels = c("GG","GA","AA")
    )
  )

str(tryptophan)

chr17_AFMID = left_join(tryptophan,
                        chr17_AFMID,
                        "IdInfect")

chr17_AFMID = chr17_AFMID |> 
  filter(!is.na(`rs3642`))

ggpubr::ggboxplot(
  chr17_AFMID,
  title = "Tryptophan",
  x = "rs3642",
  width = 0.6,
  y = "HMDB0000929",
  ylab = "Log2 (peak ion intensity)",
  add = c("jitter"),
  add.params = list(alpha = 0.3)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c(1,2)),  # Pairwise comparisons
    method = "wilcox.test",
    label = "p.format"
  )

chr17_AFMID |> 
  tidyplot(x=rs3642, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" , "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

chr17_AFMID$rs3642 = factor(chr17_AFMID$rs3642,
  levels = c("GG","GA")
)

library(gtsummary)
gtsummary::tbl_summary(chr17_AFMID |> select(
  status,
  rs3642
), by="rs3642") |> 
  gtsummary::add_p()
?add_p
ggsave("Tryptophan proportion_dom.png")


chr17_AFMID_v2 = chr17_AFMID |> 
  mutate(
    `rs2018564` = factor(
      `rs2018564`,
      levels = c("CC","CT","TT"),
      labels = c("cc","TC+CC","TC+CC")
    )
  )

ggpubr::ggboxplot(
  chr22_PRODH_v2,
  title = "Proline",
  x = "rs13056032",
  y = "HMDB0000162",
  ylab = "MFI"
) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )

write.csv(chr17_AFMID,"chr17_AFMID.csv")

# violin plot metabolites INFECT HHC -------------------------------------------------

data2 = data_plot |> 
  filter(strict_0.15_conversion != "Uncertain")

# Create individual plots
proline_plot <- create_metabolite_plot(data2, "HMDB0000162", "Proline")
ltb4_plot <- create_metabolite_plot(data2, "HMDB0001085", "Leukotriene B4")
trp_plot <- create_metabolite_plot(data2, "HMDB0000929", "Tryptophan")
gly_plot <- create_metabolite_plot(data2, "HMDB0000123", "Glycine")
hode_plot <- create_metabolite_plot(data2, "HMDB0004667", "13-HODE")
pgf2_plot <- create_metabolite_plot(data2, "HMDB0001139", "Prostaglandin F2α")
indoleactici_plot <- create_metabolite_plot(data2, "HMDB0000197", "Indoleacetic acid")
# Combine plots in one row
combined_plot <- gly_plot+ proline_plot  + trp_plot   +indoleactici_plot+
  ltb4_plot +hode_plot + pgf2_plot +  
  plot_layout(nrow = 4)

export = 
  data2 |> 
  select(IdInfect,
         strict_0.15_conversion,
         HMDB0000162,
         HMDB0001085,
         HMDB0000929,
         HMDB0000123,
         HMDB0004667,
         HMDB0001139,
         HMDB0000197
         )

write.csv(export,"INFECT metabolomics_fig5D.csv")

# Display combined plot
combined_plot

ggsave("Fig 5D_combined violin.png",
       width = 20,
       height = 40,
       units = "cm")

?ggsave
