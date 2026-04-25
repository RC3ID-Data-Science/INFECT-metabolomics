# Load data and config
source("../00_load_data.R")
create_output_dirs()

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

data = data |\u003e 
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


# Data preparation with two-line label
data_plot <- merge %\u003e%
#   mutate(status = factor(status,
#                          levels = c("Persistently_uninfected", "Converter"),
#                          labels = c("Persistently\nIGRA-negatives", "IGRA-converters")  # Added line break
#  ))
  mutate(status = factor(status,
                         levels = c("Converter", "Persistently_uninfected"),
                         labels = c("IGRA-converters", "Persistently\nIGRA-negatives")  # Added line break
  )) #%\u003e%
  #filter(strict_0.15_conversion != "Uncertain")
data_plot$sex_contact = factor(data_plot$sex_contact,
                               levels = c("Female","Male"))


# Load genotype data
chr22 = read_rds(CHR22_DOSAGE)
chr17 = read_rds(CHR17_DOSAGE)

proline = data_plot |> 
  select(
    IdInfect,
    status,
    HMDB0000162
  )

proline$IdInfect = as.character(proline$IdInfect)
str(proline)


chr22_PRODH = chr22[rownames(chr22) == "22:18938193:T:C", ] 
rownames(chr22_PRODH) = "rs13056032"
str(chr22_PRODH)

chr22_PRODH = as.data.frame(t(chr22_PRODH)) |> 
  rownames_to_column("IdInfect") |> 
  mutate(IdInfect = gsub("X","",IdInfect)) |> 
  mutate(
    `rs13056032` = if_else(
      `rs13056032` < 0.5,
      "TT", if_else(
        `rs13056032` > 1.5,
        "CC",
        "TC"
      )
    )
    ) |> 
  mutate(
    `rs13056032` = factor(
      `rs13056032`,
      levels = c("TT","TC","CC")
    )
  )

str(proline)

chr22_PRODH = left_join(proline,
                        chr22_PRODH,
                        "IdInfect")



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



chr22_PRODH |> 
  tidyplot(x=rs13056032, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" , "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

ggsave(file.path(FIGURE_DIR, "PRODH proportion.png"))

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

ggsave(file.path(FIGURE_DIR, "PRODH proportion_dom.png"))

write.csv(chr22_PRODH, file.path(OUTPUT_DIR, "chr22_PRODH.csv"))


# proline second rsid -----------------------------------------------------



chr22_PRODH_v3 = chr22[rownames(chr22) == "22:18938190:G:A", ] 
rownames(chr22_PRODH_v3) = "rs13054729"
str(chr22_PRODH_v3)


chr22_PRODH_v3 = as.data.frame(t(chr22_PRODH_v3)) |> 
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

chr22_PRODH_v3 = left_join(proline,
                           chr22_PRODH_v3,
                           "IdInfect")
chr22_PRODH_v3 = chr22_PRODH_v3 |> 
  dplyr::filter(!is.na(`rs13054729`))

ggpubr::ggboxplot(
  chr22_PRODH_v3,
  title = "Proline",
  x = "rs13054729",
  width = 0.6,
  y = "HMDB0000162",
  ylab = "Log2 (peak ion intensity)",
  add = c("jitter"),
  add.params = list(alpha = 0.3)
) +
  ggpubr::stat_compare_means(
    comparisons = list(c(1,2), c(1,3), c(2,3)), 
    method = "wilcox.test",
    label = "p.format"
  )

library(tidyplots)

chr22_PRODH_v3 |> 
  tidyplot(x=rs13054729, color=status) |>
  add_barstack_relative() |> 
  adjust_colors(new_colors = c("#FF2323" , "#2323FF")) |> 
  adjust_legend_title("") |> 
  adjust_legend_position("top") |> 
  adjust_y_axis_title("Proportion") |> 
  adjust_size(height = 100, width = 60) %>% 
  adjust_font(family = "sans", fontsize = 12)

ggsave(file.path(FIGURE_DIR, "PRODH proportion_dom.png"))


# tryptophan --------------------------------------------------------------

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

ggsave(file.path(FIGURE_DIR, "Tryptophan proportion_dom.png"))


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

write.csv(chr17_AFMID, file.path(OUTPUT_DIR, "chr17_AFMID.csv"))

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

write.csv(export, file.path(OUTPUT_DIR, "INFECT_metabolomics_fig5D.csv"))

# Display combined plot
combined_plot

ggsave(file.path(FIGURE_DIR, "Fig 5D_combined violin.png"),
       width = 20,
       height = 40,
       units = "cm")

?ggsave



