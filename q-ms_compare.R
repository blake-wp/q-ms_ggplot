rm(list = ls())


library(tidyverse)
quant <- read_csv("Quantitative_comparisons_2021.csv")
glimpse(quant)
names(quant) <- make.names(names(quant))
intens <- read_csv("intensity_comparisons_2021.csv")
glimpse(intens)
names(intens) <- make.names(names(intens))

library(janitor)
compare_df_cols(quant, intens)

# Merge on Comparison, Protein.Description, Protein/Protein.Name
intens <- rename(intens, Protein = Protein.Name)
df <- merge(quant, intens, by = c("Comparison", "Protein.Description", "Protein"))
glimpse(df)
rm(intens, quant)

# Clean up descriptions
to_remove <- c(" \\[Apis mellifera\\]", " precursor", " preproprotein", " \\[FAD, quinone\\]")
df <- 
  df %>% 
  mutate(Protein.Description = str_remove_all(Protein.Description, paste(to_remove, collapse = "|"))) %>% 
  mutate(Protein.Description = str_replace_all(Protein.Description, "uncharacterized protein", "UP"))

# Calculate coordinates for significance symbols, define labels
lab_pos_adjust <- 1.5
df$lab_pos <- ifelse(df$Log2.Fold.Change > 0, df$Log2.Max + lab_pos_adjust, df$Log2.Min - lab_pos_adjust)
df$sig_pos <- 
  case_when(
    df$lab_pos < 0 ~ case_when(
      df$Adjusted.P.Value < 0.001 ~ "***",
      df$Adjusted.P.Value < 0.01 ~ " **",
      df$Adjusted.P.Value < 0.05 ~ "  *",
      .default = ""
    ),
    .default = case_when(
      df$Adjusted.P.Value < 0.001 ~ "***",
      df$Adjusted.P.Value < 0.01 ~ "** ",
      df$Adjusted.P.Value < 0.05 ~ "*  ",
      .default = ""
    )
  )

comparison_names <- c('M_vs_R' = "Manuka\nvs\nRoyal jelly",
                      'C_vs_R' = "Clover\nvs\nRoyal jelly",
                      'MRJPs' = "MRJPs",
                      'Glycosidases' = "Glycosidases",
                      'Oxidoreductases' = "Oxidoreductases",
                      'Esterases' = "Esterases",
                      'Peptidases' = "Peptidases",
                      'Peptidase inhibitors' = "Peptidase inhibitors",
                      'Immune response' = "Immune response",
                      'Lipid transport' = "Lipid transport",
                      'Allergen' = "Allergen",
                      'Uncharacterized' = "Uncharacterized")
df$Comparison <- factor(df$Comparison, levels = c("M_vs_R", "C_vs_R"))
df$function. <- factor(df$function., levels = c("MRJPs", "Glycosidases", "Oxidoreductases", "Esterases", "Peptidases", "Peptidase inhibitors", "Immune response", "Lipid transport", "Allergen", "Uncharacterized"))
glimpse(df)

# Specifically edit the value for 'M_vs_R' 'lysozyme' 'Mean' as it is 100x less than any other and stretches the legend down to below 0.001%
df <- df %>% mutate(Mean = ifelse(Comparison == "M_vs_R" & Protein.Description == 'lysozyme', 0.00041, Mean))

# Plot
fill_colour <- "green"

# Titles and caption.
title <- "Comparative analyses of mānuka and clover honey to royal jelly by SWATH-MS"
subtitle <- "With accompanying R code for rendering various aspects of the graph"
caption <- 
  "Figure. Quantitative comparisons of proteins identified by both analysis methods, expressed as the log2-fold change in sum 
  of peak intensities for each honey protein compared to the corresponding protein found in royal jelly. Apisimin was added
  from the SWATH-MS analysis, giving 32 proteins in total. The shading of the bars represents the mean peptide intensity of
  protein in the honey samples (mānuka, n = 12; clover, n = 6, royal jelly n = 8). Significant differences areindicated by
  * P < 0.05, ** P < 0.01, *** P < 0.001. https://doi.org/10.1371/journal.pone.0272898.g002"

# Plot
ggplot(data = df,
       mapping = aes(x = Log2.Fold.Change, y = reorder(Protein.Description, desc(group_order)))) +
  geom_bar(data = df,
           aes(fill = Mean),
           stat = "identity",
           color = "black",
           width = 0.8) +
  geom_errorbar(width = 0.4,
                aes(xmin = Log2.Min, xmax = Log2.Max)) +
  scale_fill_gradient(guide = "legend",
                      low = "white",
                      high = fill_colour,
                      trans = "log10",
                      breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.30),
                      labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")) +
  scale_x_continuous(limits = c(-9.5, 9.5),
                     breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
  labs(title = title,
       subtitle = subtitle,
       caption = caption, 
       x = expression(paste(Log[2], " fold-change in honey protein")),
       y = NULL) +
  theme(panel.border = element_rect(color = "grey60", fill = NA),
        axis.line.x = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey90"),
        panel.grid.major.x = element_line(linewidth = 0.25),
        legend.position = "right",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot",
        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.text.y.right = element_text(angle = 0),
        strip.background = element_rect(color = "grey60", fill = "white")) +
  guides(fill = guide_colorbar("Protein\nabundance\nin honey",
                               barwidth = 0.5, barheight = 10,
                               ticks.colour = "black",
                               frame.colour = "black",
                               frame.linewidth = 0.8)) +
  geom_hline(yintercept = seq(1.5, length(unique(df$Protein.Description))-0.5, 1),
             lwd = 0.2,
             color = "grey90") +
  geom_vline(xintercept = 0) +
  geom_text(aes(x = lab_pos, label = sig_pos),
            position = position_dodge(0.9),
            vjust = 0.7,
            size = 3,
            color = "grey30") +
  facet_grid(function. ~ Comparison,
             scales = "free_y",
             space = "free_y",
             labeller = as_labeller(comparison_names))

