################################################################################
# Differential Abundance Testing with ALDEx2

{
  library(tidyverse)
  library(phyloseq)
  library(ALDEx2)
}

load("Jhae.microbiome.RData", verbose = TRUE)

# Extract and format the count matrix
counts_july <- round(as(otu_table(ps_july), "matrix"))

# Extract the comparison vector
conditions_july <- as.character(sample_data(ps_july)$host_plant)

set.seed(123) 

# mc.samples = 128 is the recommended default for Monte Carlo sampling
# denom = "all" uses the whole community as the reference frame
aldex_host <- aldex(
  reads = counts_july, 
  conditions = conditions_july, 
  mc.samples = 128, 
  test = "t", 
  effect = TRUE, 
  denom = "all"
)

# Merge results with taxonomy and filter for significance
taxa_july <- as.data.frame(as(tax_table(ps_july), "matrix"))

aldex_host_sig <- aldex_host |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_july |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  # Filter by Wilcoxon Benjamini-Hochberg corrected p-value < 0.05
  filter(wi.eBH < 0.05) |>
  # Sort by the largest absolute effect sizes
  arrange(desc(abs(effect))) |> 
  dplyr::select(ASV_ID, diff.btw, effect, wi.eBH, Family, Genus)

aldex_host_sig  

# ALDEx2 result dataframe with taxonomy
aldex_host_plot_df <- aldex_host |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_july |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  mutate(
    # Categorize by significance AND strong effect size (> 1 or < -1)
    Significance = case_when(
      wi.eBH < 0.05 & effect > 1 ~ "Enriched in Koelreuteria",
      wi.eBH < 0.05 & effect < -1 ~ "Enriched in Cardiospermum",
      TRUE ~ "Not Significant"
    ),
    # Create a label for only the extreme drivers (effect > |2|) to prevent plot clutter
    Taxa_Label = case_when(
      wi.eBH < 0.05  & !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_aldex_host <- ggplot(aldex_host_plot_df, aes(x = effect, y = -log10(wi.ep), color = Significance)) +
  geom_point(alpha = 0.85, size = 2.5) +
  geom_text_repel(
    aes(label = Taxa_Label), size = 2.5, 
    box.padding = 0.5, point.padding = 0.3, 
    force = 4, max.overlaps = Inf, 
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Enriched in Cardiospermum" = color.by.plant[1], 
    "Enriched in Koelreuteria" = color.by.plant[2], 
    "Not Significant" = "grey50"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  expand_limits(x = c(-6, 4), y = c(0, 9)) +
  labs(
    x = "Effect Size",
    y = expression(-log[10]("p-value"))
  )

plot_aldex_host

##########
# By month
counts_grt <- round(as(otu_table(ps_grt), "matrix"))
conditions_grt <- as.character(sample_data(ps_grt)$month)

set.seed(123) 
aldex_month <- aldex(
  reads = counts_grt, 
  conditions = conditions_grt, 
  mc.samples = 128, 
  test = "t", 
  effect = TRUE, 
  denom = "all"
)

# Merge results with taxonomy and filter for significance
taxa_grt <- as.data.frame(as(tax_table(ps_grt), "matrix"))

aldex_month_sig <- aldex_month |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_grt |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  # Filter by Wilcoxon Benjamini-Hochberg corrected p-value < 0.05
  filter(wi.eBH < 0.05) |>
  # Sort by the largest absolute effect sizes
  arrange(desc(abs(effect))) |> 
  dplyr::select(ASV_ID, diff.btw, effect, wi.eBH, Family, Genus)

aldex_month_sig  

aldex_month_plot_df <- aldex_month |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_grt |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  mutate(
    # Categorize by significance AND strong effect size (> 1 or < -1)
    Significance = case_when(
      wi.eBH < 0.05 & effect > 0 ~ "Enriched in July",
      wi.eBH < 0.05 & effect < 0 ~ "Enriched in January",
      TRUE ~ "Not Significant"
    ),
    # Create a label for only the extreme drivers (effect > |1|) to prevent plot clutter
    Taxa_Label = case_when(
      wi.eBH < 0.05  & !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  & !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      wi.eBH < 0.05  ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_aldex_month <- ggplot(aldex_month_plot_df, aes(x = effect, y = -log10(wi.ep), color = Significance)) +
  geom_point(alpha = 0.85, size = 2.5) +
  geom_text_repel(
    aes(label = Taxa_Label), size = 2.5, 
    box.padding = 0.5, point.padding = 0.3, 
    force = 4, max.overlaps = Inf, 
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Enriched in July" = color.by.month[2], 
    "Enriched in January" = color.by.month[1], 
    "Not Significant" = "grey50"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  expand_limits(x = c(-6, 4), y = c(0, 9)) +
  labs(
    x = "Effect Size",
    y = expression(-log[10]("p-value"))
  )

plot_aldex_month


aldex_figure <- (
  plot_aldex_host | plot_aldex_month
) +
  plot_annotation(tag_levels = 'A') 
aldex_figure

ggsave(
  filename = "plots/FigX.aldex.volcanos.png", 
  plot = aldex_figure,
  width = 6.5, height = 3.5, scale = 1.65, 
  dpi = 600, bg = "white"
)


################################################################################
# Save the key objects
save(
  counts_july, conditions_july, aldex_host, taxa_july, aldex_host_sig, aldex_host_plot_df, plot_aldex_host,
  counts_grt, conditions_grt, aldex_month, taxa_grt, aldex_month_sig, aldex_month_plot_df, plot_aldex_month,
  aldex_figure,
  file = "Jhae.microbiome.aldex.RData"
)
