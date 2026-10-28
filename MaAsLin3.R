################################################################################
# MaAsLin3  (Microbiome Multivariable Associations with Linear Models)

{
  library(phyloseq)
  library(tidyverse)
  library(maaslin3)
}

load("Jhae.microbiome.RData", verbose = TRUE)


taxa_info <- as.data.frame(as(tax_table(ps_main), "matrix")) |>
  rownames_to_column(var = "feature")

# Extract count matrix
counts_main <- as(otu_table(ps_main), "matrix")
counts_main <- as.data.frame(counts_main)

# Full Model
maaslin_full <- maaslin3(
  input_data = counts_main,
  input_metadata = metadata_main,
  output = "maaslin.full.model",
  formula = "~ host_plant + month + sex + green_space",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)

maaslin_full_sig <- maaslin_full$fit_data_prevalence$results %>%
  left_join(taxa_info, by = "feature") |>
  filter(qval_individual < 0.05) %>%
  arrange(qval_individual)

maaslin_full_sig

################################################################################
# Host plant model

counts_july <- as(otu_table(ps_july), "matrix")
counts_july <- as.data.frame(counts_july)

maaslin_host <- maaslin3(
  input_data = counts_july,
  input_metadata = metadata_july,
  output = "maaslin.host.model",
  formula = "~ host_plant",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)

maaslin_host_sig <- maaslin_host$fit_data_prevalence$results %>%
  left_join(taxa_info, by = "feature") |>
  filter(qval_individual < 0.05) %>%
  arrange(qval_individual)

maaslin_host_sig

maaslin_host_plot_df <- maaslin_host$fit_data_prevalence$results %>%
  left_join(taxa_info, by = "feature") |>
  mutate(
    Significance = case_when(
      qval_individual < 0.05 & coef > 1 ~ "Enriched in K. paniculata",
      qval_individual < 0.05 & coef < -1 ~ "Enriched in C. corindum",
      TRUE ~ "Not Significant"
    ),
    Taxa_Label = case_when(
      qval_individual < 0.05  & !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(feature, 1, uuid.abbrev), ")"),
      qval_individual < 0.05  & !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(feature, 1, uuid.abbrev), ")"),
      qval_individual < 0.05  & !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(feature, 1, uuid.abbrev), ")"),
      qval_individual < 0.05  & !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(feature, 1, uuid.abbrev), ")"),
      qval_individual < 0.05  & !is.na(Phylum) & Phylum != "NA" ~ paste0(Phylum, "\n(", substr(feature, 1, uuid.abbrev), ")"),
      qval_individual < 0.05  ~ paste0("Unclassified\n(", substr(feature, 1, uuid.abbrev), ")")
    )
  ) 

plot_maaslin_host <- ggplot(maaslin_host_plot_df, aes(x = coef, y = -log10(pval_individual), color = Significance)) +
  geom_point(alpha = 0.85, size = 2.5) +
  geom_text_repel(
    aes(label = Taxa_Label), size = 2.5,
    box.padding = 0.5, point.padding = 0.3,
    force = 4, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Enriched in C. corindum" = color.by.plant[1],
    "Enriched in K. paniculata" = color.by.plant[2],
    "Not Significant" = "grey50"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  # expand_limits(x = c(-6, 4), y = c(0, 6.25)) +
  labs(
    x = "Effect Size",
    y = expression(-log[10]("p-value"))
  )

plot_maaslin_host



################################################################################
# Host plant model with random effect for collection site

counts_july <- as(otu_table(ps_july), "matrix")
counts_july <- as.data.frame(counts_july)

maaslin_host_mixed <- maaslin3(
  input_data = counts_july,
  input_metadata = metadata_july,
  output = "maaslin.host.model",
  formula = "~ host_plant + (1 | site_name)",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)

maaslin_host_mixed_sig <- maaslin_host$fit_data_prevalence$results %>%
  left_join(taxa_info, by = "feature") |>
  filter(qval_individual < 0.05) %>%
  arrange(qval_individual)

maaslin_host_mixed_sig
# Same taxa as the fixed effect model!


################################################################################
# Month model

counts_grt <- as(otu_table(ps_grt), "matrix")
counts_grt <- as.data.frame(counts_grt)

maaslin_month <- maaslin3(
  input_data = counts_grt,
  input_metadata = metadata_grt,
  output = "maaslin.month.model",
  formula = "~ month",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)

maaslin_month_sig <- maaslin_month$fit_data_prevalence$results %>%
  left_join(taxa_info, by = "feature") |>
  filter(qval_individual < 0.05) %>%
  arrange(qval_individual)

maaslin_month_sig

################################################################################
# Green space / Urbanicity

maaslin_green_space <- maaslin3(
  input_data = counts_main,
  input_metadata = metadata_main,
  output = "maaslin.green_space.model",
  formula = "~ green_space",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)
# All associations had errors or were insignificant.

maaslin_urbanicity <- maaslin3(
  input_data = counts_main,
  input_metadata = metadata_main,
  output = "maaslin.urbanicity.model",
  formula = "~ urbanicity",
  normalization = "TSS",    # Total Sum Scaling (default)
  transform = "LOG",        # Log transformation (default)
  max_significance = 0.05   # FDR threshold (q-value)
)
# All associations had errors or were insignificant.


################################################################################
# Combined results

maaslin_results <- bind_rows(
  maaslin_host_sig,
  maaslin_month_sig
) |>
  dplyr::select(-name, -null_hypothesis, -Kingdom)

write.csv(maaslin_results, "tables/TableS9.MaAsLin3.csv", row.names = FALSE)

################################################################################
# Save the key objects
save(
  maaslin_host_sig, maaslin_month_sig, maaslin_results,
  file = "Jhae.microbiome.MaAsLin3.RData"
)

