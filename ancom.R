################################################################################
# Differential Abundance Testing with ANCOM-BC2

{
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
}

load("Jhae.microbiome.RData", verbose = TRUE)

# Strip out any remaining unclassified Phyla or Families before testing
ps_phylum <- subset_taxa(ps_phylum, !is.na(Phylum))
ps_family <- subset_taxa(ps_family, !is.na(Family))

#####
# ANCOM-BC2 for Phylum level
ancom_phylum <- ancombc2(
  data = ps_phylum,
  fix_formula = "month + host_plant",  # Tests both significant drivers
  p_adj_method = "holm",               # Adjusts p-values to prevent false positives
  group = "month",                     # Specifies the primary experimental group
  struc_zero = TRUE,                   # Detects if a Phylum is 100% missing in a specific plant
  neg_lb = TRUE,                       # Helps classify true differences
  alpha = 0.05
)

# Join the taxonomy names onto the ANCOM-BC2 results
taxa <- as.data.frame(tax_table(ps_phylum)) %>%
  rownames_to_column(var = "taxon") 

phylum_results <- ancom_phylum$res %>%
  left_join(taxa, by = "taxon") %>%
  # Filter for significant results
  filter(if_any(starts_with("diff"), ~ . == TRUE)) %>%
  relocate(Phylum, .after = taxon) %>%
  # Structural Zeros?
  left_join(ancom_phylum$zero_ind, by = "taxon")

phylum_results

#####
# ANCOM-BC2 for family level
ancom_family <- ancombc2(
  data = ps_family,
  fix_formula = "month + host_plant",  # Tests both significant drivers
  p_adj_method = "holm",               # Adjusts p-values to prevent false positives
  group = "month",                     # Specifies the primary experimental group
  struc_zero = TRUE,                   # Detects if a family is 100% missing in a specific plant
  neg_lb = TRUE,                       # Helps classify true differences
  alpha = 0.05
)

taxa <- as.data.frame(tax_table(ps_family)) %>%
  rownames_to_column(var = "taxon") 

family_results <- ancom_family$res %>%
  left_join(taxa, by = "taxon") %>%
  # Filter for significant results
  filter(if_any(starts_with("diff"), ~ . == TRUE)) %>%
  relocate(Family, .after = taxon) %>%
  # Structural Zeros?
  left_join(ancom_family$zero_ind, by = "taxon")

family_results

#####
# ANCOM-BC2 for the ASV level
ancom_asv <- ancombc2(
  data = ps_asv,
  fix_formula = "month + host_plant",  # Tests both significant drivers
  p_adj_method = "holm",               # Adjusts p-values to prevent false positives
  group = "host_plant",                # Specifies the primary experimental group
  struc_zero = TRUE,                   # Detects if a asv is 100% missing in a specific plant
  neg_lb = TRUE,                       # Helps classify true differences
  alpha = 0.05
)

# Join the taxonomy names onto the ANCOM-BC2 results
taxa <- as.data.frame(tax_table(ps_asv)) %>%
  rownames_to_column(var = "taxon") 

asv_results <- ancom_asv$res %>%
  left_join(taxa, by = "taxon") %>%
  # Filter for significant results
  filter(if_any(starts_with("diff"), ~ . == TRUE)) %>%
  relocate(Genus, Species, .after = taxon) %>%
  # Structural Zeros?
  left_join(ancom_asv$zero_ind, by = "taxon")

asv_results

#####
# Visualization for the ANCOM results

plot_ancom_phylum_by_month <- phylum_results %>%
  filter(diff_monthJuly == TRUE) %>%
  arrange(desc(lfc_monthJuly)) %>%
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) %>%
  mutate(Shift = ifelse(lfc_monthJuly > 0, "Enriched in July", "Enriched in January")) %>%
  ggplot(aes(x = lfc_monthJuly, y = Phylum, fill = Shift)) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(xmin = lfc_monthJuly - se_monthJuly,xmax = lfc_monthJuly + se_monthJuly),
    width = 0.2, linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, color = "gray30") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits = c(-6,5), breaks = c(-6:5)) +
  # scale_fill_manual(values = c("#4575b4","#FF8559")) +
  scale_fill_manual(values = color.by.month[2]) +
  labs(
    title = "Phylum-level seasonal shifts",
    x = "LFC (January to July)", y = NULL, 
    fill = "Seasonal Shift"
  ) 

plot_ancom_phylum_by_month

plot_ancom_family_by_host_plant <- family_results %>%
  filter(diff_host_plantGRT == TRUE) %>% 
  arrange(desc(lfc_host_plantGRT)) %>%
  mutate(Family = factor(Family, levels = unique(Family))) %>%
  mutate(Shift = ifelse(lfc_host_plantGRT > 0, "Enriched in GRT", "Enriched in BV")) %>% 
  ggplot(aes(x = lfc_host_plantGRT, y = Family, fill = Shift)) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(xmin = lfc_host_plantGRT - se_host_plantGRT,xmax = lfc_host_plantGRT + se_host_plantGRT),
    width = 0.2, linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, color = "gray30") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits = c(-6,5), breaks = c(-6:5)) +
  scale_fill_manual(values = color.by.plant) +
  labs(
    title = "Family-level host shifts",
    x = "LFC (BV to GRT)", y = NULL, 
    fill = "Host Plant Bias"
  ) 

plot_ancom_family_by_host_plant

plot_ancom_family_by_month <- family_results %>%
  filter(diff_monthJuly == TRUE) %>% 
  arrange(desc(lfc_monthJuly)) %>%
  mutate(Family = factor(Family, levels = unique(Family))) %>%
  mutate(Shift = ifelse(lfc_monthJuly > 0, "Enriched in July", "Enriched in January")) %>% 
  ggplot(aes(x = lfc_monthJuly, y = Family, fill = Shift)) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(xmin = lfc_monthJuly - se_monthJuly,xmax = lfc_monthJuly + se_monthJuly),
    width = 0.2, linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, color = "gray30") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits = c(-6,5), breaks = c(-6:5)) +
  scale_fill_manual(values = color.by.month) +
  labs(
    title = "Family-level seasonal shifts",
    x = "LFC (January to July)", y = NULL, 
    fill = "Seasonal Shift"
  ) 

plot_ancom_family_by_month

plot_ancom_asv_by_host_plant <- asv_results %>%
  filter(diff_host_plantGRT == TRUE) %>% 
  arrange(desc(lfc_host_plantGRT)) %>%
  mutate(
    Granular_Taxon = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom),
    Granular_Taxon = str_replace(
      Granular_Taxon, 
      "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 
      "Rhizobium complex"
    ),
    Plot_Name = paste0(Granular_Taxon, "\n(", substr(taxon, 1, 6), "...)")
  ) %>% 
  mutate(Plot_Name = factor(Plot_Name, levels = unique(Plot_Name))) %>%
  mutate(Shift = ifelse(lfc_host_plantGRT > 0, "Enriched in GRT", "Enriched in BV")) %>% 
  ggplot(aes(x = lfc_host_plantGRT, y = Plot_Name, fill = Shift)) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(xmin = lfc_host_plantGRT - se_host_plantGRT,xmax = lfc_host_plantGRT + se_host_plantGRT),
    width = 0.2, linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, color = "gray30") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits = c(-6,5), breaks = c(-6:5)) +
  scale_fill_manual(values = color.by.plant) +
  labs(
    title = "ASV-level host shifts",
    x = "LFC (BV to GRT)", y = NULL, 
    fill = "Host Plant Bias"
  ) 

plot_ancom_asv_by_host_plant

plot_ancom_asv_by_month <- asv_results %>%
  filter(diff_monthJuly == TRUE) %>% 
  arrange(desc(lfc_monthJuly)) %>%
  mutate(
    Granular_Taxon = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom),
    Granular_Taxon = str_replace(
      Granular_Taxon, 
      "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 
      "Rhizobium complex"
    ),
    Plot_Name = paste0(Granular_Taxon, "\n(", substr(taxon, 1, 6), "...)")
  ) %>% 
  mutate(Plot_Name = factor(Plot_Name, levels = unique(Plot_Name))) %>%
  mutate(Shift = ifelse(lfc_monthJuly > 0, "Enriched in July", "Enriched in January")) %>% 
  ggplot(aes(x = lfc_monthJuly, y = Plot_Name, fill = Shift)) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(xmin = lfc_monthJuly - se_monthJuly, xmax = lfc_monthJuly + se_monthJuly),
    width = 0.2, linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, color = "gray30") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits = c(-6,5), breaks = c(-6:5)) +
  scale_fill_manual(values = color.by.month) +
  labs(
    title = "ASV-level seasonal shifts",
    x = "LFC (January to July)", y = NULL, 
    fill = "Seasonal Shift"
  ) 

plot_ancom_asv_by_month

ancom_figure <- (
  (plot_ancom_phylum_by_month | plot_spacer()) /
  (plot_ancom_family_by_month | plot_ancom_family_by_host_plant) /
    (plot_ancom_asv_by_month | plot_ancom_asv_by_host_plant)
) +
# plot_annotation(tag_levels = 'A') +
plot_layout(heights = c(1, 9, 3), widths = c(1, 1), guides = "collect")

ancom_figure

ggsave(
  filename = "plots/ANCOM.png", 
  plot = ancom_figure,
  width = 9, height = 5.6, scale = 1, 
  dpi = 600, bg = "white"
)

################################################################################
# Save the key objects
save(
  ps_phylum, ps_family, 
  ancom_phylum, phylum_results,
  ancom_family, family_results,
  ancom_asv, asv_results,
  plot_ancom_phylum_by_month,
  plot_ancom_family_by_month, plot_ancom_family_by_host_plant,
  plot_ancom_asv_by_month, plot_ancom_asv_by_host_plant,
  ancom_figure,
  file = "Jhae.microbiome.ancom.RData"
)
