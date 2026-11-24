################################################################################
# Differential Abundance Testing with ANCOM-BC2

setwd("~/Documents/3.research/Jhae.microbiome")

# Clear memory
# {
#   rm(list=ls())
#   .rs.restartR()
#   gc()
# }

{
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
}

load("Jhae.microbiome.RData", verbose = TRUE)

# Collapse the data to different  taxonomic levels for each subset
{
  # Convert area to hectares 
  sample_data(ps_main)$intensive_development_ha <- sample_data(ps_main)$intensive_development / 10000
  
  ps_main_phylum <- tax_glom(ps_main, taxrank = "Phylum")
  ps_main_phylum <- subset_taxa(ps_main_phylum, !is.na(Phylum))
  ps_main_family <- tax_glom(ps_main, taxrank = "Family")
  ps_main_family <- subset_taxa(ps_main_family, !is.na(Family))
  ps_main_asv <- ps_main

  ps_july_phylum <- tax_glom(ps_july, taxrank = "Phylum")
  ps_july_phylum <- subset_taxa(ps_july_phylum, !is.na(Phylum))
  ps_july_family <- tax_glom(ps_july, taxrank = "Family")
  ps_july_family <- subset_taxa(ps_july_family, !is.na(Family))
  ps_july_asv <- ps_july

  ps_grt_phylum <- tax_glom(ps_grt, taxrank = "Phylum")
  ps_grt_phylum <- subset_taxa(ps_grt_phylum, !is.na(Phylum))
  ps_grt_family <- tax_glom(ps_grt, taxrank = "Family")
  ps_grt_family <- subset_taxa(ps_grt_family, !is.na(Family))
  ps_grt_asv <- ps_grt
}

# Function to clean and standardize output tables
clean_ancom_table <- function(df, level_name, feature_col, contrast_name) {

  # Define dynamic column names
  lfc_col  <- paste0("lfc_", contrast_name)
  se_col   <- paste0("se_", contrast_name)
  w_col    <- paste0("W_", contrast_name)
  p_col    <- paste0("p_", contrast_name)
  q_col    <- paste0("q_", contrast_name)
  diff_col <- paste0("diff_", contrast_name)

  df %>%
    mutate(
      Analysis_Level = level_name,
      Feature_ID = !!sym(feature_col),
      Contrast = contrast_name
    ) %>%
    dplyr::select(
      Analysis_Level,
      Feature_ID,
      Contrast,
      # Taxa columns (if they exist)
      any_of(c("Phylum", "Class", "Order", "Family", "Genus", "Species")),

      # Core statistics mapped dynamically
      LogFoldChange = !!sym(lfc_col),
      StdError = !!sym(se_col),
      W_stat = !!sym(w_col),
      p_val = !!sym(p_col),
      q_val_FDR = !!sym(q_col),
      Significant = !!sym(diff_col),

      # Structural zero columns if they exist.
      starts_with("structural_zero")
    ) |>
    relocate(any_of(c("Phylum", "Class", "Order", "Family", "Genus", "Species")), .after = last_col())
}

# Explanation of `ancombc2` arguments
#   p_adj_method = "holm", # Adjusts p-values to prevent false positives
#   struc_zero = TRUE,     # Detects if a taxon is missing in one group
#   neg_lb = TRUE,         # Call structural zero using asymptotic lower bound

# Host plant
{
  # Host Plant / Phylum
  ancom_phylum_host <- ancombc2(
    data = ps_july_phylum,
    fix_formula = "host_plant", group = "host_plant",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  # Join the taxonomy names onto the ANCOM-BC2 results
  taxa <- as.data.frame(tax_table(ps_july_phylum)) |>
    rownames_to_column(var = "taxon")

  ancom_phylum_host_results <- ancom_phylum_host$res |>
    left_join(taxa, by = "taxon") |>
    # Filter for significant results
    filter(q_host_plantGRT <= 0.05) |>
    relocate(Phylum, .after = taxon) |>
    # Structural Zeros?
    left_join(ancom_phylum_host$zero_ind, by = "taxon") |>
    dplyr::select(-taxon, -Kingdom)

  ##########
  # Host Plant / Family
  ancom_family_host <- ancombc2(
    data = ps_july_family,
    fix_formula = "host_plant", group = "host_plant",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  # Join the taxonomy names onto the ANCOM-BC2 results
  taxa <- as.data.frame(tax_table(ps_july_family)) |>
    rownames_to_column(var = "taxon")

  ancom_family_host_results <- ancom_family_host$res |>
    left_join(taxa, by = "taxon") |>
    # Filter for significant results
    filter(q_host_plantGRT <= 0.05) |>
    relocate(Family, .after = taxon) |>
    # Structural Zeros?
    left_join(ancom_family_host$zero_ind, by = "taxon") |>
    dplyr::select(-taxon, -Kingdom)

  ##########
  # Host Plant / ASV
  ancom_asv_host <- ancombc2(
    data = ps_july_asv,
    fix_formula = "host_plant", group = "host_plant",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  # Join the taxonomy names onto the ANCOM-BC2 results
  taxa <- as.data.frame(tax_table(ps_july_asv)) |>
    rownames_to_column(var = "taxon")

  ancom_asv_host_results <- ancom_asv_host$res |>
    left_join(taxa, by = "taxon") |>
    # Filter for significant results
    filter(q_host_plantGRT <= 0.05) |>
    # Structural Zeros?
    left_join(ancom_asv_host$zero_ind, by = "taxon") |>
    dplyr::select(-Kingdom)

  # Assemble all tables
  table_ancom_host <- bind_rows(
    clean_ancom_table(ancom_phylum_host_results, "Phylum", "Phylum", "host_plantGRT"),
    clean_ancom_table(ancom_family_host_results, "Family", "Family", "host_plantGRT"),
    clean_ancom_table(ancom_asv_host_results, "ASV", "taxon", "host_plantGRT")
  ) %>% arrange(Analysis_Level, q_val_FDR)

  table_ancom_host
}

# Month
{
  ##########
  # Month / Pylum
  ancom_phylum_month <- ancombc2(
    data = ps_grt_phylum,
    fix_formula = "month", group = "month",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  taxa <- as.data.frame(tax_table(ps_grt_phylum)) |>
    rownames_to_column(var = "taxon")

  ancom_phylum_month_results <- ancom_phylum_month$res |>
    left_join(taxa, by = "taxon") |>
    filter(q_monthJuly <= 0.05) |>
    relocate(Phylum, .after = taxon) |>
    left_join(ancom_phylum_month$zero_ind, by = "taxon") |>
    dplyr::select(-taxon, -Kingdom)

  ##########
  # Month / Family
  ancom_family_month <- ancombc2(
    data = ps_grt_family,
    fix_formula = "month", group = "month",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  taxa <- as.data.frame(tax_table(ps_grt_family)) |>
    rownames_to_column(var = "taxon")

  ancom_family_month_results <- ancom_family_month$res |>
    left_join(taxa, by = "taxon") |>
    filter(q_monthJuly <= 0.05) |>
    relocate(Family, .after = taxon) |>
    left_join(ancom_family_month$zero_ind, by = "taxon") |>
    dplyr::select(-taxon, -Kingdom)

  ##########
  # Month / ASV
  ancom_asv_month <- ancombc2(
    data = ps_grt_asv,
    fix_formula = "month", group = "month",
    p_adj_method = "holm",
    struc_zero = TRUE, neg_lb = TRUE,
    alpha = 0.05
  )

  taxa <- as.data.frame(tax_table(ps_grt_asv)) |>
    rownames_to_column(var = "taxon")

  ancom_asv_month_results <- ancom_asv_month$res |>
    left_join(taxa, by = "taxon") |>
    filter(q_monthJuly <= 0.05) |>
    relocate(Family, .after = taxon) |>
    left_join(ancom_asv_month$zero_ind, by = "taxon") |>
    dplyr::select(-Kingdom)

  # Assemble all tables
  table_ancom_month <- bind_rows(
    clean_ancom_table(ancom_phylum_month_results, "Phylum", "Phylum", "monthJuly"),
    clean_ancom_table(ancom_family_month_results, "Family", "Family", "monthJuly"),
    clean_ancom_table(ancom_asv_month_results, "ASV", "taxon", "monthJuly")
  ) %>% arrange(Analysis_Level, q_val_FDR)

  table_ancom_month
}

# intensive_development
{
  ##########
  # intensive_development / Pylum
  # Don't work for a continuous factor, since the number of taxa 
  # at the phylum level is too small for ANCOM 
  
  ##########
  # intensive_development / Family
  ancom_family_intensive_development <- ancombc2(
    data = ps_main_family,
    fix_formula = "intensive_development_ha", 
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    p_adj_method = "holm",
    alpha = 0.05
  )
  
  taxa <- as.data.frame(tax_table(ps_main_family)) |>
    rownames_to_column(var = "taxon")
  
  ancom_family_intensive_development_results <- ancom_family_intensive_development$res |>
    left_join(taxa, by = "taxon") |>
    filter(q_intensive_development_ha <= 0.05) |>
    relocate(Family, .after = taxon) |>
    dplyr::select(-taxon, -Kingdom)
  
  ##########
  # intensive_development / ASV
  ancom_asv_intensive_development <- ancombc2(
    data = ps_main_asv,
    fix_formula = "intensive_development_ha", 
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    p_adj_method = "holm",
    alpha = 0.05
  )
  
  taxa <- as.data.frame(tax_table(ps_main_asv)) |>
    rownames_to_column(var = "taxon")
  
  ancom_asv_intensive_development_results <- ancom_asv_intensive_development$res |>
    left_join(taxa, by = "taxon") |>
    filter(q_intensive_development_ha <= 0.05) |>
    relocate(Family, .after = taxon) |>
    dplyr::select(-Kingdom)
  
  # Assemble all tables
  table_ancom_intensive_development <- bind_rows(
    # clean_ancom_table(ancom_phylum_intensive_development_results, "Phylum", "Phylum", "intensive_development_ha"),
    clean_ancom_table(ancom_family_intensive_development_results, "Family", "Family", "intensive_development_ha"),
    clean_ancom_table(ancom_asv_intensive_development_results, "ASV", "taxon", "intensive_development_ha")
  ) %>% arrange(Analysis_Level, q_val_FDR)
  
  table_ancom_intensive_development
}

ancom.results <- bind_rows(
  table_ancom_host,
  table_ancom_month,
  table_ancom_intensive_development
) |> 
  select(
    Contrast, Analysis_Level, Feature_ID, LogFoldChange, StdError, W_stat, p_val, q_val_FDR,
    Phylum, Class, Order, Family, Genus, Species
  )

write.csv(ancom.results, "tables/TableS8.ANCOM.results.csv", row.names = FALSE)

################################################################################
# Save the key objects
save(
  ps_main_phylum, ps_main_phylum, ps_main_family, ps_main_family, ps_main_asv,
  ps_july_phylum, ps_july_phylum, ps_july_family, ps_july_family, ps_july_asv,
  ps_grt_phylum, ps_grt_phylum, ps_grt_family, ps_grt_family, ps_grt_asv,
  clean_ancom_table,
  ancom_phylum_host, ancom_phylum_host_results,
  ancom_family_host, ancom_family_host_results,
  ancom_asv_host, ancom_asv_host_results, table_ancom_host,
  ancom_phylum_month, ancom_phylum_month_results,
  ancom_family_month, ancom_family_month_results,
  ancom_asv_month, ancom_asv_month_results, table_ancom_month,
  # ancom_phylum_intensive_development, ancom_phylum_intensive_development_results,
  ancom_family_intensive_development, ancom_family_intensive_development_results,
  ancom_asv_intensive_development, ancom_asv_intensive_development_results, table_ancom_intensive_development,
  ancom.results,
  file = "ancom.RData"
)
