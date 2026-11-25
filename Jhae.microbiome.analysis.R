################################################################################
# Jadera microbiome analysis
################################################################################

# Analysis from the manuscript
#   "Host specificity, seasonality, and urbanicity shape the gut microbiome 
#    of a synanthropic insect"
#    
# by
#   Hanna Noyes, Nur-E-Jannat Anika, Josefine Just, Julia Combs, Tamheeda Khan, 
#   Rowan Graham Angelini, Su Fitzpatrick, and David R. Angelini
# 
# Colby College, Department of Biology, 5734 Mayflower Hill, Waterville, ME 04901, USA
# * Corresponding author email: david.r.angelini@gmail.com 
# 
# https://github.com/aphanotus/Jhae.microbiome

# Analysis begun 18 April 2026
#   with data processed by Hanna Noyes in Fall 2025

# Clear memory
# {
#   rm(list=ls())
#   .rs.restartR()
#   gc()
# }

setwd("~/Documents/3.research/Jhae.microbiome")

# Load packages
{
  library(tidyverse)
  library(ggrepel)
  library(phyloseq)
  library(patchwork)
  library(viridis)
}

# Load the image of all executed code below?
# load("Jhae.microbiome.RData", verbose = TRUE)

################################################################################
# Graphical settings
{
  theme_set(theme_bw())
  update_theme(legend.position = "none")
  update_theme(axis.ticks = element_blank())
  update_theme(strip.text = element_text(size = 12))
  
  color.by.sex   <- c("#F8766D","#02BFC3")
  color.by.morph <- c("darkgreen","darkred")
  color.by.month <- c("#4575b4","#FF8559")
  color.by.plant <- c("#9acd32","#DAA520")
  color.by.urban <- c("forestgreen","gray40")
  
  uuid.abbrev <- 8
}


################################################################################
# Load the data and perform initial processing
{
  # Load feature table
  {
    feature.table <- read.delim("qiime/feature.table.tsv", header = FALSE)
    colnames(feature.table) <- feature.table[2,]
    feature.table <- feature.table[-c(1:2),]
    colnames(feature.table)[1] <- "feature.id"
  }
  
  # Load the taxonomy assignments
  {
    taxonomy.lookup <- read.delim("qiime/silva138.taxonomy.tsv", header = TRUE)
    colnames(taxonomy.lookup) <- c("feature.id","taxonomy","confidence")
  }
  
  # Load the sample metadata
  {
    metadata <- read.csv("qiime/curated.metadata.csv", header = TRUE)
    
    # Manual correction of one sample name so that feature table and metadata match
    metadata <- metadata |> 
      mutate(sample.id = str_replace(sample.id, "FJ190108-017$", "FJ190108-017_plate1"))
    
    all((metadata$sample.id %in% colnames(feature.table)))
    all((colnames(feature.table)[-1] %in% metadata$sample.id))
    
    # Remove samples from the feature table for which we don't have metadata
    i <- which(!(colnames(feature.table) %in% metadata$sample.id))
    i <- i[-1]
    feature.table <- feature.table[,-i]
    
    # Calculate local intensive development
    metadata$intensive_development <-  with(metadata,
       developed_medium_intensity_m2 + 
       developed_high_intensity_m2
    )
    
    # Calculate local green space
    metadata$green_space <- with(metadata,
       developed_open_space_m2 + 
       scrub_m2 + grassland_m2 + 
       cropland_m2 + pasture_m2 + 
       herbaceous_wetlands_m2
    )
    
  }
}


################################################################################
# PhyloSeq formatting
{
  rownames(feature.table) <- NULL
  feature.table <- feature.table |>
    column_to_rownames(var = "feature.id") |>
    mutate(across(everything(), as.numeric)) |>
    as.matrix()
  
  ASV_obj <- otu_table(feature.table, taxa_are_rows = TRUE)
  
  rownames(metadata) <- NULL
  metadata <- metadata |> 
    column_to_rownames(var = "sample.id")
  
  META_obj <- sample_data(metadata)
  
  rownames(taxonomy.lookup) <- NULL
  cleaned_tax_matrix <- taxonomy.lookup |>
    # Normalize spacing (remove spaces after semicolons)
    mutate(taxonomy = str_replace_all(taxonomy, ";\\s*", ";")) |>
    
    # Split the single column into 7 ranks
    separate_wider_delim(
      cols = taxonomy, 
      delim = ";", 
      names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      too_few = "align_start", # If a row is missing ranks entirely, pad with NA
      too_many = "drop"        # If there are extra sub-ranks, drop them
    ) |>
    
    # Clean up the "d__", "p__", etc. prefixes using regular expressions
    mutate(across(Kingdom:Species, ~ str_remove(., "^[a-z]__"))) |>
    
    # Convert empty unassigned classifications (like old "s__") to NAs
    mutate(across(Kingdom:Species, ~ na_if(., ""))) |>
    
    # Move ASV IDs to row names and convert to a character matrix
    column_to_rownames(var = "feature.id") |> 
    dplyr::select(-confidence) |> 
    as.matrix()
  
  TAX_obj <- tax_table(cleaned_tax_matrix)
  
  # Combine the data
  ps <- phyloseq(ASV_obj, META_obj, TAX_obj)
  
  # Clean up
  rm(ASV_obj, META_obj, TAX_obj, cleaned_tax_matrix, metadata, feature.table)
}

# Filtration
{
  # Subset only adult Jadera haematoloma
  ps <- subset_samples(ps, ((stage == "adult") & (species == "Jhae")))
  ps 
  # 180 samples 
  
  # Include only samples with known sex
  ps <- subset_samples(ps, sex %in% c("female","male"))
  # 177 samples 

  # Record initial total read count
  sample_data(ps)$initial.reads <- sample_sums(ps)
  
  # Record the Wolbachia read counts
  ps_wolbachia <- subset_taxa(ps, Genus == "Wolbachia")
  sample_data(ps)$Wolbachia.reads <- sample_sums(ps_wolbachia)
  
  # Filter out Chloroplasts, Mitochondria, and other potential contaminants 
  # Note that the !(%in%) logic here is important to avoid rejecting ASVs 
  # with NA values at the taxonomic ranks being tested.
  ps <- subset_taxa(ps, Kingdom == "Bacteria")
  ps <- subset_taxa(
    ps, 
    !(Order %in% c("Chloroplast", "Mitochondria")) &
    !(Class %in% c("Chloroplast", "Mitochondria")) &
    !(Family %in% c("Chloroplast", "Mitochondria", "Propionibacteriaceae")) &
    !(Genus %in% "Wolbachia")
  )

  # Remove samples with fewer than 1000 total reads
  read.count.threshold <- 1000
  ps <- prune_samples(sample_sums(ps) >= read.count.threshold, ps)
  
  # Remove ASVs with zero counts across samples
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  ps_main <- ps
  
  # Extract metadata as a data frame
  metadata_main <- as(sample_data(ps_main), "data.frame")
}

# Save the collection site metadata
collection.sites <- metadata_main |> 
  filter(species=="Jhae", stage == "adult") |> 
  mutate(host_plant = str_replace(host_plant, "BV", "C. corindum")) |> 
  mutate(host_plant = str_replace(host_plant, "GRT", "K. paniculata")) |> 
  group_by(site_name) |> 
  summarise(
    Municipality = unique(city),
    Latitude = round(unique(latitude)[1],4),
    Longitude = round(unique(longitude)[1],4),
    Host_plant = unique(host_plant),
    Sample_size_in_January = sum(month == "January", na.rm = TRUE),
    Sample_size_in_July = sum(month == "July", na.rm = TRUE),
    Intensive_development_m2 = unique(intensive_development),
    Green_space_m2 = unique(green_space)
  ) |> 
  arrange(desc(Intensive_development_m2))

write.csv(collection.sites, "tables/TableS1.collection.sites.csv", row.names = FALSE)


# Read depth check by seasonality
{
  # All reads
  boxplot(log10(initial.reads) ~ month, data = metadata_main)
  wilcox.test(initial.reads ~ month, data = metadata_main)
  # W = 2145, p-value = 0.3168
  
  # Wolbachia only
  boxplot(log10(Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test(Wolbachia.reads ~ month, data = metadata_main)
  # W = 3302.5, p-value = 0.0001235
  
  # Reads without Wolbachia
  boxplot(log10(initial.reads - Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test((initial.reads - Wolbachia.reads) ~ month, data = metadata_main)
  # W = 1886, p-value = 0.0372
  # Interpretation: Wolbachia is significantly more abundant in January samples
  # although there is no difference in total initial read depth.
}

# Rarefaction
{ 
  rarefaction.replicates <- 1000
  
  # Rarefaction loop
  # Time difference of 8.18136 mins
  {
    start.time <- Sys.time()
    alpha_list <- lapply(1:rarefaction.replicates, function(i) {
      ps.i <- suppressMessages(
        rarefy_even_depth(ps_main, sample.size = read.count.threshold, replace = FALSE, rngseed = FALSE)
      )
      x <- estimate_richness(ps.i, measures = c("Observed", "Shannon", "Simpson"))
      return(x)
    })
    print(Sys.time() - start.time)
  }
  
  # Average across replicates
  # Reduce adds all the dataframes in the list together
  mean_alpha_diversity <- Reduce("+", alpha_list) / rarefaction.replicates
  
  # Add the alpha metrics to the main phyloseq object
  sample_data(ps_main)$observed_ASVs <- mean_alpha_diversity$Observed
  sample_data(ps_main)$shannon <- mean_alpha_diversity$Shannon
  sample_data(ps_main)$simpson <- mean_alpha_diversity$Simpson
  
  # Update the metadata to include the alpha metrics
  metadata_main <- as(sample_data(ps_main), "data.frame")
}


################################################################################
# Sub-setting

{
  # BV bugs were sampled in July so... 
  # - comparisons of month should be made in GRT samples
  # - comparisons of host_plant should be made in July samples
  ps_grt <- subset_samples(ps_main, host_plant == "GRT")
  ps_grt <- prune_taxa(taxa_sums(ps_grt) > 0, ps_grt)
  metadata_grt <- as(sample_data(ps_grt), "data.frame")

  ps_bv <- subset_samples(ps_main, host_plant == "BV")
  ps_bv <- prune_taxa(taxa_sums(ps_bv) > 0, ps_bv)
  metadata_bv <- as(sample_data(ps_bv), "data.frame")
  
  ps_july <- subset_samples(ps_main, month == "July")
  ps_july <- prune_taxa(taxa_sums(ps_july) > 0, ps_july)
  metadata_july <- as(sample_data(ps_july), "data.frame")
}


################################################################################
# Composition Bar Plot

{
  ps_main_family <- tax_glom(ps_main, taxrank = "Family")
  ps_rel <- transform_sample_counts(ps_main_family, function(x) x / sum(x, na.rm = TRUE))
  
  # Melt the PhyloSeq object into a tidy data frame for ggplot
  df_rel <- psmelt(ps_rel)
  
  # Find the top 9 most abundant Families overall
  top_families <- df_rel |>
    group_by(Family) |>
    summarize(Mean_Abund = mean(Abundance, na.rm = TRUE)) |>
    arrange(desc(Mean_Abund)) |>
    slice(1:9) |>
    pull(Family)
  
  # Collapse all other families into an "Other" category
  df_plot <- df_rel |>
    mutate(Family_Clean = ifelse(Family %in% top_families, as.character(Family), "Other")) |>
    # Ensure "Other" sits at the very bottom of the plot legend
    mutate(Family_Clean = factor(Family_Clean, levels = c(top_families, "Other")))
  
  # Arrange the data frame by metadata for visual order
  plot_composition <- df_plot |>
    mutate(host_plant = str_replace(host_plant,"BV","C. corindum")) |> 
    mutate(host_plant = str_replace(host_plant,"GRT","K. paniculata")) |> 
    arrange(month) |>
    mutate(Sample = factor(Sample, levels = unique(Sample))) |> 
    ggplot(aes(x = Sample, y = Abundance, fill = Family_Clean)) +
    geom_col(width = 1, linewidth = 0) + 
    facet_grid(~ host_plant * month, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = "Relative Abundance (%)",
      fill = "Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      legend.position = "bottom",
      # strip.background = element_rect(fill = "gray95", color = "gray80"),
      panel.spacing = unit(0.2, "lines") 
    )
  
  # plot_composition
  
  # Sum the 'Other' families within each individual sample first
  df_sample_agg <- df_plot |>
    group_by(Sample, host_plant, month, intensive_development, Family_Clean) |>
    summarize(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')
  
  # Group by host_plant, then calculate the average relative abundance
  plot_composition_grouped_by_host_plant <- df_sample_agg |>
    group_by(host_plant, Family_Clean) |>
    summarize(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = 'drop') |> 
    ggplot(aes(x = host_plant, y = Mean_Abundance, fill = Family_Clean)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.1) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = NULL,
      fill = "Bacterial Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      legend.position = "none",
      panel.spacing = unit(0.2, "lines") 
    )
  
  # plot_composition_grouped_by_host_plant
  
  # Group by month, then calculate the average relative abundance
  plot_composition_grouped_by_month <- df_sample_agg |>
    group_by(month, Family_Clean) |>
    summarize(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = 'drop') |> 
    ggplot(aes(x = month, y = Mean_Abundance, fill = Family_Clean)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.1) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = NULL,
      fill = "Bacterial Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      legend.position = "none",
      panel.spacing = unit(0.2, "lines") 
    )
  
  # plot_composition_grouped_by_month
  
  # Determine the median value to split samples by intensive_development
  median_dev <- median(df_sample_agg$intensive_development, na.rm = TRUE)
  
  # Assign "High" or "Low" intensive_development based on the median, summarize, and plot
  plot_composition_grouped_by_dev <- df_sample_agg |>
    mutate(dev_bin = ifelse(intensive_development > median_dev, "High", "Low")) |>
    mutate(dev_bin = factor(dev_bin, levels = c("Low", "High"))) |>
    group_by(dev_bin, Family_Clean) |>
    summarize(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = 'drop') |> 
    ggplot(aes(x = dev_bin, y = Mean_Abundance, fill = Family_Clean)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.1) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = NULL,
      fill = "Bacterial Family"
    ) +
    theme_minimal() +
    theme(
      # Keeping axis text here so "Low" and "High" are visible
      axis.text.x = element_text(angle = 0, hjust = 0.5), 
      panel.grid.major.x = element_blank(), 
      legend.position = "none",
      panel.spacing = unit(0.2, "lines") 
    )
  
  plot_composition_grouped_by_dev
  
  composition_figure <- (plot_composition | plot_composition_grouped_by_host_plant | plot_composition_grouped_by_month | plot_composition_grouped_by_dev) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(widths = c(9, 1, 1, 1), axes = "collect")
  
  print(composition_figure)
  
  ggsave(
    filename = "plots/Fig1.composition.barplots.png", 
    plot = composition_figure,
    width = 9, height = 3, scale = 1.75, 
    dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = "plots/Fig1.composition.barplots.pdf", 
    plot = composition_figure,
    width = 9, height = 3, scale = 1.75, 
    dpi = 300, bg = "white"
  )
  
  # An equivalent table
  ps_main_family <- tax_glom(ps_main, taxrank = "Family")
  ps_rel_family <- transform_sample_counts(ps_main_family, function(x) x / sum(x))
  family_mean_abundances <- taxa_sums(ps_rel_family) / nsamples(ps_rel_family)
  
  tax_table_df <- as(tax_table(ps_main_family), "matrix") |> 
    as.data.frame() |> 
    rownames_to_column(var = "ASV_ID")
  
  top_families_table <- data.frame(
    ASV_ID = names(family_mean_abundances),
    Mean_Relative_Abundance_Pct = family_mean_abundances * 100
  ) |>
    left_join(tax_table_df, by = "ASV_ID") |>
    arrange(desc(Mean_Relative_Abundance_Pct)) |>
    mutate(Mean_Relative_Abundance_Pct = round(Mean_Relative_Abundance_Pct, 2)) |>
    dplyr::select(Phylum, Class, Order, Family, Mean_Relative_Abundance_Pct)
  
  head(top_families_table, 10)
  
}


################################################################################
# Defining core microbiota
{
  # Calculate Mean Relative Abundance
  ps_bv_rel <- transform_sample_counts(ps_bv, function(x) x / sum(x, na.rm = TRUE))
  mean_abundance <- taxa_sums(ps_bv_rel) / nsamples(ps_bv_rel)
  
  # Calculate Prevalence
  # Convert raw counts to presence/absence (1 or 0), then calculate the percentage of samples
  ps_bv_pa <- transform_sample_counts(ps_bv, function(x) ifelse(x > 0, 1, 0))
  
  # Combine 
  bv_summary <- data.frame(
    ASV_ID = taxa_names(ps_bv),
    Mean_Abundance = mean_abundance,
    Prevalence = taxa_sums(ps_bv_pa) / nsamples(ps_bv_pa),
    as(tax_table(ps_bv), "matrix") # Attaches the taxonomy
  ) |> 
    dplyr::select(-Kingdom)
  
  # Sort the table by highest Mean Abundance (and then by Prevalence)
  bv_summary_sorted <- bv_summary[order(-bv_summary$Mean_Abundance, -bv_summary$Prevalence), ]
  
  # View the top 10 most dominant ASVs in the ancestral diet
  head(bv_summary_sorted, 10)
  
  # Mean Abundance and Prevalence
  ps_grt_rel <- transform_sample_counts(ps_grt, function(x) x / sum(x, na.rm = TRUE))
  mean_abund_grt <- taxa_sums(ps_grt_rel) / nsamples(ps_grt_rel)
  
  ps_grt_pa <- transform_sample_counts(ps_grt, function(x) ifelse(x > 0, 1, 0))
  
  grt_summary <- data.frame(
    ASV_ID = taxa_names(ps_grt),
    Mean_Abundance = mean_abund_grt,
    Prevalence = taxa_sums(ps_grt_pa) / nsamples(ps_grt_pa),
    as(tax_table(ps_grt), "matrix")
  ) |> 
    dplyr::select(-Kingdom)
  
  grt_summary_sorted <- grt_summary[order(-grt_summary$Mean_Abundance, -grt_summary$Prevalence), ]
  
  head(grt_summary_sorted, 10)

  # Define the core ASVs
  bv_summary_sorted$host <- "C. corindum"
  grt_summary_sorted$host <- "K. paniculata"
  
  core.microbiota <- bind_rows(bv_summary_sorted, grt_summary_sorted) |> 
    filter(Prevalence > 0.8, Mean_Abundance > 0.1) |> 
    remove_rownames() |> 
    relocate(host) 
  
  write.csv(core.microbiota, "tables/TableS3.core.microbiota.csv", row.names = FALSE)
}

##########
# Host shift 
{
  combined_summary <- bind_rows(bv_summary_sorted, grt_summary_sorted)
  
  # Extract the Top 10 from BV and GRT
  target_asvs <- unique(c(
    head(bv_summary_sorted$ASV_ID, 5), 
    head(grt_summary_sorted$ASV_ID, 5)
  ))
  
  df_for_plot_host_shift <- combined_summary |>
    filter(ASV_ID %in% target_asvs) |>
    mutate(
      Host = factor(host, levels = c("C. corindum", "K. paniculata")),
      Taxa_Label = case_when(
        !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
        !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
        !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
        !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
        !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
        TRUE ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
      )
    )
  
  plot_host_shift <- df_for_plot_host_shift |> 
    ggplot(aes(x = host, y = Mean_Abundance, group = Taxa_Label, color = Taxa_Label)) +
    geom_hline(yintercept = 0.1, linetype = "dotted", linewidth = 0.5, color = "grey35") +
    geom_line(linewidth = 1.5, alpha = 0.8) +
    geom_point(size = 4) +
    geom_text_repel(
      data = df_for_plot_host_shift |> filter(Host == "C. corindum"),
      aes(label = Taxa_Label), size = 2.5,
      nudge_x = -0.15, hjust = 1, direction = "y", fontface = "italic"
    ) +
    geom_text_repel(
      data = df_for_plot_host_shift |> filter(Host == "K. paniculata"),
      aes(label = Taxa_Label), size = 2.5,
      nudge_x = 0.15, hjust = 0, direction = "y", fontface = "italic"
    ) +
    scale_y_continuous(
      name = "Mean Relative Abundance",
      labels = scales::percent_format(), limits = c(0, NA)
    ) +
    scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, direction = -1) +
    # scale_color_viridis(option = "E", end = 0.7, discrete = TRUE, direction = -1) +
    # scale_color_brewer(palette = "Dark2") + 
    scale_x_discrete(name = NULL, expand = expansion(mult = c(1, 1))) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "italic", size = 12)
    )
  
  plot_host_shift
}

{
  # Are the highly abundant Alphaproteobacteria and Rhizobiaceae 
  # really endemic to BV and GRT?
  
  # First test this approach on an ASV, like the Enterobacteriaceae, 
  # that appear in both
  bind_rows(bv_summary_sorted, grt_summary_sorted) |> 
    filter(ASV_ID == bv_summary_sorted$ASV_ID[1])
  
  # Yes, the Alphaproteobacteria is endemic to BV
  bind_rows(bv_summary_sorted, grt_summary_sorted) |> 
    filter(ASV_ID == bv_summary_sorted$ASV_ID[2])
  
  # Yes, the Rhizobiaceae is endemic to GRT
  bind_rows(bv_summary_sorted, grt_summary_sorted) |> 
    filter(ASV_ID == grt_summary_sorted$ASV_ID[1])
}

################################################################################
# Prevalence vs. abundance for each host plant
{
  # Transform to relative abundance for the plot
  ps_july_rel <- transform_sample_counts(ps_july, function(x) x / sum(x, na.rm = TRUE))
  
  # Extract sample names for each host plant
  bv_samples <- sample_names(subset_samples(ps_july_rel, host_plant == "BV"))
  grt_samples <- sample_names(subset_samples(ps_july_rel, host_plant == "GRT"))
  
  # Extract the ASV matrix
  otu_matrix <- as(otu_table(ps_july_rel), "matrix")
  
  # Function to calculate prevalence and mean abundance
  prevalence.and.abundance <- function(mat, samples, host_name) {
    sub_mat <- mat[, samples, drop = FALSE]
    
    data.frame(
      ASV = rownames(sub_mat),
      Mean_Abundance = rowMeans(sub_mat),
      Prevalence = rowSums(sub_mat > 0) / length(samples),
      host_plant = host_name
    ) |>
      filter(Prevalence > 0) # Drop ASVs completely absent in this host
  }
  
  # Calculate metrics for both groups and combine
  df_bv <- prevalence.and.abundance(otu_matrix, bv_samples, "BV")
  df_grt <- prevalence.and.abundance(otu_matrix, grt_samples, "GRT")
  plot_df <- bind_rows(df_bv, df_grt) |> 
    mutate(host_plant = str_replace(host_plant, "BV", "C. corindum")) |> 
    mutate(host_plant = str_replace(host_plant, "GRT", "K. paniculata")) 
  
  # Extract taxonomy to create clean labels
  taxa_matrix <- as.data.frame(as(tax_table(ps_july_rel), "matrix")) |>
    rownames_to_column("ASV")
  
  plot_df <- plot_df |>
    left_join(taxa_matrix, by = "ASV") |>
    mutate(Taxa_Label = case_when(
      !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV, 1, uuid.abbrev), ")"),
      !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV, 1, uuid.abbrev), ")"),
      !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV, 1, uuid.abbrev), ")"),
      !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV, 1, uuid.abbrev), ")"),
      !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV, 1, uuid.abbrev), ")"),
      TRUE ~ paste0("Unclassified\n(", substr(ASV, 1, uuid.abbrev), ")")
    ))
  
  # Identify the top ASVs per host plant based on mean abundance
  top_asvs <- plot_df |>
    group_by(host_plant) |>
    slice_max(order_by = Mean_Abundance, n = 10) |>
    ungroup()
  
  # Color code the top ASVs
  plot_df <- plot_df |>
    mutate(highlight = ifelse(ASV %in% top_asvs$ASV, "Highlight", "Other"))
  
  # Build the plot
  pa_by_host_plot <- plot_df |> 
    ggplot(aes(x = Mean_Abundance, y = Prevalence)) +
    annotate(
      "rect", 
      xmin = 0.1, xmax = Inf, ymin = 0.8, ymax = Inf,
      fill = "lightgoldenrod1", alpha = 0.35
    ) +
    geom_point(aes(color = highlight, alpha = highlight), size = 2.5) +
    scale_color_manual(values = c("Highlight" = "#d95f02", "Other" = "grey60")) +
    scale_alpha_manual(values = c("Highlight" = 1, "Other" = 0.5)) +
    facet_wrap(~ host_plant) +
    scale_x_log10(labels = scales::percent_format(accuracy = 0.01)) +
    # scale_y_log10(labels = scales::percent_format(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_text_repel(
      data = top_asvs,
      aes(label = Taxa_Label),
      size = 2.5,
      fontface = "italic",
      box.padding = 0.5, point.padding = 0.3, 
      force = 4, max.overlaps = Inf, 
      color = "black", segment.color = "grey50"
    ) +
    theme(
      legend.position = "none",
      strip.background = element_blank(), 
      strip.text = element_text(size = 12, face = "italic")
    ) +
    labs(
      x = "Mean Relative Abundance (log10 scale)",
      y = "Prevalence (% of Samples)"
    ) +
    coord_flip()
  
  pa_by_host_plot
  
  figure_host_shift <- (pa_by_host_plot + plot_host_shift) +
    plot_annotation(tag_levels = 'A') 
  
  figure_host_shift
  
  ggsave(
    filename = "plots/Fig2.host.shift.png", 
    plot = figure_host_shift,
    width = 6.5, height = 3, scale = 1.6, 
    dpi = 600, bg = "white"
  )

}


################################################################################
# Save the image
save.image("Jhae.microbiome.RData")

# load("Jhae.microbiome.RData", verbose = TRUE)


################################################################################
# Follow-up scripts
# - "field.site.map.R"
# - "alpha.analysis.R"
# - "beta.analysis.R"
# - "GLM-PCA.R"
# - "ancom.R"
# - "aldex.R"
# - "MaAsLin3.R"
# - "consensus.DA.plot.R"
# - "permdisp.R"
