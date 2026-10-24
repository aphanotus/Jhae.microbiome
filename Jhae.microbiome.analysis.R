# Jadera microbiome analysis
# 
# Begun 18 April 2026
#   with data processed by Hanna Noyes in Fall 2025
# Completed TBD

# Clear memory
{
  rm(list=ls())
  .rs.restartR()
  gc()
}

setwd("~/Documents/3.research/Jhae.microbiome")

# Load packages
{
  library(tidyverse)
  library(ggrepel)
  library(phyloseq)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(patchwork)
  library(viridis)
}

################################################################################
# Load the data and process metrics
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
    # metadata <- read.delim("qiime/curated.metadata.tsv", header = TRUE)
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
    
    # Calculate the proportion of local green space
    metadata$green_space <- (
      metadata$developed_open_space_m2 +
        metadata$deciduous_forest_m2 +
        metadata$evergreen_forest_m2 +
        metadata$mixed_forest_m2 +
        metadata$cropland_m2
    ) / (1e6 - metadata$open_water_m2)
    
    # Determine an urban index from the land surface data
    {
      # Focus on Jadera and extract environmental variables
      env_data <- metadata[which(metadata$species=="Jhae"),c(9,21:35)]
      env_data <- unique(env_data)
      
      # Set the row names to the Site IDs (this makes the biplot labels much cleaner)
      rownames(env_data) <- env_data$site_name
      env_data <- env_data[,-1]
      
      # Drop any columns that have zero variance (e.g., 0% everywhere) before running PCA
      env_data <- env_data[, apply(env_data, 2, var) != 0]
      
      # Run PCA
      pca_env <- prcomp(env_data, scale. = TRUE)
      
      # Extract the PC1 scores for each site
      site_pc_scores <- data.frame(
        site_name = rownames(pca_env$x),
        urbanicity = pca_env$x[, "PC1"]
      )
      
      # Extract the PC1 loadings (rotation) into a new dataframe for plotting
      pc1_loadings <- data.frame(
        feature = rownames(pca_env$rotation),
        loading = pca_env$rotation[, "PC1"]
      )
      
      # # scale = 0 ensures the arrows and data points are plotted on the same relative scale
      # biplot(pca_env, scale = 0, cex = c(0.5, 0.8))
      # 
      # # Plot the loadings as a sorted bar chart
      # ggplot(pc1_loadings, aes(x = reorder(feature, loading), y = loading)) +
      #   geom_col() +
      #   coord_flip() +
      #   theme_bw()
      
      # Map the site-level PC1 scores back to the individual bugs in the metadata
      # Use match() to ensure the row order
      # Change the direction so higher values are more urban
      metadata$urbanicity <- -1*site_pc_scores$urbanicity[match(metadata$site_name, site_pc_scores$site_name)]
    }
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
    # Normalize spacing (removes any accidental spaces after semicolons)
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
  rm(ASV_obj, META_obj, TAX_obj)
  
  # Subset only adult Jadera haematoloma
  ps <- subset_samples(ps, ((stage == "adult") & (species == "Jhae")))
  
  # Include only samples with known sex
  ps <- subset_samples(ps, sex %in% c("female","male"))

  # Record total read count
  sample_data(ps)$read.count <- sample_sums(ps)
  
  # Retain the Wolbachia read counts
  ps_wolbachia <- subset_taxa(ps, Genus == "Wolbachia")
  sample_data(ps)$Wolbachia.reads <- sample_sums(ps_wolbachia)
  
  # Filter out Chloroplasts, Mitochondria, Wolbachia and other potential contaminants 
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
  
  # Remove ASVs with zero counts across all samples.
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Remove samples with zero counts across all ASVs
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  # Alpha diversity metrics
  richness_df <- estimate_richness(ps, measures = c("Observed","Shannon","Simpson"))
  
  # Add the Observed ASV number 
  sample_data(ps)$observed_ASVs <- richness_df$Observed
  sample_data(ps)$shannon <- richness_df$Shannon
  sample_data(ps)$simpson <- richness_df$Simpson
  
  # Extract the metadata into a data frame
  ps_main <- ps
  metadata_main <- as(sample_data(ps_main), "data.frame")
}

# Rarefaction
{ 
  # {
  #   # Extract the ASV table
  #   # (Ensure taxa are columns and samples are rows)
  #   otu_table_ps <- as(otu_table(ps), "matrix")
  #   if (taxa_are_rows(ps)) {
  #     otu_table_ps <- t(otu_table_ps)
  #   }
  # 
  #   # Plot the rarefaction curve, zooming in on the 0-500 read range
  #   rarecurve(
  #     otu_table_ps,
  #     step = 50,
  #     xlim = c(0, 2000), ylim = c(0,100),
  #     main = "Rarefaction Curves (Wolbachia Removed)",
  #     xlab = "Sequencing Depth", ylab = "Number of ASVs"
  #   )
  #   abline(v=c(100,500,1000, 2000), lty = 2, col = "darkred")
  # }
  # 
  # x <- prune_samples(sample_sums(ps) >= 100, ps)
  # prune_taxa(taxa_sums(x) > 0, x)
  # x <- prune_samples(sample_sums(ps) >= 500, ps)
  # prune_taxa(taxa_sums(x) > 0, x)
  # x <- prune_samples(sample_sums(ps) >= 1000, ps)
  # prune_taxa(taxa_sums(x) > 0, x)
  # x <- prune_samples(sample_sums(ps) >= 2000, ps)
  # prune_taxa(taxa_sums(x) > 0, x)
  
  # | read min | taxa |       | samples |       |
  # | -------- | ---- | ----- | ------- | ----- |
  # | NA       | 2148 |    NA | 175     |    NA |
  # | 100      | 2133 | 99.3% | 163     | 93.1% |
  # | 500      | 2101 | 97.8% | 153     | 87.4% |
  # | 1000     | 2064 | 96.0% | 144     | 82.3% |
  # | 2000     | 2055 | 95.7% | 138     | 78.9% |
  
  # Prune out samples with fewer than 1000 total reads
  ps_rare <- prune_samples(sample_sums(ps) >= 1000, ps)

  # Clean up any ASVs that now have zero counts because their samples were dropped
  ps_rare <- prune_taxa(taxa_sums(ps_rare) > 0, ps_rare)

  # Extract the metadata into a standard data frame
  metadata_rare <- as(sample_data(ps_rare), "data.frame")
}

# Read depth check by seasonality
{
  boxplot(log10(read.count+1) ~ month, data = metadata_main)
  wilcox.test(read.count ~ month, data = metadata_main)
  # W = 2865, p-value = 0.002457
  
  boxplot(log10(Wolbachia.reads+1) ~ month, data = metadata_main)
  wilcox.test(Wolbachia.reads ~ month, data = metadata_main)
  # W = 5291, p-value = 3.865e-05
  
  boxplot(log10(read.count + Wolbachia.reads+1) ~ month, data = metadata_main)
  wilcox.test((read.count + Wolbachia.reads) ~ month, data = metadata_main)
  # W = 3157, p-value = 0.03004
  # Interpretation: Wolbachia is significantly more abundant in January samples
  # although there is only a marginal difference in total read counts.
}


################################################################################
# Subsetting

{
  # Different taxonomic levels
  ps_phylum <- tax_glom(ps_main, taxrank = "Phylum")
  ps_family <- tax_glom(ps_main, taxrank = "Family")
  ps_asv <- ps_main
  
  # BV bugs were sampled in July so... 
  # - comparisons of month should be made in GRT samples
  # - comparisons of host_plant should be made in July samples
  ps_grt <- subset_samples(ps_main, host_plant == "GRT")
  ps_grt <- prune_taxa(taxa_sums(ps_grt) > 0, ps_grt)
  metadata_grt <- as(sample_data(ps_grt), "data.frame")

  ps_rare_grt <- subset_samples(ps_rare, host_plant == "GRT")
  ps_rare_grt <- prune_taxa(taxa_sums(ps_rare_grt) > 0, ps_rare_grt)
  metadata_rare_grt <- as(sample_data(ps_rare_grt), "data.frame")
  
  ps_bv <- subset_samples(ps_main, host_plant == "BV")
  ps_bv <- prune_taxa(taxa_sums(ps_bv) > 0, ps_bv)
  metadata_bv <- as(sample_data(ps_bv), "data.frame")
  
  ps_july <- subset_samples(ps_main, month == "July")
  ps_july <- prune_taxa(taxa_sums(ps_july) > 0, ps_july)
  metadata_july <- as(sample_data(ps_july), "data.frame")

  ps_rare_july <- subset_samples(ps_rare, month == "July")
  ps_rare_july <- prune_taxa(taxa_sums(ps_rare_july) > 0, ps_rare_july)
  metadata_rare_july <- as(sample_data(ps_rare_july), "data.frame")
  
  # Isolate the July Goldenrain Tree samples
  ps_grt_july <- subset_samples(ps_main, host_plant == "GRT" & month == "July")
  ps_grt_july <- prune_taxa(taxa_sums(ps_grt_july) > 0, ps_grt_july)
  
}


################################################################################
# Graphical settings
{
  theme_set(theme_bw())
  update_theme(legend.position = "none")
  update_theme(axis.ticks = element_blank())
  update_theme(strip.text = element_text(size = 12))

  color.by.sex = c("#F8766D","#02BFC3")
  color.by.morph = c("darkgreen","darkred")
  color.by.month = c("#4575b4","#FF8559")
  color.by.plant = c("#9acd32","#DAA520")
  color.by.urban = c("gray40","forestgreen")
  
  uuid.abbrev <- 8
  
}


################################################################################
# Composition Bar Plot

{
  ps_rel <- transform_sample_counts(ps_family, function(x) x / sum(x, na.rm = TRUE))
  
  # Melt the phyloseq object into a tidy data frame for ggplot
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
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      panel.spacing = unit(0.2, "lines") 
    )
  
  # plot_composition
  
  # Group by host_plant, then calculate the average relative abundance
  plot_composition_grouped_by_host_plant <- df_plot |>
    group_by(host_plant, Family_Clean) |>
    summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') |> 
    ggplot(aes(x = host_plant, y = Mean_Abundance, fill = Family_Clean)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.1) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = NULL,
      fill = "Bacterial Family"
    ) +
    theme_minimal() 
  
  # plot_composition_grouped_by_host_plant
  
  # Group by month, then calculate the average relative abundance
  plot_composition_grouped_by_month <- df_plot |>
    group_by(month, Family_Clean) |>
    summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') |> 
    ggplot(aes(x = month, y = Mean_Abundance, fill = Family_Clean)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.1) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
    scale_fill_viridis(discrete = TRUE) +
    labs(
      x = NULL, y = NULL,
      fill = "Bacterial Family"
    ) +
    theme_minimal() 
  
  # plot_composition_grouped_by_month
  
  composition_figure <- (plot_composition | plot_composition_grouped_by_host_plant | plot_composition_grouped_by_month) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(widths = c(9, 1, 1), guides = "collect", axes = "collect")
  
  composition_figure
  
  ggsave(
    filename = "plots/FigX.composition.barplots.png", 
    plot = composition_figure,
    width = 9, height = 3, scale = 1.75, 
    dpi = 300, bg = "white"
  )
  
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
  ps_grt_july_rel <- transform_sample_counts(ps_grt_july, function(x) x / sum(x, na.rm = TRUE))
  mean_abund_grt_july <- taxa_sums(ps_grt_july_rel) / nsamples(ps_grt_july_rel)
  
  ps_grt_july_pa <- transform_sample_counts(ps_grt_july, function(x) ifelse(x > 0, 1, 0))
  
  grt_july_summary <- data.frame(
    ASV_ID = taxa_names(ps_grt_july),
    Mean_Abundance = mean_abund_grt_july,
    Prevalence = taxa_sums(ps_grt_july_pa) / nsamples(ps_grt_july_pa),
    as(tax_table(ps_grt_july), "matrix")
  ) |> 
    dplyr::select(-Kingdom)
  
  grt_july_summary_sorted <- grt_july_summary[order(-grt_july_summary$Mean_Abundance, -grt_july_summary$Prevalence), ]
  
  head(grt_july_summary_sorted, 10)
}

##########
# Host shift 
{
  bv_summary_sorted$host <- "C. corindum"
  grt_july_summary_sorted$host <- "K. paniculata"
  
  combined_july_summary <- bind_rows(bv_summary_sorted, grt_july_summary_sorted)
  
  # Extract the Top 10 from BV and July GRT
  target_asvs_july <- unique(c(
    head(bv_summary_sorted$ASV_ID, 5), 
    head(grt_july_summary_sorted$ASV_ID, 5)
  ))
  
  slope_df_july <- combined_july_summary |>
    filter(ASV_ID %in% target_asvs_july) |>
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
  
  plot_host_shift <- slope_df_july |> 
    ggplot(aes(x = host, y = Mean_Abundance, group = Taxa_Label, color = Taxa_Label)) +
    geom_line(linewidth = 1.5, alpha = 0.8) +
    geom_point(size = 4) +
    geom_text_repel(
      data = slope_df_july |> filter(Host == "C. corindum"),
      aes(label = Taxa_Label),
      nudge_x = -0.15, hjust = 1, direction = "y", fontface = "italic", size = 4
    ) +
    geom_text_repel(
      data = slope_df_july |> filter(Host == "K. paniculata"),
      aes(label = Taxa_Label),
      nudge_x = 0.15, hjust = 0, direction = "y", fontface = "italic", size = 4
    ) +
    scale_y_continuous(
      name = "Mean Relative Abundance",
      labels = scales::percent_format(), limits = c(0, NA)
    ) +
    scale_color_viridis(option = "D", end = 0.7, discrete = TRUE, direction = -1) +
    # scale_color_viridis(option = "E", end = 0.7, discrete = TRUE, direction = -1) +
    # scale_color_brewer(palette = "Dark2") + 
    scale_x_discrete(name = NULL, expand = expansion(mult = c(1, 1))) +
    # scale_x_discrete(name = NULL) + 
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "italic", size = 12)
    )
  
  plot_host_shift
  
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
      size = 3.5,
      fontface = "italic",
      box.padding = 0.5, point.padding = 0.3, 
      force = 4, max.overlaps = Inf, 
      color = "black",
      segment.color = "grey50"
    ) +
    theme(
      legend.position = "none",
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
    filename = "plots/FigX.host.shift.png", 
    plot = figure_host_shift,
    width = 6.5, height = 3.5, scale = 2.25, 
    dpi = 600, bg = "white"
  )

}


################################################################################
# Alpha diversity

{
  alpha.columns <- c("observed_ASVs", "simpson") # "shannon", 
  alpha.labels <- c("Observed ASVs", "Simpson Index") # "Shannon Index", 
  names(alpha.labels) <- alpha.columns
  
  # Alpha diversity by sex
  plot_alpha_sex <- metadata_rare |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=sex, y=value, fill = sex)) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.sex) +
    labs(x = NULL, y = "alpha diversity index", title = "Sex") 
  
  plot_alpha_sex
  
  # Alpha diversity by morph
  plot_alpha_morph <- metadata_rare |> 
    filter(morph %in% c("LW","SW")) |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=morph, y=value, fill = morph)) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.morph) +
    labs(x = NULL, y = "alpha diversity index", title = "Wing morph") 
  
  plot_alpha_morph
  
  # Alpha diversity by host plant   
  plot_alpha_host <- metadata_rare |> 
    filter(month == "July") |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=host_plant, y=value, fill = host_plant)) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.plant) +
    labs(x = NULL, y = "alpha diversity index", title = "Host Plant") 
  
  plot_alpha_host
  
  # Alpha diversity by month
  plot_alpha_month <- metadata_rare |> 
    filter(host_plant == "GRT") |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=month, y=value, fill = month)) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.month) +
    labs(x = NULL, y = "alpha diversity index", title = "Month") 
  
  plot_alpha_month
  
  # Alpha diversity by urbanicity
  plot_alpha_urban_Obs <- metadata_rare |> 
    ggplot(aes(x = urbanicity, y = observed_ASVs)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = urbanicity), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = "Urbanicity PC1", y = "Observed ASVs", title = "Urbanicity") 
  
  plot_alpha_urban_Obs
  
  plot_alpha_urban_SI <- metadata_rare |> 
    ggplot(aes(x = urbanicity, y = simpson)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = urbanicity), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = "Urbanicity PC1", y = "Simpson's index", title = "Urbanicity") 
  
  plot_alpha_urban_SI
  
  # Alpha diversity by green space
  plot_alpha_green_Obs <- metadata_rare |> 
    ggplot(aes(x = green_space, y = observed_ASVs)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = green_space), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    xlim(c(0,0.8)) +
    labs(x = "Local green space (m² / km²)", y = "Observed ASVs", title = "Green space") 
  
  plot_alpha_green_Obs

    plot_alpha_green_SI <- metadata_rare |> 
    ggplot(aes(x = green_space, y = simpson)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = green_space), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    xlim(c(0,0.8)) +
    labs(x = "Local green space (m² / km²)", y = "Simpson's index", title = "Green space") 
  
  plot_alpha_green_SI
}

##########
# Alpha diversity stats
{
  # Using Wilcoxon rank sum test
  
  # sex 
  wilcox.test(observed_ASVs ~ sex, data = metadata_rare)
  # W = 1921, p-value = 0.007405
  wilcox.test(simpson ~ sex, data = metadata_rare)
  # NS
  
  # morph 
  x <- metadata_rare |> filter(morph %in% c("LW","SW"))
  wilcox.test(observed_ASVs ~ morph, data = x)
  # NS
  wilcox.test(simpson ~ morph, data = x)
  # NS
  
  # host plant
  x <- metadata_rare |> filter(month == "July")
  wilcox.test(observed_ASVs ~ host_plant, data = x)
  # NS
  wilcox.test(simpson ~ host_plant, data = x)
  # W = 359, p-value = 0.01012
  
  # month (January vs July) ***
  x <- metadata_rare |> filter(host_plant == "GRT")
  wilcox.test(observed_ASVs ~ month, data = x)
  # W = 1442, p-value = 0.0002286
  wilcox.test(simpson ~ month, data = x)
  # W = 1316, p-value = 2.353e-05
  
  # LMM for urbanicity
  urban_model <- lmer(observed_ASVs ~ urbanicity + (1 | site_name), data = metadata_rare)
  summary(urban_model)
  # NS
  urban_model <- lmer(simpson ~ urbanicity + (1 | site_name), data = metadata_rare)
  summary(urban_model)
  # NS
  
  # LMM for green space
  green_space_model <- lmer(observed_ASVs ~ green_space + (1 | site_name), data = metadata_rare)
  summary(green_space_model)
  # NS
  green_space_model <- lmer(simpson ~ green_space + (1 | site_name), data = metadata_rare)
  summary(green_space_model)
  # NS  
}


##########
# Assembling figures

{
  alpha_figure_NS <- (
    plot_alpha_morph /
    (plot_alpha_urban_Obs | plot_alpha_urban_SI) /
    (plot_alpha_green_Obs | plot_alpha_green_SI)
  ) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect")
  alpha_figure_NS
  
  ggsave(
    filename = "plots/FigX.alpha.diversity.NS.png", 
    plot = alpha_figure_NS,
    width = (2.25)*2, height = (2.0)*3, scale = 1.8, 
    dpi = 600, bg = "white"
  )
  
  alpha_figure <- (
    plot_alpha_sex | plot_alpha_host | plot_alpha_month 
  ) +
    plot_annotation(tag_levels = 'A')
  alpha_figure
  
  ggsave(
    filename = "plots/FigX.alpha.diversity.png", 
    plot = alpha_figure,
    width = (2.5)*3, height = 2.5, scale = 1.8, 
    dpi = 600, bg = "white"
  )
  
}


################################################################################
# Distance comparisons

# Distance matrices and ordinations
{
  # Bray-Curtis distance
  dist_rare_bray <- distance(ps_rare, method = "bray")
  dist_july_bray <- distance(ps_rare_july, method = "bray")
  dist_grt_bray <- distance(ps_rare_grt, method = "bray")
  
  ord_july_bray <- ordinate(ps_rare_july, method = "PCoA", distance = dist_july_bray)
  ord_grt_bray <- ordinate(ps_rare_grt, method = "PCoA", distance = dist_grt_bray)
  
  # UniFrac distances
  asv.tree <- phytools::read.newick("qiime/16s.Jhae.rooted-tree.tre")
  ps_rare_tree <- merge_phyloseq(ps_rare, asv.tree)
  ps_july_tree <- merge_phyloseq(ps_rare_july, asv.tree)
  ps_grt_tree <- merge_phyloseq(ps_rare_grt, asv.tree)
  
  dist_rare_uwUF <- distance(ps_rare_tree, method = "unifrac")
  dist_july_uwUF <- distance(ps_july_tree, method = "unifrac")
  dist_grt_uwUF <- distance(ps_grt_tree, method = "unifrac")
  
  dist_rare_wUF <- distance(ps_rare_tree, method = "wunifrac")
  dist_july_wUF <- distance(ps_july_tree, method = "wunifrac")
  dist_grt_wUF <- distance(ps_grt_tree, method = "wunifrac")
  
  ord_rare_uwUF <- ordinate(ps_rare_tree, method = "PCoA", distance = dist_rare_uwUF)
  ord_july_uwUF <- ordinate(ps_july_tree, method = "PCoA", distance = dist_july_uwUF)
  ord_grt_uwUF <- ordinate(ps_grt_tree, method = "PCoA", distance = dist_grt_uwUF)
  
  ord_rare_wUF <- ordinate(ps_rare_tree, method = "PCoA", distance = dist_rare_wUF)
  ord_july_wUF <- ordinate(ps_july_tree, method = "PCoA", distance = dist_july_wUF)
  ord_grt_wUF <- ordinate(ps_grt_tree, method = "PCoA", distance = dist_grt_wUF)
}

# PCoA Plots
{
  pcoa.bray.host <- plot_ordination(ps_july, ord_july_bray, color = "host_plant") +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(name = "Host plant", values = color.by.plant) +
    coord_fixed()
  pcoa.bray.host
  
  pcoa.uwUF.host <- plot_ordination(ps_july_tree, ord_july_uwUF, color = "host_plant") +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = color.by.plant) +
    coord_fixed()
  pcoa.uwUF.host
  
  pcoa.wUF.host <- plot_ordination(ps_july_tree, ord_july_wUF, color = "host_plant") +
    theme_classic() +
    theme(legend.position = "right") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = color.by.plant) +
    coord_fixed()
  pcoa.wUF.host

  figure_PCoA_host <- (
    (pcoa.bray.host | pcoa.uwUF.host) / pcoa.wUF.host
  ) +
    plot_annotation(tag_levels = 'A') 

  figure_PCoA_host
  
  ggsave(
    filename = "plots/FigX.PCoA.by.host.png", 
    plot = figure_PCoA_host,
    width = (2.5)*2, height = (2.5)*2, scale = 1.25, 
    dpi = 600, bg = "white"
  )
  
}


##########
# PERMANVOA Assumptions

# Testing dispersion of each factor
dispersion_test <- betadisper(dist_rare_bray, metadata_rare$sex)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_rare_bray, metadata_rare$morph)
permutest(dispersion_test) # 0.013 *
dispersion_test <- betadisper(dist_july_bray, metadata_rare_july$host_plant)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_grt_bray, metadata_rare_grt$month)
permutest(dispersion_test) # 0.052 .

dispersion_test <- betadisper(dist_rare_uwUF, metadata_rare$sex)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_rare_uwUF, metadata_rare$morph)
permutest(dispersion_test) # 0.001 ***
dispersion_test <- betadisper(dist_july_uwUF, metadata_rare_july$host_plant)
permutest(dispersion_test) # 0.001 ***
dispersion_test <- betadisper(dist_grt_uwUF, metadata_rare_grt$month)
permutest(dispersion_test) # NS

# Look in more detail ar dispersion by host_plant
dispersion_by_host_plant <- betadisper(dist_july_uwUF, metadata_rare_july$host_plant)
dispersion_by_host_plant <- data.frame(
  Group = dispersion_by_host_plant$group,
  Distance = dispersion_by_host_plant$distances
)
plot.dispersion_by_host_plant <- ggplot(dispersion_by_host_plant, aes(x = Group, y = Distance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) + # Shows the actual samples
  labs(
    x = "Host plant",
    y = "Distance to Centroid (Beta Dispersion)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Hide legend since x-axis has the labels

plot.dispersion_by_host_plant

dispersion_test <- betadisper(dist_rare_wUF, metadata_rare$sex)
permutest(dispersion_test) # 0.02 *
dispersion_test <- betadisper(dist_rare_wUF, metadata_rare$morph)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_july_wUF, metadata_rare_july$host_plant)
permutest(dispersion_test) # 0.023 *
dispersion_test <- betadisper(dist_grt_wUF, metadata_rare_grt$month)
permutest(dispersion_test) # NS

# Look in more detail ar dispersion by host_plant
dispersion_by_host_plant <- betadisper(dist_july_wUF, metadata_rare_july$host_plant)
dispersion_by_host_plant <- data.frame(
  Group = dispersion_by_host_plant$group,
  Distance = dispersion_by_host_plant$distances
)
plot.dispersion_by_host_plant <- ggplot(dispersion_by_host_plant, aes(x = Group, y = Distance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) + # Shows the actual samples
  labs(
    x = "Host plant",
    y = "Distance to Centroid (Beta Dispersion)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Hide legend since x-axis has the labels

plot.dispersion_by_host_plant


##########
# PERMANOVA Models 
# Partial sums of squares (the equivalent of Type III ANOVA)
iterations <- 1e3

# Full Model
perm_bray_full <- adonis2(
  dist_rare_bray ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_rare, 
  by = "margin",  
  permutations = iterations
)

perm_bray_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1    0.362 0.00821 1.2484 0.218781    
# morph        2    0.599 0.01358 1.0324 0.373626    
# host_plant   1    1.965 0.04457 6.7744 0.000999 ***
# month        1    0.953 0.02161 3.2851 0.008991 ** 
# urbanicity   1    0.562 0.01276 1.9393 0.080919 .  
# Residual   137   39.729 0.90125                    
# Total      143   44.083 1.00000   

perm_uwUF_full <- adonis2(
  dist_rare_uwUF ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_rare, 
  by = "margin",  
  permutations = iterations
)

perm_uwUF_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1    0.550 0.00863 1.2888 0.052947 .  
# morph        2    0.917 0.01439 1.0747 0.228771    
# host_plant   1    1.809 0.02838 4.2381 0.000999 ***
# month        1    0.929 0.01458 2.1772 0.000999 ***
# urbanicity   1    0.821 0.01289 1.9246 0.000999 ***
# Residual   137   58.463 0.91752                    
# Total      143   63.719 1.00000 

perm_wUF_full <- adonis2(
  dist_rare_wUF ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_rare, 
  by = "margin",  
  permutations = iterations
)

perm_wUF_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1 0.0001189 0.00410  0.9731 0.453546    
# morph        2 0.0000802 0.00276  0.3280 0.621379    
# host_plant   1 0.0094116 0.32446 77.0139 0.000999 ***
# month        1 0.0001492 0.00515  1.2213 0.389610    
# urbanicity   1 0.0000962 0.00332  0.7872 0.481518    
# Residual   137 0.0167423 0.57718                     
# Total      143 0.0290072 1.00000    


################################################################################
# Save the image
save.image("Jhae.microbiome.RData")


################################################################################
# Table generation

# Collection site metadata
x <- metadata_main |> 
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
    Proportion_of_green_space = round(unique(green_space),4),
    Urbanicity_index = round(unique(urbanicity),3)
  ) |> 
  arrange(Urbanicity_index)

write_delim(x, "collection.sites.tsv", delim = "\t")


################################################################################
# Follow-up scripts
# - "field.site.map.R"
# - "GLM-PCA.R"
# - "ancom.R"
# - "aldex.R"
# - "fastEmu.R"
# - "MaAsLin3.R"

