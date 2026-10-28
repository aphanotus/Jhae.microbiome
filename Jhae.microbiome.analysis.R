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
# Load the data and perform inital processing
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
      
      # Clean up
      rm(env_data, pca_env, site_pc_scores, i, pc1_loadings)
      
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
  
  # Remove samples with fewer than 1000 total reads
  read.count.threshold <- 1000
  ps <- prune_samples(sample_sums(ps) >= read.count.threshold, ps)
  
  # Remove ASVs with zero counts across samples
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  ps_main <- ps
  
  # Extract metadata as a data frame
  metadata_main <- as(sample_data(ps_main), "data.frame")
}

# Read depth check by seasonality
{
  # All reads
  boxplot(log10(initial.reads) ~ month, data = metadata_main)
  wilcox.test(initial.reads ~ month, data = metadata_main)
  # W = 2219, p-value = 0.385
  
  # Wolbachia only
  boxplot(log10(Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test(Wolbachia.reads ~ month, data = metadata_main)
  # W = 3348.5, p-value = 0.0001481
  
  # Reads without Wolbachia
  boxplot(log10(initial.reads - Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test((initial.reads - Wolbachia.reads) ~ month, data = metadata_main)
  # W = 1960, p-value = 0.05269
  # Interpretation: Wolbachia is significantly more abundant in January samples
  # although there is no difference in total initial read depth.
}

# Rarefaction
{ 
  # Sampling without replacement
  ps_rare <- rarefy_even_depth(
    ps_main, 
    sample.size = read.count.threshold, 
    rngseed = 123,
    replace = FALSE,   
    trimOTUs = TRUE,
    verbose = TRUE
  )
  # 669 taxa and 144 samples
  
  summary(colSums(otu_table(ps_main)))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1161   19739   60716   85626  131397  371922 

  summary(colSums(otu_table(ps_rare)))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1000    1000    1000    1000    1000    1000 

  # Alpha diversity metrics
  richness_df <- estimate_richness(ps_rare, measures = c("Observed","Shannon","Simpson"))
  
  # Add the Observed ASV number 
  sample_data(ps_rare)$observed_ASVs <- richness_df$Observed
  sample_data(ps_rare)$shannon <- richness_df$Shannon
  sample_data(ps_rare)$simpson <- richness_df$Simpson
  
  # Extract metadata
  # This will be identical to the "main" metadata DF, but include alpha metrics
  metadata_rare <- as(sample_data(ps_rare), "data.frame")
  
}


################################################################################
# Subsetting

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

  ps_rare_grt <- subset_samples(ps_rare, host_plant == "GRT")
  metadata_rare_grt <- as(sample_data(ps_rare_grt), "data.frame")
  
  ps_rare_july <- subset_samples(ps_rare, month == "July")
  metadata_rare_july <- as(sample_data(ps_rare_july), "data.frame")
  
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
  ps_main_family <- tax_glom(ps_main, taxrank = "Family")
  ps_rel <- transform_sample_counts(ps_main_family, function(x) x / sum(x, na.rm = TRUE))
  
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
      legend.position = "bottom",
      # strip.background = element_rect(fill = "gray95", color = "gray80"),
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
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      legend.position = "none",
      panel.spacing = unit(0.2, "lines") 
    )
  
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
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(), 
      legend.position = "none",
      panel.spacing = unit(0.2, "lines") 
    )
  
  # plot_composition_grouped_by_month
  
  composition_figure <- (plot_composition | plot_composition_grouped_by_host_plant | plot_composition_grouped_by_month) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(widths = c(9, 1, 1), axes = "collect")
  
  print(composition_figure)
  
  ggsave(
    filename = "plots/Fig1.composition.barplots.png", 
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
  
  write.csv(core.microbiota, "tables/TableSX.core.microbiota.csv", row.names = FALSE)
  
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
# Alpha diversity

{
  alpha_table <- metadata_rare |> 
    rownames_to_column(var = "Sample_ID") |> 
    dplyr::select(
      Sample_ID, sex, morph, month, host_plant, city, green_space, 
      latitude, longitude, 
      observed_ASVs, shannon,simpson
    ) 
  
  write.csv(alpha_table, "tables/TableS4.alpha.diversity.csv", row.names = FALSE)
  
  
  alpha.columns <- c("observed_ASVs", "shannon", "simpson") 
  alpha.labels <- c("Observed ASVs", "Shannon Index", "Simpson Index")
  names(alpha.labels) <- alpha.columns

  {
    theme_set(theme_bw())
    update_theme(legend.position = "none")
    update_theme(axis.ticks = element_blank())
    update_theme(panel.grid.major.x = element_blank())
    update_theme(strip.text = element_text(size = 12, hjust = 0))
    update_theme(strip.background = element_blank())
  }  
  
  # Alpha diversity by sex
  plot_alpha_sex <- metadata_rare |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=sex, y=value, fill = sex)) +
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
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.morph) +
    labs(x = NULL, y = "alpha diversity index", title = "Wing morph") 
  
  plot_alpha_morph
  
  # Alpha diversity by host plant   
  plot_alpha_host <- metadata_rare_july |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=host_plant, y=value, fill = host_plant)) +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.plant) +
    labs(x = NULL, y = "alpha diversity index", title = "Host Plant") 
  
  plot_alpha_host
  
  # Alpha diversity by month
  plot_alpha_month <- metadata_rare_grt |> 
    pivot_longer(cols = all_of(alpha.columns)) |> 
    mutate(name = alpha.labels[match(name, alpha.columns)]) |> 
    ggplot(aes(x=month, y=value, fill = month)) +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.month) +
    labs(x = NULL, y = "alpha diversity index", title = "Month") 
  
  plot_alpha_month
  
  # Alpha diversity by urbanicity
  plot_alpha_urban_Obs <- metadata_rare |> 
    ggplot(aes(x = urbanicity, y = observed_ASVs)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
    geom_point(aes(color = urbanicity), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = "Urbanicity PC1", y = "Observed ASVs", title = "Urbanicity") 
  
  plot_alpha_urban_Obs
  
  plot_alpha_urban_shannon <- metadata_rare |> 
    ggplot(aes(x = urbanicity, y = shannon)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
    geom_point(aes(color = urbanicity), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = "Urbanicity PC1", y = "Simpson's index", title = "Urbanicity") 
  
  plot_alpha_urban_shannon
  
  plot_alpha_urban_SI <- metadata_rare |> 
    ggplot(aes(x = urbanicity, y = simpson)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
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
    theme(panel.grid.major.x = NULL) + # Reset defaults
    geom_point(aes(color = green_space), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    xlim(c(0,0.8)) +
    labs(x = "Local green space (m² / km²)", y = "Observed ASVs", title = "Green space") 
  
  plot_alpha_green_Obs

  plot_alpha_green_shannon <- metadata_rare |> 
    ggplot(aes(x = green_space, y = shannon)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
    geom_point(aes(color = green_space), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    xlim(c(0,0.8)) +
    labs(x = "Local green space (m² / km²)", y = "Shannon's index", title = "Green space") 
  
  plot_alpha_green_shannon  
  
  plot_alpha_green_SI <- metadata_rare |> 
    ggplot(aes(x = green_space, y = simpson)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
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
  # W = 2149.5, p-value = 0.07717
  wilcox.test(shannon ~ sex, data = metadata_rare)
  # W = 2500.5, p-value = 0.7176
  wilcox.test(simpson ~ sex, data = metadata_rare)
  # W = 2593.5, p-value = 0.9952
  
  # morph 
  x <- metadata_rare |> filter(morph %in% c("LW","SW"))
  wilcox.test(observed_ASVs ~ morph, data = x)
  # W = 2074, p-value = 0.4926
  wilcox.test(shannon ~ morph, data = x)
  # W = 2112, p-value = 0.6029
  wilcox.test(simpson ~ morph, data = x)
  # W = 2108, p-value = 0.5909
  
  # host plant
  wilcox.test(observed_ASVs ~ host_plant, data = metadata_rare_july)
  # W = 316, p-value = 0.06921
  wilcox.test(shannon ~ host_plant, data = metadata_rare_july)
  # W = 352, p-value = 0.01435
  wilcox.test(simpson ~ host_plant, data = metadata_rare_july)
  # W = 356, p-value = 0.01177
  
  # month (January vs July) ***
  wilcox.test(observed_ASVs ~ month, data = metadata_rare_grt)
  # W = 1491, p-value = 0.0005031
  wilcox.test(shannon ~ month, data = metadata_rare_grt)
  # W = 1306.5, p-value = 1.957e-05
  wilcox.test(simpson ~ month, data = metadata_rare_grt)
  # W = 1337.5, p-value = 3.538e-05
  
  # LMM for urbanicity
  urban_model <- lmer(observed_ASVs ~ urbanicity + (1 | site_name), data = metadata_rare)
  summary(urban_model)
  # t value Pr(>|t|)    
  # -0.846    0.423    
  urban_model <- lmer(shannon ~ urbanicity + (1 | site_name), data = metadata_rare)
  summary(urban_model)
  # -1.544 0.182308   
  urban_model <- lmer(simpson ~ urbanicity + (1 | site_name), data = metadata_rare)
  summary(urban_model)
  # -1.553 0.174039    
  
  # LMM for green space
  green_space_model <- lmer(observed_ASVs ~ green_space + (1 | site_name), data = metadata_rare)
  summary(green_space_model)
  # 0.534  0.60365  
  green_space_model <- lmer(shannon ~ green_space + (1 | site_name), data = metadata_rare)
  summary(green_space_model)
  # 0.534  0.60365     
  green_space_model <- lmer(simpson ~ green_space + (1 | site_name), data = metadata_rare)
  summary(green_space_model)
  # -1.369 0.198404      
}

##########
# Assembling figures

{
  alpha_figure_NS <- (
    plot_alpha_sex /
    plot_alpha_morph /
    (plot_alpha_urban_Obs | plot_alpha_urban_shannon | plot_alpha_urban_SI) /
    (plot_alpha_green_Obs | plot_alpha_green_shannon | plot_alpha_green_SI)
  ) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect")
  alpha_figure_NS
  
  ggsave(
    filename = "plots/FigS3.alpha.diversity.NS.png", 
    plot = alpha_figure_NS,
    width = 6.5, height = (2.0)*4, scale = 1.5, 
    dpi = 600, bg = "white"
  )
  
  alpha_figure <- (
    plot_alpha_host | plot_alpha_month 
  ) +
    plot_annotation(tag_levels = 'A')
  alpha_figure
  
  ggsave(
    filename = "plots/Fig3.alpha.diversity.png", 
    plot = alpha_figure,
    width = (2.5)*2, height = 2, scale = 2, 
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
  
  ord_rare_bray <- ordinate(ps_rare, method = "PCoA", distance = dist_rare_bray)
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
  
  {
    theme_set(theme_bw())
    update_theme(legend.position = "none")
    update_theme(axis.ticks = element_blank())
  }  
  
  {
    # Host plant / global
    pcoa.global.bray.host_plant <- plot_ordination(ps_rare, ord_rare_bray, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed() +
      labs(title="Bray-Curtis dissimilarity")
    pcoa.global.bray.host_plant
    
    pcoa.global.uwUF.host_plant <- plot_ordination(ps_rare_tree, ord_rare_uwUF, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed() +
      labs(title="unweighted UniFrac")
    pcoa.global.uwUF.host_plant
    
    pcoa.global.wUF.host_plant <- plot_ordination(ps_rare_tree, ord_rare_wUF, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed() +
      labs(title="weighted UniFrac")
    pcoa.global.wUF.host_plant
    
    
    # Month / global
    pcoa.global.bray.month <- plot_ordination(ps_rare, ord_rare_bray, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.global.bray.month
    
    pcoa.global.uwUF.month <- plot_ordination(ps_rare_tree, ord_rare_uwUF, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.global.uwUF.month
    
    pcoa.global.wUF.month <- plot_ordination(ps_rare_tree, ord_rare_wUF, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.global.wUF.month
    
    
    # Host plant / July
    pcoa.july.bray.host_plant <- plot_ordination(ps_july, ord_july_bray, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed()
    pcoa.july.bray.host_plant
    
    pcoa.july.uwUF.host_plant <- plot_ordination(ps_july_tree, ord_july_uwUF, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed()
    pcoa.july.uwUF.host_plant
    
    pcoa.july.wUF.host_plant <- plot_ordination(ps_july_tree, ord_july_wUF, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed()
    pcoa.july.wUF.host_plant
    
    
    # Month / GRT
    pcoa.grt.bray.month <- plot_ordination(ps_grt, ord_grt_bray, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.grt.bray.month
    
    pcoa.grt.uwUF.month <- plot_ordination(ps_grt_tree, ord_grt_uwUF, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.grt.uwUF.month
    
    pcoa.grt.wUF.month <- plot_ordination(ps_grt_tree, ord_grt_wUF, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.grt.wUF.month
   
    # green_space / global
    pcoa.global.bray.green_space <- plot_ordination(ps_rare, ord_rare_bray, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[1], high = color.by.urban[2]) +
      coord_fixed()
    pcoa.global.bray.green_space
    
    pcoa.global.uwUF.green_space <- plot_ordination(ps_rare_tree, ord_rare_uwUF, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[1], high = color.by.urban[2]) +
      coord_fixed()
    pcoa.global.uwUF.green_space
    
    pcoa.global.wUF.green_space <- plot_ordination(ps_rare_tree, ord_rare_wUF, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[1], high = color.by.urban[2]) +
      coord_fixed()
    pcoa.global.wUF.green_space
    
  }
  
  figure_PCoA <- (
    (pcoa.global.bray.host_plant | pcoa.global.uwUF.host_plant | pcoa.global.wUF.host_plant) /
      (pcoa.global.bray.month | pcoa.global.uwUF.month | pcoa.global.wUF.month) /
      (pcoa.july.bray.host_plant | pcoa.july.uwUF.host_plant | pcoa.july.wUF.host_plant) /
      (pcoa.grt.bray.month | pcoa.grt.uwUF.month | pcoa.grt.wUF.month)
  ) +
    plot_annotation(tag_levels = 'A') 
  
  ggsave(
    filename = "plots/FigS4.PCoA.png", 
    plot = figure_PCoA,
    width = (2.5)*3, height = (2.5)*4, scale = 1.25, 
    dpi = 600, bg = "white"
  )
  
}


##########
# PERMANVOA Assumptions
{
  # Testing dispersion of each factor
  dispersion_test <- betadisper(dist_rare_bray, metadata_rare$sex)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_rare_bray, metadata_rare$morph)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_july_bray, metadata_rare_july$host_plant)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_grt_bray, metadata_rare_grt$month)
  permutest(dispersion_test) # NS
  
  dispersion_test <- betadisper(dist_rare_uwUF, metadata_rare$sex)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_rare_uwUF, metadata_rare$morph)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_july_uwUF, metadata_rare_july$host_plant)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_grt_uwUF, metadata_rare_grt$month)
  permutest(dispersion_test) # 0.002 **
  
  # Look in more detail ar dispersion by month
  dispersion_by_month <- betadisper(dist_grt_uwUF, metadata_rare_grt$month)
  dispersion_by_month <- data.frame(
    Group = dispersion_by_month$group,
    Distance = dispersion_by_month$distances
  )
  plot.dispersion_by_month <- ggplot(dispersion_by_month, aes(x = Group, y = Distance)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) + # Shows the actual samples
    labs(y = "Distance to Centroid (Beta Dispersion)") 
  
  plot.dispersion_by_month
  
  dispersion_test <- betadisper(dist_rare_wUF, metadata_rare$sex)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_rare_wUF, metadata_rare$morph)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_july_wUF, metadata_rare_july$host_plant)
  permutest(dispersion_test) # NS
  dispersion_test <- betadisper(dist_grt_wUF, metadata_rare_grt$month)
  permutest(dispersion_test) # NS
}

##########
# PERMANOVA Models 
# Partial sums of squares (the equivalent of Type III ANOVA)

# Full Models
{
  iterations <- 1e3
  perm_bray_full <- adonis2(
    dist_rare_bray ~ sex + morph + host_plant + month + green_space, 
    data = metadata_rare, 
    by = "margin",  
    permutations = iterations
  )
  perm_bray_full
  #              Df SumOfSqs      R2       F   Pr(>F)    
  # sex           1   0.1469 0.00491  0.7919 0.500500    
  # morph         2   0.2861 0.00956  0.7714 0.554446    
  # host_plant    1   2.0693 0.06914 11.1585 0.000999 ***
  # month         1   0.3815 0.01275  2.0574 0.090909 .  
  # green_space   1   0.8789 0.02937  4.7394 0.001998 ** 
  # Residual    137  25.4065 0.84892                     
  # Total       143  29.9280 1.00000   
  
  perm_uwUF_full <- adonis2(
    dist_rare_uwUF ~ sex + morph + host_plant + month + green_space, 
    data = metadata_rare, 
    by = "margin",  
    permutations = iterations
  )
  perm_uwUF_full
  #              Df SumOfSqs      R2      F   Pr(>F)    
  # sex           1    0.412 0.00937 1.4688 0.080919 .  
  # morph         2    0.741 0.01687 1.3215 0.105894    
  # host_plant    1    1.819 0.04138 6.4842 0.000999 ***
  # month         1    1.542 0.03509 5.4979 0.000999 ***
  # green_space   1    0.994 0.02260 3.5418 0.000999 ***
  # Residual    137   38.434 0.87437                    
  # Total       143   43.956 1.00000  
  
  perm_wUF_full <- adonis2(
    dist_rare_wUF ~ sex + morph + host_plant + month + green_space, 
    data = metadata_rare, 
    by = "margin",  
    permutations = iterations
  )
  perm_wUF_full
  #              Df SumOfSqs      R2      F  Pr(>F)    
  # sex           1  0.00017 0.00009 0.0137 0.94805  
  # morph         2  0.01561 0.00870 0.6337 0.56543  
  # host_plant    1  0.02669 0.01488 2.1675 0.10090  
  # month         1  0.00852 0.00475 0.6918 0.46853  
  # green_space   1  0.04483 0.02500 3.6411 0.04296 *
  # Residual    137  1.68680 0.94068                 
  # Total       143  1.79317 1.00000        
}

# One-factor models
{
  perm_bray_host <- adonis2(
    dist_july_bray ~ host_plant, 
    data = metadata_rare_july, 
    by = "margin",  
    permutations = iterations
  )
  perm_bray_host
  #            Df SumOfSqs      R2      F   Pr(>F)    
  # host_plant  1   2.6586 0.14451 14.865 0.000999 ***
  # Residual   88  15.7391 0.85549                    
  # Total      89  18.3977 1.00000   
  
  perm_uwUF_host <- adonis2(
    dist_july_uwUF ~ host_plant, 
    data = metadata_rare_july, 
    by = "margin",  
    permutations = iterations
  )
  perm_uwUF_host
  #            Df SumOfSqs      R2      F   Pr(>F)    
  # host_plant  1   1.7687 0.07362 6.9929 0.000999 ***
  # Residual   88  22.2572 0.92638                    
  # Total      89  24.0258 1.00000    
  
  perm_wUF_host <- adonis2(
    dist_july_wUF ~ host_plant, 
    data = metadata_rare_july, 
    by = "margin",  
    permutations = iterations
  )
  perm_wUF_host
  #            Df SumOfSqs      R2      F  Pr(>F)  
  # host_plant  1  0.03196 0.03447 3.1414 0.04595 *
  # Residual   88  0.89516 0.96553                 
  # Total      89  0.92712 1.00000        

  perm_bray_month <- adonis2(
    dist_grt_bray ~ month, 
    data = metadata_rare_grt, 
    by = "margin",  
    permutations = iterations
  )
  perm_bray_month
  #           Df SumOfSqs      R2      F Pr(>F)
  # month      1   0.4096 0.01536 2.1375 0.1039
  # Residual 137  26.2553 0.98464              
  # Total    138  26.6649 1.00000    
  
  perm_uwUF_month <- adonis2(
    dist_grt_uwUF ~ month, 
    data = metadata_rare_grt, 
    by = "margin",  
    permutations = iterations
  )
  perm_uwUF_month
  #           Df SumOfSqs      R2      F   Pr(>F)    
  # month      1    1.679 0.04059 5.7962 0.000999 ***
  # Residual 137   39.683 0.95941                    
  # Total    138   41.361 1.00000
  
  perm_wUF_month <- adonis2(
    dist_grt_wUF ~ month, 
    data = metadata_rare_grt, 
    by = "margin",  
    permutations = iterations
  )
  perm_wUF_month
  #           Df SumOfSqs      R2      F Pr(>F)
  # month      1   0.0096 0.00548 0.7547 0.4236
  # Residual 137   1.7420 0.99452              
  # Total    138   1.7516 1.00000 
  
}


################################################################################
# Save the image
save.image("Jhae.microbiome.RData")

# load("Jhae.microbiome.RData", verbose = TRUE)

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
# - "MaAsLin3.R"
# - "permdisp.R

