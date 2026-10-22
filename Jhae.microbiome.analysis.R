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
    select(-confidence) |> 
    as.matrix()
  
  TAX_obj <- tax_table(cleaned_tax_matrix)
  
  # Combine the data
  ps <- phyloseq(ASV_obj, META_obj, TAX_obj)
  
  # Subset only adult Jadera haematoloma
  ps <- subset_samples(ps, ((stage == "adult") & (species == "Jhae")))
  
  # Exclude any samples with unspecified sex
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
  
  # Alpha diversity metrics
  richness_df <- estimate_richness(ps, measures = c("Observed","Shannon","Simpson"))
  
  # Add the Observed OTU number 
  sample_data(ps)$observed_OTUs <- richness_df$Observed
  sample_data(ps)$shannon <- richness_df$Shannon
  sample_data(ps)$simpson <- richness_df$Simpson
  
  # Extract the metadata into a standard data frame
  metadata <- as(sample_data(ps), "data.frame")
  
  # # Rarefaction
  # {
  #   # Extract the OTU table
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
  #     xlab = "Sequencing Depth", ylab = "Number of OTUs"
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
  # | NA       | 2148 |    NA | 177     |    NA |
  # | 100      | 2133 | 99.3% | 163     | 92.1% |
  # | 500      | 2101 | 97.8% | 153     | 86.4% |
  # | 1000     | 2064 | 96.0% | 144     | 81.4% |
  # | 2000     | 2055 | 95.7% | 138     | 78.0% |
  
  # Prune out samples with fewer than 1000 total reads
  ps_main <- prune_samples(sample_sums(ps) >= 1000, ps)

  # Clean up any ASVs that now have zero counts because their samples were dropped
  ps_main <- prune_taxa(taxa_sums(ps_main) > 0, ps_main)

  # Extract the metadata into a standard data frame
  metadata_main <- as(sample_data(ps_main), "data.frame")
}

# Read depth check by seasonality
{
  boxplot(log10(read.count) ~ month, data = metadata_main)
  wilcox.test(read.count ~ month, data = metadata_main)
  # NS
  
  boxplot(log10(Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test(Wolbachia.reads ~ month, data = metadata_main)
  # W = 3348.5, p-value = 0.0001481
  
  boxplot(log10(read.count + Wolbachia.reads) ~ month, data = metadata_main)
  wilcox.test((read.count + Wolbachia.reads) ~ month, data = metadata_main)
  # W = 2360, p-value = 0.7743
  # Interpretation: Wolbachia is significantly more abundant in January samples
  # (two-tailed Wilcoxon rank sum test, W = 3348.5, p-value = 0.0001481),
  # although these is not different in total read counts.
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
  ps_bv <- subset_samples(ps_main, host_plant == "BV")
  ps_bv <- prune_taxa(taxa_sums(ps_bv) > 0, ps_bv)
  metadata_bv <- as(sample_data(ps_bv), "data.frame")
  
  ps_july <- subset_samples(ps_main, month == "July")
  ps_july <- prune_taxa(taxa_sums(ps_july) > 0, ps_july)
  metadata_july <- as(sample_data(ps_july), "data.frame")
  
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
}


################################################################################
# Composition Bar Plot

{
  ps_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
  
  # Melt the phyloseq object into a tidy data frame for ggplot
  df_rel <- psmelt(ps_rel)
  
  # Find the top 10 most abundant Families overall
  top_families <- df_rel |>
    group_by(Family) |>
    summarize(Mean_Abund = mean(Abundance)) |>
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
    summarize(Mean_Abundance = mean(Abundance), .groups = 'drop') |> 
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
      strip.background = element_rect(fill = "gray95", color = "gray80")
    )
  
  # plot_composition_grouped_by_host_plant
  
  # Group by month, then calculate the average relative abundance
  plot_composition_grouped_by_month <- df_plot |>
    group_by(month, Family_Clean) |>
    summarize(Mean_Abundance = mean(Abundance), .groups = 'drop') |> 
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
      strip.background = element_rect(fill = "gray95", color = "gray80")
    )
  
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
  ps_bv_rel <- transform_sample_counts(ps_bv, function(x) x / sum(x))
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
    select(-Kingdom)
  
  # Sort the table by highest Mean Abundance (and then by Prevalence)
  bv_summary_sorted <- bv_summary[order(-bv_summary$Mean_Abundance, -bv_summary$Prevalence), ]
  
  # View the top 10 most dominant ASVs in the ancestral diet
  head(bv_summary_sorted, 10)
  
  # Mean Abundance and Prevalence
  ps_grt_july_rel <- transform_sample_counts(ps_grt_july, function(x) x / sum(x))
  mean_abund_grt_july <- taxa_sums(ps_grt_july_rel) / nsamples(ps_grt_july_rel)
  
  ps_grt_july_pa <- transform_sample_counts(ps_grt_july, function(x) ifelse(x > 0, 1, 0))
  
  grt_july_summary <- data.frame(
    ASV_ID = taxa_names(ps_grt_july),
    Mean_Abundance = mean_abund_grt_july,
    Prevalence = taxa_sums(ps_grt_july_pa) / nsamples(ps_grt_july_pa),
    as(tax_table(ps_grt_july), "matrix")
  ) |> 
    select(-Kingdom)
  
  grt_july_summary_sorted <- grt_july_summary[order(-grt_july_summary$Mean_Abundance, -grt_july_summary$Prevalence), ]
  
  head(grt_july_summary_sorted, 10)
}
##########
# Host shift 
{
  bv_summary_sorted$host <- "Cardiospermum"
  grt_july_summary_sorted$host <- "Koelreuteria"
  
  combined_july_summary <- bind_rows(bv_summary_sorted, grt_july_summary_sorted)
  
  # Extract the Top 10 from BV and July GRT
  target_asvs_july <- unique(c(
    head(bv_summary_sorted$ASV_ID, 5), 
    head(grt_july_summary_sorted$ASV_ID, 2)
  ))
  
  uuid.abbrev <- 8
  
  slope_df_july <- combined_july_summary |>
    filter(ASV_ID %in% target_asvs_july) |>
    mutate(
      Host = factor(host, levels = c("Cardiospermum", "Koelreuteria")),
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
      data = slope_df_july |> filter(Host == "Cardiospermum"),
      aes(label = Taxa_Label),
      nudge_x = -0.15, hjust = 1, direction = "y", fontface = "italic", size = 4
    ) +
    geom_text_repel(
      data = slope_df_july |> filter(Host == "Koelreuteria"),
      aes(label = Taxa_Label),
      nudge_x = 0.15, hjust = 0, direction = "y", fontface = "italic", size = 4
    ) +
    scale_y_continuous(
      name = "Mean Relative Abundance",
      labels = scales::percent_format(), limits = c(0, NA)
    ) +
    scale_color_brewer(palette = "Dark2") + 
    scale_x_discrete(name = NULL, expand = expansion(mult = c(1, 1))) +
    # scale_x_discrete(name = NULL) + 
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "italic", size = 12)
    )
  
  plot_host_shift
  
  ggsave(
    filename = "plots/FigX.core_microbiome_shift.png", 
    plot = plot_host_shift,
    width = 4, height = 3, scale = 1.5, 
    dpi = 600, bg = "white"
  )
  
}


################################################################################
# Alpha diversity

{
  alpha.columns <- c("observed_OTUs", "simpson") # "shannon", 
  alpha.labels <- c("Observed OTUs", "Simpson Index") # "Shannon Index", 
  names(alpha.labels) <- alpha.columns
  
  # Alpha diversity by sex
  plot_alpha_sex <- metadata_main |> 
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
  plot_alpha_morph <- metadata_main |> 
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
  plot_alpha_host <- metadata_july |> 
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
  plot_alpha_month <- metadata_grt |> 
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
  plot_alpha_urban <- metadata_main |> 
    ggplot(aes(x = urbanicity, y = simpson)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = urbanicity), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = "Urbanicity PC1", y = "Simpson's index", title = "Urbanicity") 
  
  plot_alpha_urban
  
  # Alpha diversity by green space
  plot_alpha_green <- metadata_main |> 
    ggplot(aes(x = green_space, y = observed_OTUs)) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_point(aes(color = green_space), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    xlim(c(0,0.8)) +
    labs(x = "Local green space (m² / km²)", y = "Observed OTUs", title = "Green space") 
  
  plot_alpha_green
}

##########
# Alpha diversity stats
{
  # Using Wilcoxon rank sum test
  
  # sex 
  wilcox.test(observed_OTUs ~ sex, data = metadata_main)
  # W = 1921, p-value = 0.007405
  wilcox.test(simpson ~ sex, data = metadata_main)
  # NS
  
  # morph 
  x <- metadata_main |> filter(morph %in% c("LW","SW"))
  wilcox.test(observed_OTUs ~ morph, data = x)
  # NS
  wilcox.test(simpson ~ morph, data = x)
  # NS
  
  # host plant
  wilcox.test(observed_OTUs ~ host_plant, data = metadata_july)
  # NS
  wilcox.test(simpson ~ host_plant, data = metadata_july)
  # W = 359, p-value = 0.01012
  
  # month (January vs July) ***
  wilcox.test(observed_OTUs ~ month, data = metadata_grt)
  # W = 1442, p-value = 0.0002286
  wilcox.test(simpson ~ month, data = metadata_grt)
  # W = 1316, p-value = 2.353e-05
  
  # LMM for urbanicity
  urban_model <- lmer(observed_OTUs ~ urbanicity + (1 | site_name), data = metadata_main)
  summary(urban_model)
  # NS
  urban_model <- lmer(simpson ~ urbanicity + (1 | site_name), data = metadata_main)
  summary(urban_model)
  # NS
  
  # LMM for green space
  green_space_model <- lmer(observed_OTUs ~ green_space + (1 | site_name), data = metadata_main)
  summary(green_space_model)
  # NS
  green_space_model <- lmer(simpson ~ green_space + (1 | site_name), data = metadata_main)
  summary(green_space_model)
  # NS  
}


##########
# Assembling figures

{
  alpha_figure_NS <- (
    plot_alpha_morph | plot_alpha_urban | plot_alpha_green
  ) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect")
  alpha_figure_NS
  
  ggsave(
    filename = "plots/FigX.alpha.diversity.NS.png", 
    plot = alpha_figure_NS,
    width = (2.25)*3, height = 2.5, scale = 1.8, 
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
  dist_main_bray <- distance(ps_main, method = "bray")
  dist_july_bray <- distance(ps_july, method = "bray")
  dist_grt_bray <- distance(ps_july, method = "bray")
  
  ord_july_bray <- ordinate(ps_july, method = "PCoA", distance = dist_july_bray)
  ord_grt_bray <- ordinate(ps_july, method = "PCoA", distance = dist_grt_bray)
  
  # UniFrac distances
  asv.tree <- phytools::read.newick("qiime/16s.Jhae.rooted-tree.tre")
  ps_main_tree <- merge_phyloseq(ps_main, asv.tree)
  ps_july_tree <- merge_phyloseq(ps_july, asv.tree)
  ps_grt_tree <- merge_phyloseq(ps_grt, asv.tree)
  
  dist_main_uwUF <- distance(ps_main_tree, method = "unifrac")
  dist_july_uwUF <- distance(ps_july_tree, method = "unifrac")
  dist_grt_uwUF <- distance(ps_grt_tree, method = "unifrac")
  
  dist_main_wUF <- distance(ps_main_tree, method = "wunifrac")
  dist_july_wUF <- distance(ps_july_tree, method = "wunifrac")
  dist_grt_wUF <- distance(ps_grt_tree, method = "wunifrac")
  
  ord_main_uwUF <- ordinate(ps_main_tree, method = "PCoA", distance = dist_main_uwUF)
  ord_july_uwUF <- ordinate(ps_july_tree, method = "PCoA", distance = dist_july_uwUF)
  ord_grt_uwUF <- ordinate(ps_grt_tree, method = "PCoA", distance = dist_grt_uwUF)
  
  ord_main_wUF <- ordinate(ps_main_tree, method = "PCoA", distance = dist_main_wUF)
  ord_july_wUF <- ordinate(ps_july_tree, method = "PCoA", distance = dist_july_wUF)
  ord_grt_wUF <- ordinate(ps_grt_tree, method = "PCoA", distance = dist_grt_wUF)
}

# PCoA Plots
pcoa.bray.host <- plot_ordination(ps_july, ord_july_bray, color = "host_plant") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(name = "Host plant", values = color.by.plant) +
  coord_fixed()
pcoa.bray.host

pcoa.uwUF.host <- plot_ordination(ps_july_tree, ord_july_uwUF, color = "host_plant") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color.by.plant) +
  coord_fixed()
pcoa.uwUF.host

pcoa.wUF.host <- plot_ordination(ps_july_tree, ord_july_wUF, color = "host_plant") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color.by.plant) +
  coord_fixed()
pcoa.wUF.host


##########
# PERMANVOA Assumptions

# Testing dispersion of each factor
dispersion_test <- betadisper(dist_main_bray, metadata_main$sex)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_main_bray, metadata_main$morph)
permutest(dispersion_test) # 0.008 **
dispersion_test <- betadisper(dist_main_bray, metadata_main$host_plant)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_main_bray, metadata_main$month)
permutest(dispersion_test) # NS

dispersion_test <- betadisper(dist_main_uwUF, metadata_main$sex)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_main_uwUF, metadata_main$morph)
permutest(dispersion_test) # 0.001 ***
dispersion_test <- betadisper(dist_main_uwUF, metadata_main$host_plant)
permutest(dispersion_test) # 0.001 ***
dispersion_test <- betadisper(dist_main_uwUF, metadata_main$month)
permutest(dispersion_test) # NS

# Look in more detail ar dispersion by host_plant
dispersion_by_host_plant <- betadisper(dist_main_uwUF, metadata_main$host_plant)
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

dispersion_test <- betadisper(dist_main_wUF, metadata_main$morph)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_main_wUF, metadata_main$morph)
permutest(dispersion_test) # NS
dispersion_test <- betadisper(dist_main_wUF, metadata_main$host_plant)
permutest(dispersion_test) # 0.001 ***
dispersion_test <- betadisper(dist_main_wUF, metadata_main$month)
permutest(dispersion_test) # 0.041 *

# Look in more detail ar dispersion by host_plant
dispersion_by_host_plant <- betadisper(dist_main_wUF, metadata_main$host_plant)
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
  dist_main_bray ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_main, 
  by = "margin",  
  permutations = iterations
)

perm_bray_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1    0.362 0.00821 1.2484 0.220779    
# morph        2    0.599 0.01358 1.0324 0.373626    
# host_plant   1    1.965 0.04457 6.7744 0.000999 ***
# month        1    0.953 0.02161 3.2851 0.008991 ** 
# urbanicity   1    0.562 0.01276 1.9393 0.063936 .  
# Residual   137   39.729 0.90125                    
# Total      143   44.083 1.00000     

perm_uwUF_full <- adonis2(
  dist_main_uwUF ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_main, 
  by = "margin",  
  permutations = iterations
)

perm_uwUF_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1    0.550 0.00863 1.2888 0.041958 *  
# morph        2    0.917 0.01439 1.0747 0.224775    
# host_plant   1    1.809 0.02838 4.2381 0.000999 ***
# month        1    0.929 0.01458 2.1772 0.000999 ***
# urbanicity   1    0.821 0.01289 1.9246 0.000999 ***
# Residual   137   58.463 0.91752                    
# Total      143   63.719 1.00000       

perm_wUF_full <- adonis2(
  dist_main_wUF ~ sex + morph + host_plant + month + urbanicity, 
  data = metadata_main, 
  by = "margin",  
  permutations = iterations
)

perm_wUF_full
#             Df SumOfSqs      R2      F   Pr(>F)    
# sex          1 0.0001189 0.00410  0.9731 0.450549    
# morph        2 0.0000802 0.00276  0.3280 0.607393    
# host_plant   1 0.0094116 0.32446 77.0139 0.000999 ***
# month        1 0.0001492 0.00515  1.2213 0.384615    
# urbanicity   1 0.0000962 0.00332  0.7872 0.446553    
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
    Sample_size_in_January = sum(month == "January"),
    Sample_size_in_July = sum(month == "July"),
    Proportion_of_green_space = round(unique(green_space),4),
    Urbanicity_index = round(unique(urbanicity),3)
  ) |> 
  arrange(Urbanicity_index)

write_delim(x, "collection.sites.tsv", delim = "\t")


################################################################################
# Follow-up scripts
# - "ancom.R"
# - "aldex.R"
# - "GLM-PCA.R"
# - "fastEmu.R"
# - "MaAsLin3.R"
# - "Map.of.field.sites.R"

