################################################################################
# Differential Abundance Testing with ALDEx2

{
  library(tidyverse)
  library(phyloseq)
  library(ALDEx2)
}

load("Jhae.microbiome.RData", verbose = TRUE)

################################################################################
# Host plant
{
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
    mutate(contrast = "host_plant") |> 
    dplyr::select(ASV_ID, contrast, diff.btw, effect, wi.ep, wi.eBH, Phylum, Class, Order, Family, Genus, Species)
  
  aldex_host_sig  
  
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
    expand_limits(x = c(-6, 4), y = c(0, 6.25)) +
    labs(
      x = "Effect Size",
      y = expression(-log[10]("p-value"))
    )
  
  plot_aldex_host
}

################################################################################
# Month
{
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
    mutate(contrast = "month") |> 
    dplyr::select(ASV_ID, contrast, diff.btw, effect, wi.ep, wi.eBH, Phylum, Class, Order, Family, Genus, Species)
  
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
    expand_limits(x = c(-6, 4), y = c(0, 6.25)) +
    labs(
      x = "Effect Size",
      y = expression(-log[10]("p-value"))
    )
  
  plot_aldex_month
}

################################################################################
# Green Space
{
  # Standardize the green_space variable (Calculate Z-scores)
  # ((value - mean) / standard deviation)
  sample_data(ps_main)$green_space_z <- scale(sample_data(ps_main)$green_space)
  metadata_main <- as(sample_data(ps_main), "data.frame")
  
  counts_main <- round(as(otu_table(ps_main), "matrix"))
  
  # Build the model matrix 
  model_matrix <- model.matrix(~ green_space_z, data = metadata_main)
  
  taxa_main <- as.data.frame(as(tax_table(ps_main), "matrix"))
  
  # ALDEx2 CLR transformation
  set.seed(123) 
  aldex_clr <- aldex.clr(
    reads = counts_main, 
    conds = model_matrix, 
    mc.samples = 128, 
    denom = "all"
  )
  
  # Fits the continuous model across all Monte Carlo instances
  aldex_glm_results <- aldex.glm(aldex_clr, model_matrix)
  
  aldex_green_sig <- aldex_glm_results |>
    rownames_to_column(var = "ASV_ID") |>
    left_join(taxa_main |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
    # Rename the specific GLM columns to standardized names
    rename(
      LogFoldChange = `green_space_z:Est`,
      p_val = `green_space_z:pval`,
      q_val_FDR = `green_space_z:pval.padj`
    ) 
  
  any(aldex_green_sig$`green_space_z:pval.padj` <= 0.05)
  # FALSE
  # No significant ASVs!

  sum(aldex_green_sig$`green_space_z:pval` <= 0.05)
  # 5
  
  # Try this at the Family level
  sample_data(ps_family)$green_space_z <- scale(sample_data(ps_family)$green_space)
  metadata_family <- as(sample_data(ps_family), "data.frame")
  
  counts_family <- round(as(otu_table(ps_family), "matrix"))
  
  # Build the model matrix 
  model_matrix <- model.matrix(~ green_space_z, data = metadata_family)
  
  taxa_family <- as.data.frame(as(tax_table(ps_family), "matrix"))
  
  # ALDEx2 CLR transformation
  set.seed(123) 
  aldex_clr <- aldex.clr(
    reads = counts_family, 
    conds = model_matrix, 
    mc.samples = 128, 
    denom = "all"
  )
  
  # Fits the continuous model across all Monte Carlo instances
  aldex_glm_results <- aldex.glm(aldex_clr, model_matrix)
  
  aldex_green_sig <- aldex_glm_results |>
    rownames_to_column(var = "ASV_ID") |>
    left_join(taxa_main |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
    # Rename the specific GLM columns to standardized names
    rename(
      LogFoldChange = `green_space_z:Est`,
      p_val = `green_space_z:pval`,
      q_val_FDR = `green_space_z:pval.padj`
    ) 
  
  any(aldex_green_sig$q_val_FDR <= 0.05)
  # FALSE
  # No significant ASVs!
  
  # I also tried looking at just the GRT dataset, 
  # and there were no significant results.
}

################################################################################
# ALDEx Volcano Plot Figure

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

# Since green space is a significant factor in PERMANOVAs but 
# no individual ASV (or Family-level taxa) emerge as differentially abundant
# in ALDEx models, perhaps it's differences in variance, not means that
# drive the globsal community contrasts in PERMANOVA?
{
  # Categorize the continuous green_space variable into Tertiles (Low, Med, High)
  metadata_rare_grt <- metadata_rare_grt %>%
    mutate(
      green_space_cat = ntile(green_space, 3),
      green_space_cat = factor(
        green_space_cat, 
        levels = c(1, 2, 3), 
        labels = c("Low (Urban)", "Medium", "High (Rural)"))
    )
  
  dispersion_bray <- betadisper(dist_grt_bray, metadata_rare_grt$green_space_cat)
  permutest_bray <- permutest(dispersion_bray, permutations = 1e4-1)
  print(permutest_bray)
  #            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
  # Groups      2  1.4526 0.72629 7.3334   9999  7e-04 ***
  # Residuals 136 13.4691 0.09904  
  
  dispersion_uwUF <- betadisper(dist_grt_uwUF, metadata_rare_grt$green_space_cat)
  permutest_uwUF <- permutest(dispersion_uwUF, permutations = 1e4-1)
  print(permutest_uwUF)
  #            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
  # Groups      2 0.4246 0.212289 7.086   9999  0.001 ***
  # Residuals 136 4.0744 0.029959   
  
  dispersion_wUF <- betadisper(dist_grt_wUF, metadata_rare_grt$green_space_cat)
  permutest_wUF <- permutest(dispersion_wUF, permutations = 1e4-1)
  print(permutest_uwUF)
  #            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
  # Groups      2 0.4246 0.212289 7.086   9999  0.001 ***
  # Residuals 136 4.0744 0.029959   
  
  {
    theme_set(theme_bw())
    update_theme(legend.position = "none")
    update_theme(axis.ticks = element_blank())
    update_theme(panel.grid.major.x = element_blank())
  }  
  
  # Extract the distance-to-centroid values to see the direction of the variance
  dispersion_bray_for_plot <- data.frame(
    sample = names(dispersion_bray$distances),
    distance_to_centroid = dispersion_bray$distances,
    green_space = metadata_rare_grt$green_space_cat
  )
  
  # Plot the differences in variance of Beta Dispersion by Green Space
  plot_dispersion_bray_by_green_space <- dispersion_bray_for_plot |> 
    ggplot(aes(x = green_space, y = distance_to_centroid, fill = green_space)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
    labs(
      title = "Brasy-Curtis Distance",
      x = "Local Green Space",
      y = "Distance to Group Centroid (Community Variance)"
    ) 
  
  plot_dispersion_bray_by_green_space
  
  dispersion_uwUF_for_plot <- data.frame(
    sample = names(dispersion_uwUF$distances),
    distance_to_centroid = dispersion_uwUF$distances,
    green_space = metadata_rare_grt$green_space_cat
  )
  
  plot_dispersion_uwUF_by_green_space <- dispersion_uwUF_for_plot |> 
    ggplot(aes(x = green_space, y = distance_to_centroid, fill = green_space)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
    labs(
      title = "unweighted UniFrac",
      x = "Local Green Space",
      y = "Distance to Group Centroid (Community Variance)"
    ) 
  
  plot_dispersion_uwUF_by_green_space
  
  dispersion_wUF_for_plot <- data.frame(
    sample = names(dispersion_wUF$distances),
    distance_to_centroid = dispersion_wUF$distances,
    green_space = metadata_rare_grt$green_space_cat
  )
  
  plot_dispersion_wUF_by_green_space <- dispersion_wUF_for_plot |> 
    ggplot(aes(x = green_space, y = distance_to_centroid, fill = green_space)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
    labs(
      title = "weighted UniFrac",
      x = "Local Green Space",
      y = "Distance to Group Centroid (Community Variance)"
    ) 
  
  plot_dispersion_wUF_by_green_space
  
}

beta_dispersion_figure <- (
  plot_dispersion_bray_by_green_space | plot_dispersion_uwUF_by_green_space | plot_dispersion_wUF_by_green_space
) +
  plot_annotation(tag_levels = 'A') 
beta_dispersion_figure

ggsave(
  filename = "plots/FigX.beta.dispersion.by.green.space.png", 
  plot = beta_dispersion_figure,
  width = 6.5, height = 3.5, scale = 1.65, 
  dpi = 600, bg = "white"
)


################################################################################
# Save the key objects
save(
  counts_july, conditions_july, aldex_host, taxa_july, aldex_host_sig, aldex_host_plot_df, plot_aldex_host,
  counts_grt, conditions_grt, aldex_month, taxa_grt, aldex_month_sig, aldex_month_plot_df, plot_aldex_month,
  aldex_figure, 
  ps_main, 
  beta_dispersion_figure,
  file = "Jhae.microbiome.aldex.RData"
)
