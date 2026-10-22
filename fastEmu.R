################################################################################
# Differential Abundance Testing with fastEmu

# devtools::install_github("statdivlab/fastEmu")
{
  library(tidyverse)
  library(phyloseq)
  library(fastEmu)
}

load("Jhae.microbiome.RData", verbose = TRUE)

# Examining the effect of host_plant
Y <- as.matrix(if(taxa_are_rows(ps_july)) t(otu_table(ps_july)) else otu_table(ps_july))
X <- model.matrix(~ host_plant, data = metadata_july)
k_index <- grep("host_plant", colnames(X))

test_kj <- expand.grid(
  k = k_index, 
  j = 1:ncol(Y)
)

j <- 1:ncol(Y)

# NOTE: Time difference of 3.938068 hours
{
  start.time <- Sys.time()
  fastEmu_fit_host_plant <- fastEmuFit(
    formula = ~ host_plant, 
    data = metadata_july, 
    Y = Y, 
    test_kj = test_kj,
    reference_set = "data_driven" 
  )
  print(Sys.time() - start.time)
}

taxa_info <- as.data.frame(as(tax_table(ps_july), "matrix")) |>
  rownames_to_column(var = "ASV_ID")

fastEmu_host_results <- fastEmu_fit_host_plant$coef |>
  mutate(
    # Re-map the ASV IDs
    ASV_ID = colnames(Y)[j],
  ) |>
  left_join(taxa_info, by = "ASV_ID") |>
  # Apply Benjamini-Hochberg FDR correction
  mutate(p_adj = p.adjust(pval, method = "BH")) |>
  arrange(p_adj)

# View the top hits
fastEmu_host_results |> 
  filter(p_adj < 0.05) |> 
  dplyr::select(ASV_ID, estimate, pval, p_adj, Family, Genus) |> 
  slice_head(n=10)

# Plot estimates in the style of the fastEmu vignette
fastEmu_host_results |> 
  arrange(estimate) |> 
  mutate(ord = 1:n()) |>
  ggplot(aes(x = estimate, xmin = lower, xmax = upper, color = Phylum, y = ord)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  scale_colour_viridis(discrete = TRUE) +
  labs(y = "", x = "Estimate", color = "Phylum") + 
  theme_bw() + 
  ggtitle("Estimates and intervals for all mOTUs") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Prepare results for a volcano plot
fastEmu_host_results_for_plot <- fastEmu_host_results |>
  mutate(
    log10_pval = -log10(pval),
    Significance = ifelse(p_adj < 0.05, "Significant", "Not Significant"),
    Taxa_Label = case_when(
      p_adj >= 0.05 ~ "",   
      !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      TRUE ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_fastEmu_volcano_host <- fastEmu_host_results_for_plot |> 
  ggplot(aes(x = estimate, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      "Not Significant" = "grey75", 
      "Significant" = color.by.plant[1]
    )) +
  geom_text_repel(
    aes(label = Taxa_Label),
    size = 3.5, color = "black",
    box.padding = 0.5, point.padding = 0.3, 
    force = 4, max.overlaps = Inf, 
    show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  expand_limits(x = c(-16, 15), y = c(0, 5.5)) +
  labs(
    x = "Log-Fold Change (fastEmu Estimate)",
    y = expression(-log[10]("p-value"))
  ) 

plot_fastEmu_volcano_host



# Examining the effect of month
metadata_grt <- data.frame(sample_data(ps_grt))
Y <- as.matrix(if(taxa_are_rows(ps_grt)) t(otu_table(ps_grt)) else otu_table(ps_grt))
# Filter out ASV columns with zero reads  
Y <- Y[,-which(colSums(Y) == 0)]

X <- model.matrix(~ month, data = metadata_grt)
k_index <- grep("month", colnames(X))

test_kj <- expand.grid(
  k = k_index, 
  j = 1:ncol(Y)
)

j <- 1:ncol(Y)

# Fitting the fastEmu model should take about 2-3 hours on my laptop
{
  start.time <- Sys.time()
  fastEmu_fit_month <- fastEmuFit(
    formula = ~ month, 
    data = metadata_grt, 
    Y = Y, 
    test_kj = test_kj,
    reference_set = "data_driven" 
  )
  print(Sys.time() - start.time)
}

taxa_info <- as.data.frame(as(tax_table(ps_grt), "matrix")) |>
  rownames_to_column(var = "ASV_ID")

fastEmu_month_results <- fastEmu_fit_month$coef |>
  mutate(
    # Re-map the ASV IDs
    ASV_ID = colnames(Y)[j],
  ) |>
  left_join(taxa_info, by = "ASV_ID") |>
  # Apply Benjamini-Hochberg FDR correction
  mutate(p_adj = p.adjust(pval, method = "BH")) |>
  arrange(p_adj)

# View the top hits
fastEmu_month_results |> 
  filter(p_adj < 0.05) |> 
  dplyr::select(ASV_ID, estimate, pval, p_adj, Family, Genus) |> 
  slice_head(n=10)

# Plot estimates in the style of the fastEmu vignette
fastEmu_month_results |> 
  arrange(estimate) |> 
  mutate(ord = 1:n()) |>
  ggplot(aes(x = estimate, xmin = lower, xmax = upper, color = Phylum, y = ord)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  scale_colour_viridis(discrete = TRUE) +
  labs(y = "", x = "Estimate", color = "Phylum") + 
  theme_bw() + 
  ggtitle("Estimates and intervals for all OTUs") + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Prepare results for a volcano plot
fastEmu_month_results_for_plot <- fastEmu_month_results |>
  mutate(
    log10_pval = -log10(pval),
    Significance = ifelse(p_adj < 0.05, "Significant", "Not Significant"),
    Taxa_Label = case_when(
      p_adj >= 0.05 ~ "",   
      !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Phylum) & Class != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      TRUE ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_fastEmu_volcano_month <- ggplot(plot_data_month, aes(x = estimate, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      "Not Significant" = "grey75", 
      "Significant" = color.by.month[2]
    )) +
  geom_text_repel(
    aes(label = Taxa_Label),
    size = 3.5, color = "black",
    box.padding = 0.5, point.padding = 0.3, 
    force = 4, max.overlaps = Inf, 
    show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  expand_limits(x = c(-16, 15), y = c(0, 5.5)) +
  labs(
    x = "Log-Fold Change (fastEmu Estimate)",
    y = expression(-log[10]("p-value"))
  ) 

plot_fastEmu_volcano_month

fastEmu_figure <- (
  plot_fastEmu_volcano_host | plot_fastEmu_volcano_month
) +
  plot_annotation(tag_levels = 'A') 
fastEmu_figure

ggsave(
  filename = "plots/FigX.fastEmu.volcanos.png", 
  plot = fastEmu_figure,
  width = 6.5, height = 3.5, scale = 1.5, 
  dpi = 600, bg = "white"
)


################################################################################
# Save the key objects
save(
  ps_july, metadata_july, fastEmu_host_results, 
  fastEmu_host_results_for_plot, plot_fastEmu_volcano_host,
  ps_grt, metadata_grt, fastEmu_month_results,
  # # the `fastEmuFit` objects themselves are too large for storage on GitHub
  # fastEmu_fit_host_plant, fastEmu_fit_month,  
  fastEmu_month_results_for_plot, plot_fastEmu_volcano_month,
  fastEmu_figure,
  file = "Jhae.microbiome.fastEmu.RData"
)
