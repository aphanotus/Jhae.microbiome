################################################################################
# Differential Abundance Testing with fastEmu

# devtools::install_github("statdivlab/fastEmu")
{
  library(tidyverse)
  library(phyloseq)
  library(fastEmu)
}

load("Jhae.microbiome.RData", verbose = TRUE)

################################################################################
# Examining the effect of host_plant
Y <- as.matrix(if(taxa_are_rows(ps_july)) t(otu_table(ps_july)) else otu_table(ps_july))
X <- model.matrix(~ host_plant, data = metadata_july)
k_index <- grep("host_plant", colnames(X))

test_kj <- expand.grid(
  k = k_index,
  j = 1:ncol(Y)
)

j <- 1:ncol(Y)

# NOTE: Time difference of 4.026951 hours
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

# View significant results
x <- fastEmu_host_results |>
  filter(p_adj < 0.05) |>
  dplyr::select(ASV_ID, estimate, pval, p_adj, Phylum, Class, Order, Family, Genus, Species)
x

write.csv(x, "tables/TableSX.fastEmu.host.csv", row.names = FALSE)

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


################################################################################
# Examining the effect of month

metadata_grt <- data.frame(sample_data(ps_grt))
Y <- as.matrix(if(taxa_are_rows(ps_grt)) t(otu_table(ps_grt)) else otu_table(ps_grt))
X <- model.matrix(~ month, data = metadata_grt)
k_index <- grep("month", colnames(X))

test_kj <- expand.grid(
  k = k_index,
  j = 1:ncol(Y)
)

j <- 1:ncol(Y)

# NOTE: Time difference of 15.19104 hours
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

# View significant results
x <- fastEmu_month_results |>
  filter(p_adj < 0.05) |>
  dplyr::select(ASV_ID, estimate, pval, p_adj, Phylum, Class, Order, Family, Genus, Species)
x
  
write.csv(x, "tables/TableSX.fastEmu.month.csv", row.names = FALSE)
  
  
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

plot_fastEmu_volcano_month <- fastEmu_month_results_for_plot |> 
  ggplot(aes(x = estimate, y = log10_pval, color = Significance)) +
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


################################################################################
# Examining the effect of green space

Y <- as.matrix(if(taxa_are_rows(ps_main)) t(otu_table(ps_main)) else otu_table(ps_main))
X <- model.matrix(~ green_space, data = metadata_main)
k_index <- grep("green_space", colnames(X))

test_kj <- expand.grid(
  k = k_index,
  j = 1:ncol(Y)
)

j <- 1:ncol(Y)

# NOTE: Time difference of MANY hours
{
  start.time <- Sys.time()
  fastEmu_fit_green_space <- fastEmuFit(
    formula = ~ green_space,
    data = metadata_main,
    Y = Y,
    test_kj = test_kj,
    reference_set = "data_driven"
  )
  print(Sys.time() - start.time)
}

taxa_info <- as.data.frame(as(tax_table(ps_main), "matrix")) |>
  rownames_to_column(var = "ASV_ID")

fastEmu_green_space_results <- fastEmu_fit_green_space$coef |>
  mutate(
    # Re-map the ASV IDs
    ASV_ID = colnames(Y)[j],
  ) |>
  left_join(taxa_info, by = "ASV_ID") |>
  # Apply Benjamini-Hochberg FDR correction
  mutate(p_adj = p.adjust(pval, method = "BH")) |>
  arrange(p_adj)

# View significant results
x <- fastEmu_green_space_results |>
  filter(p_adj < 0.05) |>
  dplyr::select(ASV_ID, estimate, pval, p_adj, Phylum, Class, Order, Family, Genus, Species)
x

write.csv(x, "tables/TableSX.fastEmu.green_space.csv", row.names = FALSE)


# Prepare results for a volcano plot
fastEmu_green_space_results_for_plot <- fastEmu_green_space_results |>
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

plot_fastEmu_volcano_green_space <- fastEmu_green_space_results_for_plot |> 
  ggplot(aes(x = estimate, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      "Not Significant" = "grey75",
      "Significant" = color.by.urban[2]
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

plot_fastEmu_volcano_green_space


##########
# Summary figure

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
# Summary figure for robust DA ASVs using the fastEmu results

# Manually curated list of robust DA ASvs (those DA by 2 or more methods)
DA.ASVs.by.host <- c(
  "221659bea0a1fb86d90ab90a24cd7dfe", "3cc27a22b5adf4e21f48c1f09ee198de", "5f9403f34de2d6658df09b913aa13184", 
  "78849fe4c1e26973c2a70cf0727da7e8", "90dc413777fc22c4efdda035d0b2edd0", "feeb651ed135a05fff03c033989caded" 
)
DA.ASVs.by.month <- c(
  "08bb6389255945886b20e36acbab566f", "4bca6d2c2632165857073fbaf8e36f2c", "8084aa68823c09679d87c31776e22478", 
  "aa5738a607473674763d7e712002a471", "fcb2943f9c7b866112a146c62e824bdb"
)

DA_host_plot_df <- fastEmu_host_results |>
  mutate(
    log10_pval = -log10(pval),
    Significance = ifelse(ASV_ID %in% DA.ASVs.by.host, "Enriched in C. corindum", "Not Significant"),
    Taxa_Label = case_when(
      !(ASV_ID %in% DA.ASVs.by.host) ~ "",
      !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Phylum) & Phylum != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      TRUE ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_DA_host <- DA_host_plot_df |>
  ggplot(aes(x = estimate, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Not Significant" = "grey75",
      "Enriched in C. corindum" = color.by.plant[1]
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

plot_DA_host



DA_month_plot_df <- fastEmu_month_results |>
  mutate(
    log10_pval = -log10(pval),
    Significance = ifelse(ASV_ID %in% DA.ASVs.by.month, "Enriched in January", "Not Significant"),
    Taxa_Label = case_when(
      !(ASV_ID %in% DA.ASVs.by.month) ~ "",
      !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      !is.na(Phylum) & Phylum != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      TRUE ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  )

plot_DA_month <- DA_month_plot_df |>
  ggplot(aes(x = estimate, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Not Significant" = "grey75",
      "Enriched in January" = color.by.month[1]
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

plot_DA_month


DA_figure <- (
  plot_DA_host | plot_DA_month
) +
  plot_annotation(tag_levels = 'A')
DA_figure

ggsave(
  filename = "plots/FigX.consensus.DA.png",
  plot = DA_figure,
  width = 6.5, height = 3.0, scale = 1.65,
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
  fastEmu_month_results_for_plot, plot_fastEmu_volcano_month, # plot_fastEmu_volcano_green_space,
  fastEmu_figure,
  plot_DA_host, plot_DA_month, DA_figure,
  file = "Jhae.microbiome.fastEmu.RData"
)
