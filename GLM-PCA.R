################################################################################
# Generalized Linear Model Principal Component Analysis (GLM-PCA)
# Dimensionality reduction on raw counts using a Negative Binomial model

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
  library(glmpca)
}

load("Jhae.microbiome.RData", verbose = TRUE)
load("beta.analysis.RData", verbose = TRUE)

{
  # Extract raw counts
  # `glmpca` expects a matrix with ASVs as rows and samples as columns
  count_matrix <- as(otu_table(ps_main), "matrix")
  
  # Ensure the matrix contains integers
  count_matrix <- round(count_matrix)
  
  # Filter out zero-count rows
  if (any(rowSums(count_matrix) == 0)) {
    count_matrix <- count_matrix[-which(rowSums(count_matrix) == 0),]
  }
  
  # Run GLM-PCA
  glmpca_res <- glmpca(
    Y = count_matrix,
    L = 2,           # Extract PC1-2
    fam = "nb"       # "nb" = Negative Binomial, "poi" = Poisson
  )
  
  # Extract coordinates (`factors`)
  glmpca_coords <- metadata_main |>
    mutate(host_plant = str_replace(host_plant,"BV","C. corindum")) |> 
    mutate(host_plant = str_replace(host_plant,"GRT","K. paniculata")) |> 
    rownames_to_column(var = "SampleID") |>
    left_join(
      rownames_to_column(glmpca_res$factors, var = "SampleID"), .,
      by = "SampleID"
    )
  
}

################################################################################
# GLM-PCA Plots

{
  theme_set(theme_bw())
  update_theme(legend.position = "none")
  update_theme(axis.ticks = element_blank())
  update_theme(plot.title = element_text(size = 10))
  update_theme(axis.title = element_text(size = 10))
  update_theme(legend.title = element_text(size = 8))
  update_theme(legend.text = element_text(size = 8))
}

plot_glmpca_host <- glmpca_coords |>
  ggplot(aes(x = dim1, y = dim2, color = host_plant)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(
    aes(fill = host_plant),
    geom = "polygon", alpha = 0.075, show.legend = FALSE
  ) +
  labs(
    title = "GLM-PCA\n(host plant)",
    x = "GLM-PC1", y = "GLM-PC2",
    color = "Host plant"
  ) +
  scale_color_manual(values = color.by.plant) +
  scale_fill_manual(values = color.by.plant) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(face = "italic")
  ) +
  coord_fixed()

plot_glmpca_host

plot_glmpca_month <- glmpca_coords |>
  ggplot(aes(x = dim1, y = dim2, color = month)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  labs(
    title = "GLM-PCA\n(month)",
    x = "GLM-PC1", y = "GLM-PC2",
    color = "Month"
  ) +
  scale_color_manual(values = color.by.month) +
  scale_fill_manual(values = color.by.month) +
  theme(
    legend.position = "bottom"
  ) +
  coord_fixed()

plot_glmpca_month

# tertiles
glmpca_coords$intensive_development_tertile <- cut(
  glmpca_coords$intensive_development, 
  breaks = quantile(glmpca_coords$intensive_development, 
                    probs = seq(0, 1, 1/3), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Low", "Medium", "High")
) 

plot_glmpca_intensive_development <- glmpca_coords |>
  # Note that 10000 (1e5) m^2 is 1 hectare (ha)
  ggplot(aes(x = dim1, y = dim2, color = (intensive_development)/1e5)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_gradient(name = "Intensive Development (hectares)", low = color.by.urban[1], high = color.by.urban[2]) +
  stat_ellipse(aes(group = intensive_development_tertile, fill = intensive_development_tertile), 
               geom = "polygon", alpha = 0.2, color = "gray95") +  
  scale_fill_manual(values = c("Low" = color.by.urban[1], "Medium" = "palegreen3", "High" = color.by.urban[2])) +
  labs(
    title = "GLM-PCA\n(urbanicity)",
    x = "GLM-PC1", y = "GLM-PC2",
    fill = "Intensive Development"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  coord_fixed()

plot_glmpca_intensive_development

################################################################################
# Multi-panel plots

glmpca_figure <- (
  plot_glmpca_host | plot_glmpca_month | plot_glmpca_intensive_development
) +
  plot_annotation(tag_levels = 'A')
# glmpca_figure

# ggsave(
#   filename = "plots/GLM-PCA.plots.png",
#   plot = glmpca_figure,
#   width = 6.5, height = 3, scale = 1.5,
#   dpi = 600, bg = "white"
# )

# A figure combining GLM-PCA results with a few of the PCoA plots
pcoa.july.unifrac.host_plant <- pcoa.july.unifrac.host_plant + labs(title="unweighted UniFrac\n(host plant)")
dpcoa.july.rao.host_plant <- dpcoa.july.rao.host_plant + labs(title="Rao's dissimilarity\n(host plant)")
dpcoa.grt.rao.intensive_development <- dpcoa.grt.rao.intensive_development + labs(title="Rao's dissimilarity\n(urbanicity)")

ordination_figure <- (
  (plot_glmpca_host | plot_glmpca_month | plot_glmpca_intensive_development) /
  (pcoa.july.unifrac.host_plant | dpcoa.july.rao.host_plant | dpcoa.grt.rao.intensive_development)
) +
  plot_annotation(tag_levels = 'A') 
ordination_figure

ggsave(
  filename = "plots/Fig4.ordination.plots.png", 
  plot = ordination_figure,
  width = 6.5, height = 3.75, scale = 1.75,
  dpi = 600#, bg = "white"
)


################################################################################
# Save the key objects
save(
  count_matrix, glmpca_res, glmpca_coords,
  plot_glmpca_host, plot_glmpca_month, plot_glmpca_intensive_development,
  glmpca_figure,
  ordination_figure,
  file = "glmpca.RData"
)
