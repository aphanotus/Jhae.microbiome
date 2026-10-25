################################################################################
# Generalized Linear Model Principal Component Analysis (GLM-PCA)
# Dimensionality reduction on raw counts using a Negative Binomial model

{
  library(tidyverse)
  library(phyloseq)
  library(glmpca)
}

# load("Jhae.microbiome.RData", verbose = TRUE)

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
  rownames_to_column(var = "SampleID") |> 
  left_join(
    rownames_to_column(glmpca_res$factors, var = "SampleID"), ., 
    by = "SampleID"
  )

# Plots
plot_glmpca_host <- glmpca_coords |>
  ggplot(aes(x = dim1, y = dim2, color = host_plant)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(
    aes(fill = host_plant), 
    geom = "polygon", alpha = 0.075, show.legend = FALSE
  ) +
  labs(
    x = "GLM-PC1", y = "GLM-PC2",
    color = "Host plant"
  ) +
  scale_color_manual(values = color.by.plant) +  
  scale_fill_manual(values = color.by.plant) +  
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  coord_fixed()

plot_glmpca_host

plot_glmpca_month <- glmpca_coords |>
  ggplot(aes(x = dim1, y = dim2, color = month)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  labs(
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


glmpca_figure <- (
  plot_glmpca_host | plot_glmpca_month
) +
  plot_annotation(tag_levels = 'A') 
glmpca_figure

ggsave(
  filename = "plots/FigX.glmpca.png", 
  plot = glmpca_figure,
  width = 6.5, height = 3.5, scale = 1.5, 
  dpi = 600, bg = "white"
)


################################################################################
# Save the key objects
save(
  count_matrix, glmpca_res, glmpca_coords, 
  plot_glmpca_host, plot_glmpca_month,
  glmpca_figure,
  file = "Jhae.microbiome.glmpca.RData"
)
