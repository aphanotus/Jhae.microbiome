################################################################################
# Generalized Linear Model Principal Component Analysis (GLM-PCA)

{
  library(tidyverse)
  library(phyloseq)
  library(glmpca)
}

# load("Jhae.microbiome.RData", verbose = TRUE)

# For host_plant (in July)

# Extract raw counts
# `glmpca` expects a matrix with ASVs as rows and samples as columns
count_matrix_july <- as(otu_table(ps_july), "matrix")

# Ensure the matrix contains integers
count_matrix_july <- round(count_matrix_july)

# Filter out zero-count rows
if (any(rowSums(count_matrix_july) == 0)) {
  count_matrix_july <- count_matrix_july[-which(rowSums(count_matrix_july) == 0),]
}

# Run GLM-PCA
glmpca_res_july <- glmpca(
  Y = count_matrix_july, 
  L = 2,           # Extract PC1-2
  fam = "nb"       # "nb" = Negative Binomial, "poi" = Poisson
)

# Extract coordinates (`factors`)
pca_coords_july <- metadata_july |> 
  rownames_to_column(var = "SampleID") |> 
  left_join(rownames_to_column(glmpca_res_july$factors, var = "SampleID"), ., by = "SampleID")

# Build the plot
plot_glmpca_host <- ggplot(pca_coords_july, aes(x = dim1, y = dim2, color = host_plant)) +
  geom_point(size = 3, alpha = 0.85) +
  # Add stat_ellipse if you want 95% confidence intervals around your groups
  stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.075, show.legend = FALSE) +
  labs(
    # caption = "Dimensionality reduction on raw counts using a Negative Binomial model",
    x = "GLM-PC1", y = "GLM-PC2",
    color = "Host plant"
  ) +
  scale_color_manual(values = color.by.plant) +  
  scale_fill_manual(values = color.by.plant) +  
  theme(
    legend.position = "bottom"
  ) 

plot_glmpca_host


##########
# For month (among GRT samples)

# Extract raw counts
count_matrix_grt <- as(otu_table(ps_grt), "matrix")

# Ensure the matrix contains integers
count_matrix_grt <- round(count_matrix_grt)

# Filter out zero-count rows
if (any(rowSums(count_matrix_grt) == 0)) {
  count_matrix_grt <- count_matrix_grt[-which(rowSums(count_matrix_grt) == 0),]
}

# Run GLM-PCA
glmpca_res_grt <- glmpca(
  Y = count_matrix_grt, 
  L = 2,           
  fam = "nb"       
)

# Extract coordinates (`factors`)
pca_coords_grt <- metadata_grt |> 
  rownames_to_column(var = "SampleID") |> 
  left_join(rownames_to_column(glmpca_res_grt$factors, var = "SampleID"), ., by = "SampleID")

# Build the plot
plot_glmpca_month <- ggplot(pca_coords_grt, aes(x = dim1, y = dim2, color = month)) +
  geom_point(size = 3, alpha = 0.85) +
  # Add stat_ellipse if you want 95% confidence intervals around your groups
  stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  labs(
    # caption = "Dimensionality reduction on raw counts using a Negative Binomial model",
    x = "GLM-PC1", y = "GLM-PC2",
    color = "Month"
  ) +
  scale_color_manual(values = color.by.month) +  
  scale_fill_manual(values = color.by.month) +  
  theme(
    legend.position = "bottom"
  ) 

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
  count_matrix_july, glmpca_res_july, pca_coords_july, plot_glmpca_host,
  count_matrix_grt, glmpca_res_grt, pca_coords_grt, plot_glmpca_month,
  glmpca_figure,
  file = "Jhae.microbiome.glmpca.RData"
)
