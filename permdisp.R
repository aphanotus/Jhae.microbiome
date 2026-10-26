################################################################################
# Homogeneity of multivariate dispersions (PERMDISP)

# Since green space is a significant factor in PERMANOVAs but
# no individual ASV (or Family-level taxa) emerge as differentially abundant
# in ALDEx models, perhaps it's differences in variance, not means that
# drive the global community contrasts in PERMANOVA?

{
  library(tidyverse)
  library(phyloseq)
  library(vegan)
}

load("Jhae.microbiome.RData", verbose = TRUE)

# Categorize the continuous green_space variable into tertiles (Low, Med, High)
metadata_rare_grt <- metadata_rare_grt %>%
  mutate(
    green_space_cat = ntile(green_space, 3),
    green_space_cat = factor(
      green_space_cat,
      levels = c(1, 2, 3),
      labels = c("Low (Urban)", "Medium", "High (Rural)"))
  )


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
    title = "Bray-Curtis Distance",
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



dispersion_bray <- betadisper(dist_grt_bray, metadata_rare_grt$green_space_cat)
permutest_bray <- permutest(dispersion_bray, permutations = 1e5-1)
print(permutest_bray)
#            Df  Sum Sq Mean Sq      F N.Perm  Pr(>F)   
# Groups      2  1.4526 0.72629 7.3334  99999 0.00104 **
# Residuals 136 13.4691 0.09904   

dispersion_uwUF <- betadisper(dist_grt_uwUF, metadata_rare_grt$green_space_cat)
permutest_uwUF <- permutest(dispersion_uwUF, permutations = 1e5-1)
print(permutest_uwUF)
#            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.4246 0.212289 7.086  99999 0.00155 **
# Residuals 136 4.0744 0.029959   

dispersion_wUF <- betadisper(dist_grt_wUF, metadata_rare_grt$green_space_cat)
permutest_wUF <- permutest(dispersion_wUF, permutations = 1e5-1)
print(permutest_wUF)
#            Df  Sum Sq  Mean Sq      F N.Perm  Pr(>F)
# Groups      2 0.10743 0.053715 6.3726  99999 0.00223 **
# Residuals 136 1.14635 0.008429   

  

beta_dispersion_figure <- (
  plot_dispersion_bray_by_green_space | plot_dispersion_uwUF_by_green_space | plot_dispersion_wUF_by_green_space
) +
  plot_annotation(tag_levels = 'A')
beta_dispersion_figure

ggsave(
  filename = "plots/FigX.beta.dispersion.by.green.space.png",
  plot = beta_dispersion_figure,
  width = 6.5, height = 2.5, scale = 1.75,
  dpi = 600, bg = "white"
)


################################################################################
# Save the key objects
save(
  beta_dispersion_figure,
  file = "Jhae.microbiome.permdisp.RData"
)
