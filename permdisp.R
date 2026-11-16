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
load("beta.analysis.RData", verbose = TRUE)

# Categorize the continuous intensive_development variable into tertiles (Low, Medium, High)

metadata_grt$intensive_development_tertile <- cut(
  metadata_grt$intensive_development, 
  breaks = quantile(metadata_grt$intensive_development, 
                    probs = seq(0, 1, 1/3), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Low", "Medium", "High")
) 
sample_data(ps_grt) <- sample_data(metadata_grt)
sample_data(ps_grt_tree) <- sample_data(metadata_grt)

iterations <- 1e5-1

intensive_development_disp_bray <- betadisper(dist_grt_bray, metadata_grt$intensive_development_tertile)
intensive_development_permtest_bray <- permutest(intensive_development_disp_bray, permutations = iterations)
print(intensive_development_permtest_bray)
#            Df Sum Sq Mean Sq      F N.Perm  Pr(>F)  
# Groups      2 0.2361 0.11803 2.8379  99999 0.06286 .
# Residuals 135 5.6146 0.04159   

intensive_development_disp_unif <- betadisper(dist_grt_unifrac, metadata_grt$intensive_development_tertile)
intensive_development_permtest_unif <- permutest(intensive_development_disp_unif, permutations = iterations)
print(intensive_development_permtest_unif)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.004268 0.0021342 1.0028  99999 0.3701
# Residuals 135 0.287323 0.0021283  

intensive_development_disp_rao <- betadisper(dist_grt_rao, metadata_grt$intensive_development_tertile)
intensive_development_permtest_rao <- permutest(intensive_development_disp_rao, permutations = iterations)
print(intensive_development_permtest_rao)
#            Df Sum Sq Mean Sq      F N.Perm  Pr(>F)   
# Groups      2 0.4665 0.23325 5.5548  99999 0.00451 **
# Residuals 135 5.6687 0.04199     

{
  theme_set(theme_bw())
  update_theme(legend.position = "none")
  update_theme(axis.ticks = element_blank())
  update_theme(panel.grid.major.x = element_blank())
}

# Extract the distance-to-centroid values to see the direction of the variance
dispersion_bray_for_plot <- data.frame(
  sample = names(intensive_development_disp_bray$distances),
  distance_to_centroid = intensive_development_disp_bray$distances,
  intensive_development = metadata_grt$intensive_development_tertile
)

# Plot the differences in variance of Beta Dispersion by intensive_development_tertile
plot_intensive_development_disp_bray <- dispersion_bray_for_plot |>
  ggplot(aes(x = intensive_development, y = distance_to_centroid, fill = intensive_development)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
  ylim(c(0,1)) +
  labs(
    title = "Bray-Curtis dissimilarity",
    x = "Urban development ",
    y = "Distance to Group Centroid (Community Variance)"
  )

plot_intensive_development_disp_bray


dispersion_unifrac_for_plot <- data.frame(
  sample = names(intensive_development_disp_unif$distances),
  distance_to_centroid = intensive_development_disp_unif$distances,
  intensive_development = metadata_grt$intensive_development_tertile
)

plot_intensive_development_disp_unifrac <- dispersion_unifrac_for_plot |>
  ggplot(aes(x = intensive_development, y = distance_to_centroid, fill = intensive_development)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
  ylim(c(0,1)) +
  labs(
    title = "unweighted UniFrac",
    x = "Urban development ",
    y = "Distance to Group Centroid (Community Variance)"
  )

plot_intensive_development_disp_unifrac


dispersion_rao_for_plot <- data.frame(
  sample = names(intensive_development_disp_rao$distances),
  distance_to_centroid = intensive_development_disp_rao$distances,
  intensive_development = metadata_grt$intensive_development_tertile
)

plot_intensive_development_disp_rao <- dispersion_rao_for_plot |>
  ggplot(aes(x = intensive_development, y = distance_to_centroid, fill = intensive_development)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c(color.by.urban[1],"palegreen3",color.by.urban[2])) +
  ylim(c(0,1)) +
  labs(
    title = "Rao's dissimilarity",
    x = "Urban development ",
    y = "Distance to Group Centroid (Community Variance)"
  )

plot_intensive_development_disp_rao



beta_dispersion_figure <- (
  plot_intensive_development_disp_bray | plot_intensive_development_disp_unifrac | plot_intensive_development_disp_rao
) +
  plot_annotation(tag_levels = 'A')
beta_dispersion_figure

ggsave(
  filename = "plots/Fig6.beta.dispersion.png",
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
