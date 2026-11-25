################################################################################
# Beta diversity analysis

setwd("~/Documents/3.research/Jhae.microbiome")

# Clear memory
# {
#   rm(list=ls())
#   .rs.restartR()
#   gc()
# }

{
  library(tidyverse)
  library(patchwork)
  library(vegan)
}

load("Jhae.microbiome.RData", verbose = TRUE)
# load("beta.analysis.RData", verbose = TRUE)

################################################################################
# Distance matrices and ordinations

{
  # Save current warning level and turn warnings off
  old_warn <- getOption("warn")
  options(warn = -1)

  # Bray-Curtis distance
  dist_main_bray <- distance(ps_main, method = "bray")
  dist_july_bray <- distance(ps_july, method = "bray")
  dist_grt_bray <- distance(ps_grt, method = "bray")

  ord_main_bray <- ordinate(ps_main, method = "PCoA", distance = dist_main_bray)
  ord_july_bray <- ordinate(ps_july, method = "PCoA", distance = dist_july_bray)
  ord_grt_bray <- ordinate(ps_grt, method = "PCoA", distance = dist_grt_bray)

  # unweighted UniFrac distance
  asv.tree <- phytools::read.newick("qiime/16s.Jhae.rooted-tree.tre")
  ps_main_tree <- merge_phyloseq(ps_main, asv.tree)
  ps_july_tree <- merge_phyloseq(ps_july, asv.tree)
  ps_grt_tree <- merge_phyloseq(ps_grt, asv.tree)

  dist_main_unifrac <- distance(ps_main_tree, method = "unifrac")
  dist_july_unifrac <- distance(ps_july_tree, method = "unifrac")
  dist_grt_unifrac <- distance(ps_grt_tree, method = "unifrac")

  ord_main_unifrac <- ordinate(ps_main_tree, method = "PCoA", distance = dist_main_unifrac)
  ord_july_unifrac <- ordinate(ps_july_tree, method = "PCoA", distance = dist_july_unifrac)
  ord_grt_unifrac <- ordinate(ps_grt_tree, method = "PCoA", distance = dist_grt_unifrac)

  # Rao's distance (via DPCoA)
  ord_main_DPCoA <- DPCoA(ps_main_tree, correction = "cailliez")
  ord_july_DPCoA <- DPCoA(ps_july_tree, correction = "cailliez")
  ord_grt_DPCoA <- DPCoA(ps_grt_tree, correction = "cailliez")

  dist_main_rao <- ord_main_DPCoA$RaoDis
  dist_july_rao <- ord_july_DPCoA$RaoDis
  dist_grt_rao <- ord_grt_DPCoA$RaoDis

  # Restore original warning setting (default is 0)
  options(warn = old_warn)
}

################################################################################
# Preliminary model exploration

test.betadisp <- function (x, f, exclude.level = NULL, iter = 1e4-1, verbose = FALSE) {
  if (!is.null(exclude.level)) {
    i <- which(f == exclude.level)
    f <- f[-i]
    x <- as.dist(as.matrix(x)[-i,-i])
  }
  test <- betadisper(x, f)
  if (verbose) {
    print(permutest(test, permutations = iter))
    boxplot(test$distances ~ test$group)
  } else {
    x <- permutest(test, permutations = iter)
    return(x$tab[1,6])
  }
}


{
  categorical.factors <- c("sex", "morph","host_plant", "month")
  individual.factors <- c("sex", "morph")
  environmental.factors <- c(
    "host_plant", "month",
    "developed_open_space_m2",
    "developed_low_intensity_m2", "developed_medium_intensity_m2", "developed_high_intensity_m2",
    "barren_land_m2", "evergreen_forest_m2", "scrub_m2", "grassland_m2",
    "pasture_m2", "cropland_m2", "woody_wetlands_m2", "herbaceous_wetlands_m2"
  )
  aggregate.factors <- c("intensive_development", "green_space")
  all.factors <- c(individual.factors, environmental.factors, aggregate.factors)

  all.columns <- which(colnames(metadata_main) %in% all.factors)
  x <- rep(NA,length(all.columns))
  single.factor.adonis.results <- data.frame(
    factor = x,
    bray.Bdsp = x, bray.R2 = x, bray.F = x, bray.p = x, bray.q = x,
    unif.Bdsp = x, unif.R2 = x, unif.F = x, unif.p = x, unif.q = x,
    rao.Bdsp = x,  rao.R2 = x,  rao.F = x,  rao.p = x,  rao.q = x
  )

  iterations <- 1e4-1

  # Save current warning level and turn warnings off
  old_warn <- getOption("warn")
  options(warn = -1)

  for (i in 1:length(all.columns)) {
    col.i <- colnames(metadata_main)[all.columns[i]]
    data.i <- metadata_grt
    dist.bray.i <- dist_grt_bray
    dist.unif.i <- dist_grt_unifrac
    dist.rao.i <- dist_grt_rao
    if (col.i == "host_plant") {
      data.i <- metadata_july
      dist.bray.i <- dist_july_bray
      dist.unif.i <- dist_july_unifrac
      dist.rao.i <- dist_july_rao
    }
    factor.i <- data.i[,all.columns[i]]
    data.i$factor.i <- factor.i
    single.factor.adonis.results$factor[i] <- col.i

    if (col.i %in% categorical.factors) {
      single.factor.adonis.results$bray.Bdsp[i] <- test.betadisp(dist.bray.i, factor.i)
      single.factor.adonis.results$unif.Bdsp[i] <- test.betadisp(dist.unif.i, factor.i)
      single.factor.adonis.results$rao.Bdsp[i] <- test.betadisp(dist.rao.i, factor.i)
    }

    model.bray.i <- adonis2(
      dist.bray.i ~ factor.i, data = data.i,
      by = "margin", permutations = iterations
    )
    single.factor.adonis.results$bray.R2[i] <- model.bray.i$R2[1]
    single.factor.adonis.results$bray.F[i] <- model.bray.i$`F`[1]
    single.factor.adonis.results$bray.p[i] <- model.bray.i$`Pr(>F)`[1]

    model.unif.i <- adonis2(
      dist.unif.i ~ factor.i, data = data.i,
      by = "margin", permutations = iterations
    )
    single.factor.adonis.results$unif.R2[i] <- model.unif.i$R2[1]
    single.factor.adonis.results$unif.F[i] <- model.unif.i$`F`[1]
    single.factor.adonis.results$unif.p[i] <- model.unif.i$`Pr(>F)`[1]

    model.rao.i <- adonis2(
      dist.rao.i ~ factor.i, data = data.i,
      by = "margin", permutations = iterations
    )
    single.factor.adonis.results$rao.R2[i] <- model.rao.i$R2[1]
    single.factor.adonis.results$rao.F[i] <- model.rao.i$`F`[1]
    single.factor.adonis.results$rao.p[i] <- model.rao.i$`Pr(>F)`[1]
  }

  # Restore original warning setting (default is 0)
  options(warn = old_warn)

  # Clean up
  rm(col.i, factor.i, data.i, dist.bray.i, dist.unif.i, dist.rao.i, model.bray.i, model.unif.i, model.rao.i)

  # Apply FDR correction to individual factors
  i <- 1:length(individual.factors)
  single.factor.adonis.results$bray.q[i] <- p.adjust(single.factor.adonis.results$bray.p[i], method = "BH")
  single.factor.adonis.results$unif.q[i] <- p.adjust(single.factor.adonis.results$unif.p[i], method = "BH")
  single.factor.adonis.results$rao.q[i] <- p.adjust(single.factor.adonis.results$rao.p[i], method = "BH")

  # Apply FDR correction to environmental factors
  i <- (length(individual.factors)+1):(length(individual.factors)+length(environmental.factors))
  single.factor.adonis.results$bray.q[i] <- p.adjust(single.factor.adonis.results$bray.p[i], method = "BH")
  single.factor.adonis.results$unif.q[i] <- p.adjust(single.factor.adonis.results$unif.p[i], method = "BH")
  single.factor.adonis.results$rao.q[i] <- p.adjust(single.factor.adonis.results$rao.p[i], method = "BH")

  # Apply FDR correction to aggregate factors
  i <- (length(all.factors)-length(aggregate.factors)+1):length(all.factors)
  single.factor.adonis.results$bray.q[i] <- p.adjust(single.factor.adonis.results$bray.p[i], method = "BH")
  single.factor.adonis.results$unif.q[i] <- p.adjust(single.factor.adonis.results$unif.p[i], method = "BH")
  single.factor.adonis.results$rao.q[i] <- p.adjust(single.factor.adonis.results$rao.p[i], method = "BH")

  # Significant results?
  single.factor.adonis.results |> dplyr::filter(bray.q < 0.05) |> dplyr::select(factor, starts_with("bray")) |> arrange(desc(bray.R2))
  #       factor bray.Bdsp    bray.R2   bray.F bray.p bray.q
  # 1 host_plant    0.6669 0.07421676 7.054648  1e-04 0.0014

  single.factor.adonis.results |> dplyr::filter(unif.q < 0.05) |> dplyr::select(factor, starts_with("unif")) |> arrange(desc(unif.R2))
  #                           factor unif.Bdsp    unif.R2   unif.F unif.p     unif.q
  # 1                     host_plant    0.0001 0.05321624 4.946250 0.0001 0.00070000
  # 2                       scrub_m2    0.0077 0.01709503 2.365360 0.0002 0.00070000
  # 3                          month    0.8564 0.01538674 2.125298 0.0002 0.00070000
  # 4                   grassland_m2    0.3780 0.01501169 2.072705 0.0001 0.00070000
  # 5                    green_space    0.0001 0.01436005 1.981421 0.0006 0.00120000
  # 6                    cropland_m2    0.0950 0.01432843 1.976993 0.0003 0.00084000
  # 7         herbaceous_wetlands_m2    0.0001 0.01302953 1.795409 0.0009 0.00210000
  # 8  developed_medium_intensity_m2    0.0001 0.01266014 1.743857 0.0017 0.00340000
  # 9          intensive_development    0.0001 0.01142961 1.572399 0.0063 0.00630000
  # 10                    pasture_m2    0.0001 0.01064009 1.462614 0.0124 0.02068889
  # 11       developed_open_space_m2    0.0001 0.01059362 1.456158 0.0133 0.02068889
  # 12    developed_low_intensity_m2    0.0001 0.01037390 1.425640 0.0153 0.02142000

  single.factor.adonis.results |> dplyr::filter(rao.q < 0.05) |> dplyr::select(factor, starts_with("rao")) |> arrange(desc(rao.R2))
  #                          factor rao.Bdsp     rao.R2     rao.F  rao.p      rao.q
  # 1                    host_plant   0.3133 0.12067482 12.076744 0.0001 0.00140000
  # 2        herbaceous_wetlands_m2   0.0228 0.05208945  7.473454 0.0019 0.01330000
  # 3         intensive_development   0.0266 0.04704449  6.713903 0.0043 0.00860000
  # 4 developed_medium_intensity_m2   0.0301 0.04702744  6.711351 0.0034 0.01586667
  # 5   developed_high_intensity_m2   0.0296 0.04179378  5.931870 0.0060 0.02100000
  # 6       developed_open_space_m2   0.0257 0.03426685  4.825652 0.0141 0.03948000
  # 7                   green_space   0.0285 0.02360682  3.288151 0.0433 0.04330000

  write.csv(single.factor.adonis.results, "tables/TableS6.single.factor.permanova.results.csv", row.names = FALSE)
}


################################################################################
# PCoA / DPCoA Plots
{
  {
    theme_set(theme_bw())
    update_theme(legend.position = "none")
    update_theme(axis.ticks = element_blank())
  }

  {
   # Host plant / July
    pcoa.july.bray.host_plant <- plot_ordination(ps_july, ord_july_bray, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed() +
      labs(title="Bray-Curtis dissimilarity")
    pcoa.july.bray.host_plant

    pcoa.july.unifrac.host_plant <- plot_ordination(ps_july_tree, ord_july_unifrac, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.plant) +
      coord_fixed() +
      labs(title="unweighted UniFrac")
    pcoa.july.unifrac.host_plant

    dpcoa.july.rao.host_plant <- plot_ordination(ps_july_tree, ord_july_DPCoA, color = "host_plant") +
      stat_ellipse(aes(fill = host_plant), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.plant) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(
        name = "Host plant", values = color.by.plant,
        labels = c("BV" = "C. corindum", "GRT" = "K. paniculata")
      ) +
      coord_fixed() +
      labs(title="Rao's dissimilarity") +
      theme(
        legend.position = "right",
        legend.text = element_text(face = "italic")
      )
    dpcoa.july.rao.host_plant


    # Month / GRT
    pcoa.grt.bray.month <- plot_ordination(ps_grt, ord_grt_bray, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.grt.bray.month

    pcoa.grt.unifrac.month <- plot_ordination(ps_grt_tree, ord_grt_unifrac, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = NULL, values = color.by.month) +
      coord_fixed()
    pcoa.grt.unifrac.month

    dpcoa.grt.rao.month <- plot_ordination(ps_grt_tree, ord_grt_DPCoA, color = "month") +
      stat_ellipse(aes(fill = month), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
      scale_fill_manual(name = NULL, values = color.by.month) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(name = "Month of Ccllection", values = color.by.month) +
      coord_fixed() +
      theme(
        legend.position = "right"
      )
    
    dpcoa.grt.rao.month

    # green_space / GRT

    # tertiles
    metadata_grt$green_space_tertile <- cut(
      metadata_grt$green_space,
      breaks = quantile(metadata_grt$green_space,
                        probs = seq(0, 1, 1/3), na.rm = TRUE),
      include.lowest = TRUE, labels = c("Low", "Medium", "High")
    )
    sample_data(ps_grt) <- sample_data(metadata_grt)
    sample_data(ps_grt_tree) <- sample_data(metadata_grt)

    pcoa.grt.bray.green_space <- plot_ordination(ps_grt, ord_grt_bray, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[2], high = color.by.urban[1]) +
      stat_ellipse(aes(group = green_space_tertile, fill = green_space_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(values = c("Low" = color.by.urban[3], "Medium" = color.by.urban[2], "High" = color.by.urban[1])) +
      coord_fixed()
    pcoa.grt.bray.green_space

    pcoa.grt.unifrac.green_space <- plot_ordination(ps_grt_tree, ord_grt_unifrac, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[2], high = color.by.urban[1]) +
      stat_ellipse(aes(group = green_space_tertile, fill = green_space_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(values = c("Low" = color.by.urban[3], "Medium" = color.by.urban[2], "High" = color.by.urban[1])) +
      coord_fixed()
    pcoa.grt.unifrac.green_space

    dpcoa.grt.rao.green_space <- plot_ordination(ps_grt_tree, ord_grt_DPCoA, color = "green_space") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(
        name = "Green space (ha)", 
        breaks = seq(1e5,6e5,1e5), labels = 1:6,
        low = color.by.urban[2], high = color.by.urban[1]
        ) +
      stat_ellipse(aes(group = green_space_tertile, fill = green_space_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(
        name = "Green space tertile", 
        values = c("Low" = color.by.urban[3], "Medium" = color.by.urban[2], "High" = color.by.urban[1])
      ) +
      coord_fixed() +
      theme(
        legend.position = "right"
      ) +
      guides(
        color = guide_legend(order = 1),
        fill = guide_legend(order = 2)
      )
    
    dpcoa.grt.rao.green_space

    # intensive_development / GRT

    # tertiles
    metadata_grt$intensive_development_tertile <- cut(
      metadata_grt$intensive_development,
      breaks = quantile(metadata_grt$intensive_development,
                        probs = seq(0, 1, 1/3), na.rm = TRUE),
      include.lowest = TRUE, labels = c("Low", "Medium", "High")
    )
    sample_data(ps_grt) <- sample_data(metadata_grt)
    sample_data(ps_grt_tree) <- sample_data(metadata_grt)

    pcoa.grt.bray.intensive_development <- plot_ordination(ps_grt, ord_grt_bray, color = "intensive_development") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[1], high = color.by.urban[2]) +
      stat_ellipse(aes(group = intensive_development_tertile, fill = intensive_development_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(values = c("Low" = color.by.urban[1], "Medium" = color.by.urban[2], "High" = color.by.urban[3])) +
      coord_fixed()
    pcoa.grt.bray.intensive_development

    pcoa.grt.unifrac.intensive_development <- plot_ordination(ps_grt_tree, ord_grt_unifrac, color = "intensive_development") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = color.by.urban[1], high = color.by.urban[2]) +
      stat_ellipse(aes(group = intensive_development_tertile, fill = intensive_development_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(values = c("Low" = color.by.urban[1], "Medium" = color.by.urban[2], "High" = color.by.urban[3])) +
      coord_fixed()
    pcoa.grt.unifrac.intensive_development

    dpcoa.grt.rao.intensive_development <- plot_ordination(ps_grt_tree, ord_grt_DPCoA, color = "intensive_development") +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(
        name = "Intensive development (ha)", 
        breaks = seq(1e5,6e5,1e5), labels = 1:6,
        low = color.by.urban[1], high = color.by.urban[2]
      ) +
      
      stat_ellipse(aes(group = intensive_development_tertile, fill = intensive_development_tertile),
                   geom = "polygon", alpha = 0.2, color = "gray95") +
      scale_fill_manual(
        name = "Intensive development terile", 
        values = c("Low" = color.by.urban[1], "Medium" = color.by.urban[2], "High" = color.by.urban[3])
      ) +
      coord_fixed() +
      theme(
        legend.position = "right"
      ) +
      guides(
        color = guide_legend(order = 1),
        fill = guide_legend(order = 2)
      )
    dpcoa.grt.rao.intensive_development

  }

  figure_PCoA <- (
    (pcoa.july.bray.host_plant | pcoa.july.unifrac.host_plant | dpcoa.july.rao.host_plant) /
    (pcoa.grt.bray.month | pcoa.grt.unifrac.month | dpcoa.grt.rao.month) /
    (pcoa.grt.bray.intensive_development | pcoa.grt.unifrac.intensive_development | dpcoa.grt.rao.intensive_development) /
    (pcoa.grt.bray.green_space | pcoa.grt.unifrac.green_space | dpcoa.grt.rao.green_space) 
  ) +
    plot_annotation(tag_levels = 'A')

  ggsave(
    filename = "plots/FigS3.PCoA.png",
    plot = figure_PCoA,
    width = 6.5, height = (6.5)*(4/3), scale = 1.8,
    dpi = 600, bg = "white"
  )

}


##########
# Full PERMANOVA Models
# Partial sums of squares (the equivalent of Type III ANOVA)

{
  iterations <- 1e4-1
  perm_bray_full <- adonis2(
    dist_main_bray ~ host_plant + month + intensive_development,
    data = metadata_main,
    by = "margin",
    permutations = iterations
  )
  perm_bray_full
  #                        Df SumOfSqs      R2      F Pr(>F)
  # host_plant              1    2.120 0.04834 7.2927 0.0001 ***
  # month                   1    0.859 0.01958 2.9539 0.0072 **
  # intensive_development   1    0.691 0.01574 2.3753 0.0251 *
  # Residual              139   40.415 0.92130
  # Total                 142   43.867 1.00000

  perm_unifrac_full <- adonis2(
    dist_main_unifrac ~ host_plant + month + intensive_development,
    data = metadata_main,
    by = "margin",
    permutations = iterations
  )
  perm_unifrac_full
  #                        Df SumOfSqs      R2      F Pr(>F)
  # host_plant              1    1.901 0.03005 4.4352 0.0001 ***
  # month                   1    1.002 0.01583 2.3367 0.0001 ***
  # intensive_development   1    0.689 0.01088 1.6064 0.0026 **
  # Residual              139   59.577 0.94175
  # Total                 142   63.262 1.00000

  perm_rao_full <- adonis2(
    dist_main_rao ~ host_plant + month + intensive_development,
    data = metadata_main,
    by = "margin",
    permutations = iterations
  )
  perm_rao_full
  #                        Df SumOfSqs      R2       F Pr(>F)
  # host_plant              1   1.0758 0.08613 13.8877 0.0001 ***
  # month                   1   0.1811 0.01450  2.3377 0.0980 .
  # intensive_development   1   0.5327 0.04265  6.8773 0.0022 **
  # Residual              139  10.7676 0.86202
  # Total                 142  12.4912 1.00000
}

{
  adonis.list <- list(
    Bray_Curtis__full_model = perm_bray_full,
    UniFrac__full_model = perm_unifrac_full,
    Rao__full_model = perm_rao_full
  )
  
  # Combine into one table, adjust p-values, and add stars
  full.model.adonis.results <- bind_rows(
    lapply(adonis.list, \(x) rownames_to_column(as.data.frame(x), var = "Term")),
    .id = "listname"
  ) |>
    mutate(
      data_set = str_remove(listname, "__.*"),
      # model = str_remove(listname, "^.*__"),
      p_adj = p.adjust(`Pr(>F)`, method = "BH"),
      signif = case_when(
        p_adj <= 0.001 ~ "***",
        p_adj <= 0.01  ~ "**",
        p_adj <= 0.05  ~ "*",
        p_adj <= 0.1   ~ ".",
        TRUE           ~ ""
      )
    ) |> 
    select(-listname) |> 
    relocate(data_set)
  
  full.model.adonis.results
  
  write.csv(full.model.adonis.results, "tables/TableS7.full.model.permanova.results.csv", row.names = FALSE)
}


################################################################################
# Save the key objects
save(
  dist_main_bray, dist_july_bray, dist_grt_bray,
  ord_main_bray, ord_july_bray, ord_grt_bray,
  asv.tree,
  ps_main_tree, ps_july_tree, ps_grt_tree,
  dist_main_unifrac, dist_july_unifrac, dist_grt_unifrac,
  ord_main_unifrac, ord_july_unifrac, ord_grt_unifrac,
  ord_main_DPCoA, ord_july_DPCoA, ord_grt_DPCoA,
  dist_main_rao, dist_july_rao, dist_grt_rao,
  test.betadisp,
  single.factor.adonis.results,
  pcoa.july.bray.host_plant, pcoa.july.unifrac.host_plant, dpcoa.july.rao.host_plant,
  pcoa.grt.bray.month, pcoa.grt.unifrac.month, dpcoa.grt.rao.month,
  metadata_grt,
  pcoa.grt.bray.green_space, pcoa.grt.unifrac.green_space, dpcoa.grt.rao.green_space,
  pcoa.grt.bray.intensive_development, pcoa.grt.unifrac.intensive_development, dpcoa.grt.rao.intensive_development,
  figure_PCoA,
  perm_bray_full, perm_unifrac_full, perm_rao_full,
  full.model.adonis.results,
  file = "beta.analysis.RData"
)
