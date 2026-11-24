################################################################################
# Alpha diversity analysis

# Clear memory
# {
#   rm(list=ls())
#   .rs.restartR()
#   gc()
# }

setwd("~/Documents/3.research/Jhae.microbiome")

{
  library(tidyverse)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(rstatix)
  library(effectsize)
}

load("Jhae.microbiome.RData", verbose = TRUE)
# load("alpha.RData", verbose = TRUE)

# Save the raw alpha diversity
{
  alpha.table <- metadata_main |>
    rownames_to_column(var = "Sample_ID") |>
    mutate(
      intensive_development_ha = intensive_development / 1e5,
      green_space_ha = green_space / 1e5,
    ) |> 
    mutate(host_plant = str_replace(host_plant, "BV", "C. corindum")) |> 
    mutate(host_plant = str_replace(host_plant, "GRT", "K. paniculata")) |> 
    dplyr::select(
      Sample_ID, sex, morph, month, host_plant, city, 
      intensive_development_ha, green_space_ha,
      observed_ASVs, shannon, simpson
    )

  write.csv(alpha.table, "tables/TableS4.alpha.diversity.csv", row.names = FALSE)
}

################################################################################
# Alpha stats
{
  alpha.columns <- c("observed_ASVs", "shannon", "simpson")
  alpha.labels <- c("Observed ASVs", "Shannon Index", "Simpson Index")
  names(alpha.labels) <- alpha.columns

  categorical.factors <- c("sex", "morph", "host_plant", "month")
  continuous.factors <- c(
    "intensive_development", "green_space",
    "developed_open_space_m2",
    "developed_low_intensity_m2", "developed_medium_intensity_m2", "developed_high_intensity_m2",
    "barren_land_m2", "evergreen_forest_m2", "scrub_m2", "grassland_m2",
    "pasture_m2", "cropland_m2", "woody_wetlands_m2", "herbaceous_wetlands_m2"
  )

  # Categorical factors
  {
    categorical.columns <- which(colnames(metadata_main) %in% categorical.factors)
    x <- rep(NA,length(categorical.columns))
    alpha.wrst.results <- data.frame(
      factor = x,
      obsd.W = x, obsd.r = x, obsd.p = x, obsd.q = x,
      shan.W = x, shan.r = x, shan.p = x, shan.q = x,
      simp.W = x, simp.r = x, simp.p = x, simp.q = x
    )

    for (i in 1:length(categorical.columns)) {
      col.i <- colnames(metadata_main)[categorical.columns[i]]
      data.i <- metadata_grt
      if (col.i == "host_plant") { data.i <- metadata_july }
      factor.i <- data.i[,categorical.columns[i]]
      data.i$factor.i <- factor.i
      alpha.wrst.results$factor[i] <- col.i

      wrst.obsd.i <- wilcox.test(observed_ASVs ~ get(col.i), data = data.i)
      alpha.wrst.results$obsd.W[i] <- wrst.obsd.i$statistic
      alpha.wrst.results$obsd.p[i] <- wrst.obsd.i$p.value
      wrst.obsd.i.eff <- wilcox_effsize(formula = observed_ASVs ~ factor.i, data = data.i)
      alpha.wrst.results$obsd.r[i] <- wrst.obsd.i.eff$effsize

      wrst.shan.i <- wilcox.test(shannon ~ get(col.i), data = data.i)
      alpha.wrst.results$shan.W[i] <- wrst.shan.i$statistic
      alpha.wrst.results$shan.p[i] <- wrst.shan.i$p.value
      wrst.shan.i.eff <- wilcox_effsize(formula = shannon ~ factor.i, data = data.i)
      alpha.wrst.results$shan.r[i] <- wrst.shan.i.eff$effsize

      wrst.simp.i <- wilcox.test(simpson ~ get(col.i), data = data.i)
      alpha.wrst.results$simp.W[i] <- wrst.simp.i$statistic
      alpha.wrst.results$simp.p[i] <- wrst.simp.i$p.value
      wrst.simp.i.eff <- wilcox_effsize(formula = simpson ~ factor.i, data = data.i)
      alpha.wrst.results$simp.r[i] <- wrst.simp.i.eff$effsize
    }

    # Clean up
    rm(col.i, factor.i, data.i, wrst.obsd.i, wrst.shan.i, wrst.simp.i, wrst.obsd.i.eff, wrst.shan.i.eff, wrst.simp.i.eff)

    # Apply FDR correction
    alpha.wrst.results$obsd.q <- p.adjust(alpha.wrst.results$obsd.p, method = "BH")
    alpha.wrst.results$shan.q <- p.adjust(alpha.wrst.results$shan.p, method = "BH")
    alpha.wrst.results$simp.q <- p.adjust(alpha.wrst.results$simp.p, method = "BH")

    # Significant results?
    alpha.wrst.results |> dplyr::filter(obsd.q < 0.05) # month
    alpha.wrst.results |> dplyr::filter(shan.q < 0.05) # month, host_plant
    alpha.wrst.results |> dplyr::filter(simp.q < 0.05) # month, host_plant

  }

  # Continuous factors
  {
    continuous.columns <- which(colnames(metadata_main) %in% continuous.factors)
    x <- rep(NA,length(continuous.columns))
    alpha.lmer.results <- data.frame(
      factor = x,
      obsd.t = x, obsd.stdB = x, obsd.p = x, obsd.q = x,
      shan.t = x, shan.stdB = x, shan.p = x, shan.q = x,
      simp.t = x, simp.stdB = x, simp.p = x, simp.q = x
    )

    # Save current warning level and turn warnings off
    old_warn <- getOption("warn")
    options(warn = -1)

    for (i in 1:length(continuous.columns)) {
      col.i <- colnames(metadata_main)[continuous.columns[i]]
      data.i <- metadata_grt
      factor.i <- data.i[,continuous.columns[i]]
      data.i$factor.i <- factor.i
      alpha.lmer.results$factor[i] <- col.i

      lmer.obsd.i <- lmer(observed_ASVs ~ factor.i + (1 | site_name), data = data.i)
      summ.obsd.i <- summary(lmer.obsd.i)
      eff.obsd.i <- standardize_parameters(lmer.obsd.i)
      alpha.lmer.results$obsd.t[i] <- summ.obsd.i$coefficients[2,4]
      alpha.lmer.results$obsd.p[i] <- summ.obsd.i$coefficients[2,5]
      alpha.lmer.results$obsd.stdB[i] <- eff.obsd.i$Std_Coefficient[2]

      lmer.shan.i <- lmer(shannon ~ factor.i + (1 | site_name), data = data.i)
      summ.shan.i <- summary(lmer.shan.i)
      eff.shan.i <- standardize_parameters(lmer.shan.i)
      alpha.lmer.results$shan.t[i] <- summ.shan.i$coefficients[2,4]
      alpha.lmer.results$shan.p[i] <- summ.shan.i$coefficients[2,5]
      alpha.lmer.results$shan.stdB[i] <- eff.shan.i$Std_Coefficient[2]

      lmer.simp.i <- lmer(simpson ~ factor.i + (1 | site_name), data = data.i)
      summ.simp.i <- summary(lmer.simp.i)
      eff.simp.i <- standardize_parameters(lmer.simp.i)
      alpha.lmer.results$simp.t[i] <- summ.simp.i$coefficients[2,4]
      alpha.lmer.results$simp.p[i] <- summ.simp.i$coefficients[2,5]
      alpha.lmer.results$simp.stdB[i] <- eff.simp.i$Std_Coefficient[2]
    }

    # Restore original warning setting (default is 0)
    options(warn = old_warn)

    # Clean up
    rm(old_warn, col.i, factor.i, data.i, lmer.obsd.i, summ.obsd.i, eff.obsd.i, lmer.shan.i, summ.shan.i, eff.shan.i, lmer.simp.i, summ.simp.i, eff.simp.i)

    # Apply FDR correction to individual factors
    i <- 1:(length(alpha.lmer.results$factor)-2)
    alpha.lmer.results$obsd.q[i] <- p.adjust(alpha.lmer.results$obsd.p[i], method = "BH")
    alpha.lmer.results$shan.q[i] <- p.adjust(alpha.lmer.results$shan.p[i], method = "BH")
    alpha.lmer.results$simp.q[i] <- p.adjust(alpha.lmer.results$simp.p[i], method = "BH")

    # Apply FDR correction to the two aggregate factors
    i <- (length(alpha.lmer.results$factor)-1):length(alpha.lmer.results$factor)
    alpha.lmer.results$obsd.q[i] <- p.adjust(alpha.lmer.results$obsd.p[i], method = "BH")
    alpha.lmer.results$shan.q[i] <- p.adjust(alpha.lmer.results$shan.p[i], method = "BH")
    alpha.lmer.results$simp.q[i] <- p.adjust(alpha.lmer.results$simp.p[i], method = "BH")

    # Significant results?
    alpha.lmer.results |> dplyr::filter(obsd.q < 0.05) # scrub_m2
    alpha.lmer.results |> dplyr::filter(shan.q < 0.05) # none
    alpha.lmer.results |> dplyr::filter(simp.q < 0.05) # none

    write.csv(alpha.lmer.results, "tables/TableSX.alpha.LMER.results.csv", row.names = FALSE)
  }
}

################################################################################
# Combine the stats 

{
  alpha.wrst.results.long <- alpha.wrst.results |>
    pivot_longer(
      cols = -factor,
      names_to = c("dataset", ".value"),
      names_sep = "\\."
    ) |>
    rename(
      test_statistic = W,
      effect_size    = r,
      p_value        = p,
      q_value        = q
    ) |>
    mutate(test = "Wilcoxon (W, r)")
  
  alpha.lmer.results.long <- alpha.lmer.results |>
    pivot_longer(
      cols = -factor,
      names_to = c("dataset", ".value"),
      names_sep = "\\."
    ) |>
    rename(
      test_statistic = t,
      effect_size    = stdB,
      p_value        = p,
      q_value        = q
    ) |>
    mutate(test = "LMM (t, Std. Beta)")
  
  # Combine, clean up dataset labels, and order columns
  alpha.combined <- bind_rows(alpha.wrst.results.long, alpha.lmer.results.long) |>
    mutate(
      dataset = recode(
        dataset,
        "obsd" = "Observed ASVs",
        "shan" = "Shannon's index",
        "simp" = "Simpson's index"
      ),
      dataset = factor(dataset, levels = c("Observed ASVs", "Shannon's index", "Simpson's index"))
    ) |> 
    mutate(
      signif = case_when(
        q_value <= 0.001 ~ "***",
        q_value <= 0.01  ~ "**",
        q_value <= 0.05  ~ "*",
        q_value <= 0.1   ~ ".",
        TRUE           ~ ""
      )
    ) |> 
    dplyr::select(factor, dataset, test, test_statistic, effect_size, p_value, q_value, signif) |>
    arrange(dataset, desc(test), factor)
  
  alpha.combined |> print(n=Inf)
  
  write.csv(alpha.combined, "tables/TableS5.alpha.diversity.tests.csv", row.names = FALSE)
}


################################################################################
# Alpha plots

{
  {
    theme_set(theme_bw())
    update_theme(legend.position = "none")
    update_theme(axis.ticks = element_blank())
    update_theme(panel.grid.major.x = element_blank())
    update_theme(strip.text = element_text(size = 12, hjust = 0))
    update_theme(strip.background = element_blank())
  }

  # Alpha diversity by host plant
  plot.alpha.host <- metadata_july |>
    pivot_longer(cols = all_of(alpha.columns)) |>
    mutate(name = alpha.labels[match(name, alpha.columns)]) |>
    ggplot(aes(x=host_plant, y=value, fill = host_plant)) +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.plant) +
    labs(x = NULL, y = "alpha diversity index", title = "Host Plant")

  plot.alpha.host

  table(metadata_july$host_plant)

  # Alpha diversity by month
  plot.alpha.month <- metadata_grt |>
    pivot_longer(cols = all_of(alpha.columns)) |>
    mutate(name = alpha.labels[match(name, alpha.columns)]) |>
    ggplot(aes(x=month, y=value, fill = month)) +
    facet_wrap(.~name, scales = "free") +
    geom_boxplot(outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1) +
    scale_fill_manual(values = color.by.month) +
    labs(x = NULL, y = "alpha diversity index", title = "Month")

  plot.alpha.month

  table(metadata_grt$month)

  # Alpha diversity by scrub_m2
  plot.alpha.scrub <- metadata_main |>
    ggplot(aes(x = scrub_m2, y = observed_ASVs)) +
    theme(panel.grid.major.x = NULL) + # Reset defaults
    geom_point(aes(color = scrub_m2), alpha = 0.8, size = 2) +
    scale_color_gradient(
      low = color.by.urban[1], high = color.by.urban[2]
    ) +
    geom_smooth(method = "lm", color = "grey35", se = TRUE) +
    labs(x = expression("Scrub (" * m^2 * ")"), y = "Observed ASVs", title = "Scrub land area")

  plot.alpha.scrub

  alpha.figure <- (
    plot.alpha.host | plot.alpha.month | plot.alpha.scrub
  ) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(widths = c(1, 1, 2/3))
  alpha.figure

  ggsave(
    filename = "plots/Fig3.alpha.diversity.png",
    plot = alpha.figure,
    width = 6.5, height = 2.5, scale = 2,
    dpi = 600, bg = "white"
  )

}

################################################################################
# Save the key objects
save(
  alpha.table, 
  alpha.columns, alpha.labels,
  categorical.columns, continuous.columns,
  alpha.wrst.results, alpha.lmer.results,
  alpha.wrst.results.long, alpha.lmer.results.long,
  alpha.combined,
  plot.alpha.host, plot.alpha.month, plot.alpha.scrub, 
  alpha.figure,
  file = "alpha.RData"
)
