################################################################################
# Summary figure for robust DA ASVs using the ALDEx2 results

{
  library(phyloseq)
  library(tidyverse)
}

load("Jhae.microbiome.RData", verbose = TRUE)
load("Jhae.microbiome.aldex.RData", verbose = TRUE)

# Manually curated list of robust DA ASVs (those DA by 2 or more methods)
{
  robust.ASVs.by.host <- c(
    "3cc27a22b5adf4e21f48c1f09ee198de"
  )
  robust.ASVs.by.month <- c(
    "8084aa68823c09679d87c31776e22478", "fcb2943f9c7b866112a146c62e824bdb"
  )
  
  ancom.ASVs.by.host <- c(
    "7364a9d6bfb4a993fabfa2e22252a18d", "02502141f0a1bf07cc6d6c60b94433f5", "675a0962f1197511dde4ff1172f15cf5", 
    "78849fe4c1e26973c2a70cf0727da7e8"
  )
  ancom.ASVs.by.month <- c(
    "4bca6d2c2632165857073fbaf8e36f2c", "67fabb997fa0ee06acd42351d5489c0d", "08bb6389255945886b20e36acbab566f", 
    "aa5738a607473674763d7e712002a471", "2d4ebd9a510325bf2ec849e5de14bf09", "dab8a32307bba56fc64cb9fe24c668b6", 
    "7679539daee22c8c08343bf10bdb96c4", "646c9d0d52e1d605a7124052de9442d7", "d08be062b768afab4c0f1a77285e454d", 
    "b4b952684ab336effd02bf62380fdbcc", "e381bce7cb0b2c24e9957a72073bc457", "419af7300fc62bd827b4de4172d5ab89", 
    "461cf6a54476ed11f3feb7289cfeb2e5", "e2cae11fb0bdbbe746db5f6861af0b7d", "b5071e583d3e566ffcc3edd393ecf184", 
    "d243b2ed42daeb914fe8937b05198369", "10bfe46c36777ac10d115e9d41b740c5"
  )
  
  maaslin.ASVs.by.host <- c(
    "2dc455842240c61e9ef032e9188da70b", "41b2d4598322bf23bfe6fea245d1cc2b", "4a76266d6ae2d1c237bfcd773e7c76bf", 
    "5a2e07879cad183af5484d4a966f3dd9", "606a76fee20d4a9914edefcd5499234c", "789c66897d76d452e204430813721965", 
    "d562e3cded5f63c9a0cd3ea75ff18f77", "fca75566d53f92ea06d1e663876983e0", "4bf4f2cd62fecb0557ee4275fa2c9011", 
    "6431099aca94ed96093977b221b67a2a", "668c2e39fa41c78ac6934f9da7dd8f7e", "8e58415bd25490f244500fc718b657b3", 
    "908674b8f5b0e8bd8efc18fcf1767564", "9bbddedd973b089f52cec3e9373aa576", "cbf3919cdcc89b36c4e339775e22897a"
  )
  maaslin.ASVs.by.month <- c(
    "fcb2943f9c7b866112a146c62e824bdb", "d54e70632fabd5c40f20787ecffbf0b4", "8084aa68823c09679d87c31776e22478"
  )
  
  aldex.ASVs.by.host <- aldex_host |> 
    filter(wi.eBH < 0.05) |> 
    rownames_to_column(var = "ASV_ID") |> 
    pull(ASV_ID)
  
  aldex.ASVs.by.month <- aldex_month |> 
    filter(wi.eBH < 0.05) |> 
    rownames_to_column(var = "ASV_ID") |> 
    pull(ASV_ID)

  notable.ASVs.by.host <- unique(c(ancom.ASVs.by.host, aldex.ASVs.by.host, maaslin.ASVs.by.host))
  notable.ASVs.by.month <- unique(c(ancom.ASVs.by.month, aldex.ASVs.by.month, maaslin.ASVs.by.month))
}

{
  theme_set(theme_bw())
  update_theme(legend.position = "bottom")
  update_theme(axis.ticks = element_blank())
}

taxa_july <- as.data.frame(as(tax_table(ps_july), "matrix"))

DA_host_plot_df <- aldex_host |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_july |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  mutate(Genus = str_replace(Genus, "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium complex")) |> 
  mutate(
    Significance = case_when(
      (ASV_ID %in% robust.ASVs.by.host) ~ "Robust DA",
      (ASV_ID %in% aldex.ASVs.by.host) ~ "ALDEx2",
      (ASV_ID %in% maaslin.ASVs.by.host) ~ "MaAsLin3",
      (ASV_ID %in% ancom.ASVs.by.host) ~ "ANCOM-BC2",
      TRUE ~ " Not Significant"
    ),
    Taxa_Label = case_when(
      (ASV_ID %in% notable.ASVs.by.host)  & !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.host)  & !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.host)  & !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.host)  & !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.host)  & !is.na(Phylum) & Phylum != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.host)  ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  ) 

plot_DA_host <- DA_host_plot_df |>
  ggplot(aes(x = effect, y = -log10(wi.ep), color = Significance)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_text_repel(
    aes(label = Taxa_Label),
    size = 2, color = "black",
    box.padding = 0.5, point.padding = 0.3,
    force = 4, max.overlaps = Inf,
    show.legend = FALSE) +
  scale_color_manual(values = c(
    "Robust DA" = "#d95f02",
    "ALDEx2" = "#31688E",
    "ANCOM-BC2" = "#22A884",
    "MaAsLin3" = "#A5DB36",
    " Not Significant" = "grey50"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  expand_limits(x = c(-6, 6), y = c(0, 5.5)) +
  labs(
    title = "Host plant",
    x = "ALDEx2 Effect Size",
    y = expression(-log[10]("p-value"))
  )

plot_DA_host


taxa_grt <- as.data.frame(as(tax_table(ps_grt), "matrix"))

DA_month_plot_df <- aldex_month |>
  rownames_to_column(var = "ASV_ID") |>
  left_join(taxa_grt |> rownames_to_column(var = "ASV_ID"), by = "ASV_ID") |>
  mutate(Genus = str_replace(Genus, "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium complex")) |> 
  mutate(
    Significance = case_when(
      (ASV_ID %in% robust.ASVs.by.month) ~ "Robust DA",
      (ASV_ID %in% aldex.ASVs.by.month) ~ "ALDEx2",
      (ASV_ID %in% maaslin.ASVs.by.month) ~ "MaAsLin3",
      (ASV_ID %in% ancom.ASVs.by.month) ~ "ANCOM-BC2",
      TRUE ~ " Not Significant"
    ),
    Taxa_Label = case_when(
      (ASV_ID %in% notable.ASVs.by.month)  & !is.na(Genus) & Genus != "NA" ~ paste0(Genus, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.month)  & !is.na(Family) & Family != "NA" ~ paste0(Family, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.month)  & !is.na(Order) & Order != "NA" ~ paste0(Order, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.month)  & !is.na(Class) & Class != "NA" ~ paste0(Class, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.month)  & !is.na(Phylum) & Phylum != "NA" ~ paste0(Phylum, "\n(", substr(ASV_ID, 1, uuid.abbrev), ")"),
      (ASV_ID %in% notable.ASVs.by.month)  ~ paste0("Unclassified\n(", substr(ASV_ID, 1, uuid.abbrev), ")")
    )
  ) 

plot_DA_month <- DA_month_plot_df |>
  ggplot(aes(x = effect, y = -log10(wi.ep), color = Significance)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_text_repel(
    aes(label = Taxa_Label),
    size = 2, color = "black",
    box.padding = 0.5, point.padding = 0.3,
    force = 4, max.overlaps = Inf,
    show.legend = FALSE) +
  scale_color_manual(values = c(
    "Robust DA" = "#d95f02",
    "ALDEx2" = "#31688E",
    "ANCOM-BC2" = "#22A884",
    "MaAsLin3" = "#A5DB36",
    " Not Significant" = "grey50"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  expand_limits(x = c(-6, 6), y = c(0, 5.5)) +
  labs(
    title = "Month",
    x = "ALDEx2 Effect Size",
    y = expression(-log[10]("p-value"))
  ) +
  theme(legend.position = "none")

plot_DA_month

DA_figure <- (
  plot_DA_host | plot_DA_month
) +
  plot_annotation(tag_levels = 'A')
DA_figure

set.seed(123)
ggsave(
  filename = "plots/Fig5.consensus.DA.png",
  plot = DA_figure,
  width = 6.5, height = 3.0, scale = 1.65,
  dpi = 600, bg = "white"
)

