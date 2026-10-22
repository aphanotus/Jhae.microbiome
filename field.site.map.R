################################################################################
# Map of field sites

{
  library(tidyverse)
  library(sf)
  library(tigris)
}

florida_boundary <- states(cb = TRUE, resolution = "20m") |> 
  filter(NAME == "Florida")

sites <- read.csv("qiime/curated.metadata.csv", header = TRUE) |> 
  filter(!is.na(site_name) & site_name != "") |> 
  filter(str_detect(city, "FL$")) |> 
  mutate(
    site_name = str_replace(site_name, "Delrey", "Delray"),
    host_plant = str_replace(host_plant, "BV", "C. corindum"),
    host_plant = str_replace(host_plant, "GRT", "K. paniculata"),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) |> 
  group_by(site_name) |> 
  summarise(
    latitude = mean(latitude, na.rm = T), 
    longitude = mean(longitude, na.rm = T),
    host_plant = unique(host_plant)
  ) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

sites

fl_cities <- tibble(
  location = c("Miami", "Orlando", "Tampa", "Jacksonville", "Fort Myers","Pensacola"),
  longitude = c(-80.1918, -81.3792, -82.4572, -81.6557, -81.8556, -87.2115), 
  latitude = c(25.7617, 28.5383, 27.9506, 30.3322, 26.6233, 30.4740)     
) |>
  st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = 4326
  )

plot.field.sites <- ggplot() +
  geom_sf(data = florida_boundary, fill = "grey80", color = "white") +
  geom_sf(data = fl_cities, color = "grey40", size = 3) +
  geom_sf(aes(color = host_plant), data = sites, size = 2, alpha = 0.8) +
  scale_color_manual(name = "Host Plant", values = c("#9acd32","#DAA520")) +
  geom_sf_text(
    data = fl_cities, 
    aes(label = location), 
    nudge_y = 0.2,       # Nudges the text slightly north (in degrees)
    fontface = "bold",   # Makes the text pop
    color = "grey40"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold.italic", size = 12),
    legend.position = c(0.25, 0.25)
  )

plot.field.sites

ggsave(
  filename = "plots/FigX.field.sites.png", 
  plot = plot.field.sites,
  width = 4.5, height = 4.5, scale = 1, 
  dpi = 300, bg = "white"
)

################################################################################
# Save the key objects
save(
  sites, florida_boundary, fl_cities,
  plot.field.sites,
  file = "field.site.map.RData"
)


