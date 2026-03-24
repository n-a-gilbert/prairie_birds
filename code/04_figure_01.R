# plot a map showing the locations of the surveyed reserves in MN with an inset showing location of MN in US
library(tidyverse)
library(sf)
library(cowplot)
library(MetBrewer)
library(here)
library(maps)
library(ggspatial)

sites <- sf::st_read(here::here("data/spatial/allsites_SPICE_20201228_multipart.shp")) |> 
  dplyr::filter(!grepl("dropped", SiteName)) |> 
  sf::st_centroid()

prov <- sf::st_read(here::here("data/spatial/ecs_provinces_of_mn_v99a.shp")) |> 
  sf::st_transform( sf::st_crs(sites) ) |> 
  dplyr::mutate(PROVNAME = stringr::str_remove_all(PROVNAME, "Province"))

usa <- sf::st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) 

mn <- dplyr::filter(usa, ID == "minnesota") |> 
  sf::st_transform( st_crs(sites))

( pts <- ggplot() + 
    geom_sf(data = prov, aes(geometry = geometry, fill = PROVNAME), color = NA) +
    geom_sf(data = sites, aes(geometry = geometry), size = 2.5, pch = 21,
            fill = "black", color = "white") +
    scale_fill_manual("Province", 
                      values = MetBrewer::MetPalettes$Isfahan1[[1]][c(5,6,2,3)]) +
    theme_void() +
    theme(legend.title = element_text(size = 9, color = "black"), 
          legend.text = element_text(size = 8, color = "black"),
          plot.margin = margin(0, 10, 10, 0, unit = "pt"),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA), 
          panel.background = element_rect(fill = "white", color = NA)) +
    guides(fill = guide_legend( nrow = 2, title.position = "top", title.hjust = 0.5)) +
    ggspatial::annotation_scale(location = "tr"))

inset <- ggplot() +
  geom_sf(data = usa) +
  geom_sf(data = filter(usa, ID == "minnesota"), fill = "gray20") +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        linewidth = 0.1)

( figure_01 <- cowplot::ggdraw(pts) +
    cowplot::draw_plot(
      {
        inset
      },
      x = 0.6,
      y = 0.35, 
      width = 0.25,
      height = 0.25
    )
)

setwd(here::here("figures"))
ggsave(
  filename = "figure_01.png",
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)