library(sf)
library(here)
library(tidyverse)
library(readxl)

grass <- sf::st_read( here::here("data/spatial/SPICE_contig_grass_NoHoles_20260213.shp"))

d <- readxl::read_excel( here::here("data/MN_prairie_bird_data.xlsx"))

site_coords <- d |> 
  dplyr::select(site.point, x = easting, y = northing) |> 
  dplyr::distinct() |> 
  sf::st_as_sf(coords = c("x", "y"),
               crs = 26915) 

grass_diss <- grass |> 
  sf::st_make_valid() |> 
  dplyr::summarise()

gr_boundary <- sf::st_boundary( grass_diss )

buff_250 <- sf::st_buffer( site_coords, dist = 250 )
buff_3000 <- sf::st_buffer( site_coords, dist = 3000 )

buff_250$buffer_area <- sf::st_area( buff_250 )
buff_3000$buffer_area <- sf::st_area( buff_3000 )

grass_buff_250 <- sf::st_intersection( buff_250, grass_diss)
grass_buff_3000 <- sf::st_intersection( buff_3000, grass_diss)

grass_buff_250$grass_250 <- sf::st_area( grass_buff_250 )
grass_buff_3000$grass_3000 <- sf::st_area( grass_buff_3000 )

edge_250 <- sf::st_intersection( buff_250, gr_boundary )
edge_3000 <- sf::st_intersection( buff_3000, gr_boundary )

edge_250$el_250 <- sf::st_length(edge_250)
edge_3000$el_3000 <- sf::st_length(edge_3000)

metrics250 <- buff_250 |> 
  dplyr::left_join( sf::st_drop_geometry(grass_buff_250)) |> 
  dplyr::mutate(grass_250 = ifelse(is.na(grass_250), 0, grass_250)) |> 
  dplyr::mutate(pgrass_250 = as.numeric(grass_250 / buffer_area)) |> 
  dplyr::left_join( sf::st_drop_geometry(edge_250)) |>
  dplyr::mutate(el_250 = ifelse(is.na(el_250), 0, el_250)) |> 
  dplyr::mutate( ed_250 = el_250 / grass_250 ) |> 
  sf::st_drop_geometry() |> 
  dplyr::select(site.point, pgrass_250, ed_250)

metrics3000 <- buff_3000 |> 
  dplyr::left_join( sf::st_drop_geometry(grass_buff_3000)) |> 
  dplyr::mutate(grass_3000 = ifelse(is.na(grass_3000), 0, grass_3000)) |> 
  dplyr::mutate(pgrass_3000 = as.numeric(grass_3000 / buffer_area)) |> 
  dplyr::left_join( sf::st_drop_geometry(edge_3000)) |>
  dplyr::mutate(el_3000 = as.numeric( ifelse(is.na(el_3000), 0, el_3000))) |> 
  dplyr::mutate( ed_3000 = el_3000 / grass_3000 ) |> 
  sf::st_drop_geometry() |> 
  dplyr::select(site.point, pgrass_3000, ed_3000) 

vars <- dplyr::left_join(metrics250, metrics3000) |> 
  dplyr::mutate(across(c(ed_250, ed_3000), function(x) ifelse(is.na(x), 0, x))) |> 
  dplyr::rename(point = site.point)

write_csv(vars, here::here("data/habitat_vars.csv"))

cleandat <- readr::read_csv(here::here("data/mn_prairie_bird_data_clean.csv"))

cleandat |> 
  left_join(vars) 

