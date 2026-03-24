library(here)
library(tidyverse)
library(readxl)
library(sf)
library(daymetr)

d <- readxl::read_excel( here::here("data/MN_prairie_bird_data.xlsx"))

# maybe should do May / June to capture current hygrid conditions

usa <- sf::st_as_sf(maps::map("state", fill = TRUE, plot = FALSE)) |> 
  dplyr::filter(ID == "minnesota")

site_coords <- d |> 
  dplyr::select(x = easting, y = northing) |> 
  dplyr::distinct() |> 
  sf::st_as_sf(coords = c("x", "y"),
               crs = 26915) |> 
  sf::st_transform(crs = 4326) 

ggplot() + 
  geom_sf(data = usa, aes(geometry = geom)) +
  geom_sf(data = site_coords, aes(geometry = geometry))

coords <- site_coords |> 
  sf::st_coordinates() |> 
  cbind(
    d |> 
      dplyr::select(easting, northing, site.point) |> 
      dplyr::distinct() |> 
      dplyr::select(site.point)
  )

daymetr::download_daymet( lat = coords[[1,2]], 
                          lon = coords[[1,1]], 
                          start = 1980, 
                          end = 2023)

df <- list(list())
for(i in 1:nrow(coords)){
  df[[i]] <- daymetr::download_daymet(
    # site = site_coords[[i, "site.point"]],
    lat = coords[[i, 2]],
    lon = coords[[i, 1]],
    start = 1980, 
    end = 2023,
    # start = site_dates[[1, "min_date"]],
    # end = site_dates[[1, "max_date"]],
    path = tempdir(),
    internal = TRUE, 
    silent = FALSE, 
    force = FALSE, 
    simplify = FALSE)$data |> 
    tibble::as_tibble() |>
    janitor::clean_names() |> 
    tibble::add_column( siteID = coords[[i, "site.point"]] ) |> 
    dplyr::mutate( date = as.Date( paste( yday, year ), format = "%j %Y")) |> 
    dplyr::select( siteID, year, yday, date, dayl_s:vp_pa )
  
  print(paste("finihsed", i, "of", nrow(coords)))
}

bind_rows(df)

setwd(here::here("data"))
write_csv(bind_rows(df), "site_weather_data.csv")

avg_25yr <- bind_rows(df) |> 
  dplyr::mutate(month = lubridate::month(date)) |> 
  dplyr::filter(month == 5 | month == 6) |> 
  dplyr::filter(year <= 2005) |> 
  dplyr::group_by(siteID, year) |> 
  dplyr::summarise(prcp = sum(prcp_mm_day),
                   tmp = mean(tmax_deg_c)) |>
  dplyr::group_by(siteID) |> 
  dplyr::summarise(mean_prcp = mean(prcp), 
                   sd_prcp = sd(prcp),
                   mean_t = mean(tmp))

avg_25yr_wi <- bind_rows(df) |> 
  dplyr::mutate(month = lubridate::month(date)) |> 
  dplyr::filter(month < 5) |> 
  dplyr::filter(year <= 2005) |> 
  dplyr::group_by(siteID, year) |> 
  dplyr::summarise(prcp = sum(prcp_mm_day),
                   tmp = mean(tmax_deg_c)) |>
  dplyr::group_by(siteID) |> 
  dplyr::summarise(mean_prcp = mean(prcp), 
                   sd_prcp = sd(prcp),
                   mean_t = mean(tmp))

anom <- bind_rows(df) |> 
  dplyr::mutate(month = lubridate::month(date)) |> 
  dplyr::filter(month == 5 | month == 6) |> 
  dplyr::filter(year >= 2008) |> 
  dplyr::group_by(siteID, year) |> 
  dplyr::summarise(prcp = sum(prcp_mm_day),
                   tmp = mean(tmax_deg_c)) |> 
  dplyr::full_join(avg_25yr) |> 
  dplyr::mutate(anom = prcp - mean_prcp,
                tanom = tmp - mean_t) 

anom_wi <- bind_rows(df) |> 
  dplyr::mutate(month = lubridate::month(date)) |> 
  dplyr::filter(month < 5) |> 
  dplyr::filter(year >= 2008) |> 
  dplyr::group_by(siteID, year) |> 
  dplyr::summarise(prcp = sum(prcp_mm_day),
                   tmp = mean(tmax_deg_c)) |> 
  dplyr::full_join(avg_25yr_wi) |> 
  dplyr::mutate(anom = prcp - mean_prcp,
                tanom = tmp - mean_t) 

climate <- anom |> 
  dplyr::ungroup() |> 
  dplyr::select(siteID, year, panom_br = anom, tanom_br = tanom) |> 
  dplyr::full_join(
    anom_wi |> 
      dplyr::ungroup() |> 
      dplyr::select(siteID, year, panom_nb = anom, tanom_nb = tanom))

write_csv(climate, "weather_anomaly.csv")  
