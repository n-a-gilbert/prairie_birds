# script for a quick check of whether there is widespread support for 
# quadratic rather than linear effects of climate anomalies
# doing a simplified GLMM analysis for speed 

library(here)
library(tidyverse)
library(glmmTMB)

anomaly <- readr::read_csv(here::here("data/weather_anomaly.csv")) |> 
  dplyr::rename(point = siteID)

# read in bird data; contains columns for the focal covariates as well
d <- readr::read_csv(here::here("data/mn_prairie_bird_data_clean.csv")) |> 
  dplyr::filter(!is.na(site)) |> 
  dplyr::left_join(anomaly)

# retain species detected on at least 3% of surveys
sp_to_select <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::group_by(sp) |> 
  dplyr::summarise(prop = sum( n > 0) / sum(!is.na(n))) |> 
  dplyr::filter(prop >= 0.03) |> 
  # omit swallows (bars and tres) and ducks (mall) which are mostly flyovers, and unknown (unk) birds
  dplyr::filter(!sp %in% c("bars", "mall", "tres", "unk")) |> 
  dplyr::pull(sp)

# "z information" - aka information about species-point-year combinations 
z_info <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::filter(sp %in% sp_to_select) |> 
  dplyr::select(site, point, year, area, ratio, open, tanom_br, panom_br, sp) |> 
  # critical line here - we are narrowing down to point-species-year combinations (combos we estimate z, aka latent occurrence state, for)
  dplyr::distinct() |> 
  dplyr::group_by(sp) |> 
  # create species index
  dplyr::mutate(sp_z = dplyr::cur_group_id()) |> 
  dplyr::group_by(sp, site) |>
  # create species x site index
  dplyr::mutate(sp_site_z = dplyr::cur_group_id()) |> 
  dplyr::group_by(sp, point) |>
  # create species x point index
  dplyr::mutate(sp_point_z = dplyr::cur_group_id()) |> 
  dplyr::group_by(point, sp, year) |> 
  # create "z index" - for a grouping of species, point, and year
  dplyr::mutate(z_index = dplyr::cur_group_id()) |> 
  dplyr::arrange(z_index) |> 
  dplyr::ungroup()

# "y information" - aka the bird data and associated observation covariates
y_info <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::filter(sp %in% sp_to_select) |> 
  dplyr::select(site, point, area, ratio, open, tanom_br, panom_br, year, visit, date, obs, start, sp, n) |> 
  dplyr::group_by(sp) |> 
  # create species index
  dplyr::mutate(sp_y = dplyr::cur_group_id()) |> 
  dplyr::group_by(obs) |>
  # create observer index
  dplyr::mutate(obs_y = dplyr::cur_group_id()) |> 
  # link the bird data to the "z state"
  dplyr::full_join(
    z_info |> 
      dplyr::select(site, point, year, sp, z_index)) |> 
  dplyr::rename(z_index_y = z_index) |> 
  dplyr::mutate(y = ifelse(n > 0, 1, 0)) |> 
  dplyr::ungroup() 

# aggregate data for binomial GLMM ("x detections out of y surveys")
tmp <- y_info |> 
  dplyr::left_join(
    z_info |> 
      dplyr::select(site, point, year, sp, sp_site_z, sp_point_z)) |> 
  dplyr::group_by(site, point, area, ratio, open, tanom_br, panom_br, year, sp, sp_site_z, sp_point_z) |> 
  dplyr::summarise(y = sum(y), 
                   nvis = max(visit)) |> 
  dplyr::ungroup()

# fit models with linear and quadratic effects of precip to see which is preferred by
# AIC model selection

sp_list <- unique(tmp$sp)
res <- list(list())
for(i in 1:length(sp_list)){ # loop through and fit species-specific mdoels
  test <- tmp |> 
    filter(sp == sp_list[i]) |> 
    dplyr::mutate(psc = as.numeric(scale(panom_br)),
                  tan = as.numeric(scale(tanom_br)), 
                  open = as.numeric(scale(open)), 
                  ratio = as.numeric(scale(ratio)),
                  area = as.numeric(scale(log(area))), # log of reserve area
                  year = as.numeric(scale(year))) 
  
  # basically same model structure as the occ model
  linear <- glmmTMB::glmmTMB(
    cbind(y, nvis - y) ~ 1 + year + open + area + ratio + tan + psc +
    year:open + year:area + year:ratio + year:tan + year:psc +
    area:ratio +
    area:ratio:year + (1|site) + (1|point),
    family = binomial, 
    data = test)
  
  # now with quadratic effect of precip anomaly
  quadratic <- glmmTMB::glmmTMB(
    cbind(y, nvis - y) ~ 1 + year + open + area + ratio + tan + psc + I(psc^2) +
      year:open + year:area + year:ratio + year:tan + year:psc + year:I(psc^2) +
      area:ratio +
      area:ratio:year + (1|site) + (1|point),
    family = binomial, 
    data = test)
  
  res[[i]] <- AIC(linear, quadratic) |> 
    as_tibble(rownames = "model") |> 
    mutate(daic = AIC - min(AIC)) |> 
    add_column(sp = sp_list[i])
  
  print( paste0( round(100*(i/length(sp_list)),1), "% done" ))
}

bind_rows(res) |> 
  group_by(sp) |> 
  dplyr::mutate(max_daic = max(daic)) |> 
  dplyr::filter( daic == min(daic)) |> 
  # if models are within 2 delta AIC units of each other, no preference
  dplyr::mutate( pref = ifelse(max_daic < 2, "no pref", model)) |> 
  pull(pref) |> 
  table()


# same thing, but temp
sp_list <- unique(tmp$sp)
res_temp <- list(list())
for(i in 1:length(sp_list)){
  test <- tmp |> 
    filter(sp == sp_list[i]) |> 
    dplyr::mutate(psc = as.numeric(scale(panom_br)),
                  tan = as.numeric(scale(tanom_br)), 
                  open = as.numeric(scale(open)), 
                  ratio = as.numeric(scale(ratio)),
                  area = as.numeric(scale(log(area))), # log of reserve area
                  year = as.numeric(scale(year))) 
  
  linear <- glmmTMB::glmmTMB(
    cbind(y, nvis - y) ~ 1 + year + open + area + ratio + tan + psc +
      year:open + year:area + year:ratio + year:tan + year:psc +
      area:ratio +
      area:ratio:year + (1|site) + (1|point),
    family = binomial, 
    data = test)
  
  quadratic <- glmmTMB::glmmTMB(
    cbind(y, nvis - y) ~ 1 + year + open + area + ratio + tan + psc + I(tan^2) +
      year:open + year:area + year:ratio + year:tan + year:psc + year:I(tan^2) +
      area:ratio +
      area:ratio:year + (1|site) + (1|point),
    family = binomial, 
    data = test)
  
  res_temp[[i]] <- AIC(linear, quadratic) |> 
    as_tibble(rownames = "model") |> 
    mutate(daic = AIC - min(AIC)) |> 
    add_column(sp = sp_list[i])
  
  print( paste0( round(100*(i/length(sp_list)),1), "% done" ))
}

bind_rows(res_temp) |> 
  group_by(sp) |> 
  dplyr::mutate(max_daic = max(daic)) |> 
  dplyr::filter( daic == min(daic)) |> 
  dplyr::mutate( pref = ifelse(max_daic < 2, "no pref", model)) |> 
  pull(pref) |> 
  table()