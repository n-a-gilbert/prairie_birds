# create figure 3
# visualize community-level and species-level covariate effects and influence on trend
library(here)
library(tidyverse)
library(patchwork)
library(MCMCvis)
library(MetBrewer)
library(reshape2)

load(here::here("results/prairie_bird_trends_global_tanom_panom_br22026-02-01.RData"))

sp_key <- readr::read_csv(here::here("data/sp_key.csv"))

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

# create design matrix for covariates for occurrence
x <- z_info |> 
  dplyr::select(open, area, ratio, tanom_br, panom_br, year) |> 
  tibble::add_column(int = 1) |>  # column of repeated 1s for the intercept
  dplyr::mutate(tanom_br = as.numeric(scale(tanom_br)),
                panom_br = as.numeric(scale(panom_br)),
                open = as.numeric(scale(open)), 
                ratio = as.numeric(scale(ratio)),
                area = as.numeric(scale(log(area))), # log of reserve area
                year = as.numeric(scale(year))) |>
  # create columns for interaction terms
  dplyr::mutate( 
    tanom_br2 = tanom_br * tanom_br,
    panom_br2 = panom_br * tanom_br,
    area.ratio = area * ratio, 
    tanom_br.yr = tanom_br * year,
    tanom_br2.yr = tanom_br2 * year,
    panom_br.yr = panom_br * year,
    panom_br2.yr = panom_br2*year,
    open.yr = open * year, 
    ratio.yr = ratio * year, 
    area.yr = area * year, 
    area.ratio.yr = area * ratio * year ) |> 
  dplyr::select( int, year, tanom_br, tanom_br2, panom_br, panom_br2, area, open, ratio, 
                 area.ratio, 
                 tanom_br.yr, tanom_br2.yr, panom_br.yr, panom_br2.yr, area.yr, open.yr, ratio.yr, 
                 area.ratio.yr ) |> 
  as.matrix() |> 
  unname()

# create design matrix for covariates for detection
w <- y_info |> 
  tibble::add_column(int = 1) |> # column of repeated 1s for the intercept
  dplyr::select(int, date, start) |> 
  dplyr::mutate(date = as.numeric(scale(date)),
                start = as.numeric(scale(start))) |> 
  dplyr::mutate(date2 = date*date, 
                start2 = start*start) |> 
  as.matrix() |> 
  unname() 

area.sc <- scale(log(z_info$area))
ratio.sc <- scale(z_info$ratio)
open.sc <- scale(z_info$open)
tanom.sc <- scale(z_info$tanom_br)
panom.sc <- scale(z_info$panom_br)
year.sc <- scale(z_info$year)

# int, year, tanom_nb, panom_nb, area, open, ratio, 
# area.ratio, 
# tanom_nb.yr, panom_nb.yr, area.yr, open.yr, ratio.yr, 
# area.ratio.yr

param_key <- tibble(
  param = 1:18, 
  name = c("intercept", "yr", "tanom_br", "tanom_br2",
           "panom_br", "panom_br2",
           "area", "open", "ratio", 
           "area.ratio", 
           "tanom_br.yr", "tanom_br2.yr",
           "panom_br.yr", "panom_br2.yr",
           "area.yr", "open.yr", "ratio.yr", 
           "area.ratio.yr"))

post <- MCMCvis::MCMCpstr( out, params = "mu_beta", type = "chains")[[1]] |> 
  tibble::as_tibble(rownames = "param") |> 
  tidyr::pivot_longer(dplyr::starts_with("V"), names_to = "iter", values_to = "value") |> 
  dplyr::mutate( across( param:iter, function(x) readr::parse_number( x ))) |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

stars <- post |>
  dplyr::select(iter, yr:area.ratio.yr) |> 
  tidyr::pivot_longer(yr:area.ratio.yr, names_to = "var", values_to = "val") |> 
  dplyr::group_by(var) |> 
  dplyr::summarise( 
    mean = mean(val),
    l95 = quantile(val, c(0.025)), 
    u95 = quantile(val, c(0.975)), 
    l68 = quantile(val, c(0.160)), 
    u68 = quantile(val, c(0.840))) |> 
  mutate(sig68 = ifelse(l68 < 0 & u68 > 0, 0, 1), 
         sig95 = ifelse(l95 < 0 & u95 > 0, 0, 1)) |> 
  mutate(lab = ifelse( sig95 == 1, "**", 
                       ifelse(sig95 == 0 & sig68 == 1, "*", ""))) |> 
  dplyr::select(var, lab) 

p_open <- tidyr::expand_grid(
  yr = c(2008, 2023),
  open = seq(from = min(open.sc), to = max(open.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

p_area <- tidyr::expand_grid(
  yr = c(2008, 2023),
  area = seq(from = min(area.sc), to = max(area.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

p_ratio <- tidyr::expand_grid(
  yr = c(2008, 2023),
  ratio = seq(from = min(ratio.sc), to = max(ratio.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

p_tanom <- tidyr::expand_grid(
  yr = c(2008, 2023),
  tanom = seq(from = min(tanom.sc), to = max(tanom.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

p_panom <- tidyr::expand_grid(
  yr = c(2008, 2023),
  panom = seq(from = min(panom.sc), to = max(panom.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

post <- MCMCvis::MCMCpstr( out, params = "mu_beta", type = "chains")[[1]] |> 
  tibble::as_tibble(rownames = "param") |> 
  tidyr::pivot_longer(dplyr::starts_with("V"), names_to = "iter", values_to = "value") |> 
  dplyr::mutate( across( param:iter, function(x) readr::parse_number( x ))) |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

pdat <- tibble::tibble(
  open = seq(from = min(open.sc), to = max(open.sc), length.out = 10), 
  area = seq(from = min(area.sc), to = max(area.sc), length.out = 10), 
  ratio = seq(from = min(ratio.sc), to = max(ratio.sc), length.out = 10), 
  tanom = seq(from = min(tanom.sc), to = max(tanom.sc), length.out = 10),
  panom = seq(from = min(panom.sc), to = max(panom.sc), length.out = 10))

fit_open <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (open.x * open.y))) |> 
  dplyr::select(iter, open = open.y, psi ) |> 
  dplyr::group_by( open) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( open = open*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center")) |> 
  dplyr::mutate(open = round(open, 1)) |> 
  dplyr::rename(cov = open) |> 
  tibble::add_column( cov_name = "Percent open habitat")

fit_area <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (area.x * area.y))) |> 
  dplyr::select(iter, area = area.y, psi ) |> 
  dplyr::group_by(area) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( area = area*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center")) |> 
  dplyr::mutate(area = round(area, 1)) |> 
  dplyr::rename(cov = area) |> 
  dplyr::mutate(cov = exp(cov)) |> 
  tibble::add_column( cov_name = "Reserve area (km2)")

fit_ratio <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (ratio.x * ratio.y))) |> 
  dplyr::select(iter, ratio = ratio.y, psi ) |> 
  dplyr::group_by( ratio) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( ratio = ratio*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center")) |> 
  # dplyr::mutate(ratio = round(ratio, 1)) |> 
  dplyr::rename(cov = ratio) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio")

fit_tanom <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (tanom_br * tanom) + (tanom_br2 * tanom * tanom))) |> 
  dplyr::select(iter, tanom = tanom, psi ) |> 
  dplyr::group_by(tanom) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( tanom = tanom*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |> 
  dplyr::rename(cov = tanom) |> 
  tibble::add_column( cov_name = "Temp. anomaly (°C)")

fit_panom <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (panom_br * panom) + (panom_br2 * panom * panom))) |> 
  dplyr::select(iter, panom = panom, psi ) |> 
  dplyr::group_by(panom) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( panom = panom*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |> 
  dplyr::rename(cov = panom) |> 
  tibble::add_column( cov_name = "Precip. anomaly (mm)")

com_fit <- dplyr::full_join(
  fit_area, 
  fit_open) |> 
  dplyr::full_join(fit_ratio) |> 
  dplyr::full_join(fit_tanom) |> 
  dplyr::full_join(fit_panom) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Temp. anomaly (°C)",
                                     "Precip. anomaly (mm)")))

post_sp <- MCMCvis::MCMCpstr( out, params = "beta", type = "chains")[[1]] |>
  reshape2::melt(c("sp", "param", "iter")) |> 
  tibble::as_tibble() |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

fit_open_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (open.x * open.y))) |> 
  dplyr::select(sp, iter, open = open.y, psi ) |> 
  dplyr::group_by( sp, open) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( open = open*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center")) |> 
  dplyr::rename(cov = open) |> 
  tibble::add_column( cov_name = "Percent open habitat")

fit_area_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (area.x * area.y))) |> 
  dplyr::select(sp, iter, area = area.y, psi ) |> 
  dplyr::group_by( sp, area) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( area = area*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center")) |> 
  dplyr::rename(cov = area) |> 
  dplyr::mutate(cov = exp(cov)) |> 
  tibble::add_column( cov_name = "Reserve area (km2)")

fit_ratio_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (ratio.x * ratio.y))) |> 
  dplyr::select(sp, iter, ratio = ratio.y, psi ) |> 
  dplyr::group_by( sp, ratio) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( ratio = ratio*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center")) |> 
  dplyr::rename(cov = ratio) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio")

fit_tanom_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (tanom_br * tanom + (tanom_br2 * tanom * tanom)))) |> 
  dplyr::select(sp, iter, tanom = tanom, psi ) |> 
  dplyr::group_by( sp, tanom) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( tanom = tanom*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |> 
  dplyr::rename(cov = tanom) |> 
  tibble::add_column( cov_name = "Temp. anomaly (°C)")

fit_panom_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(psi = plogis( intercept +  (panom_br * panom) + (panom_br2 * panom * panom))) |> 
  dplyr::select(sp, iter, panom = panom, psi ) |> 
  dplyr::group_by( sp, panom) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( panom = panom*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |> 
  dplyr::rename(cov = panom) |> 
  tibble::add_column( cov_name = "Precip. anomaly (mm)")

sp_fit <- dplyr::full_join(fit_area_sp, fit_open_sp) |> 
  dplyr::full_join(fit_ratio_sp) |> 
  dplyr::full_join(fit_tanom_sp) |> 
  dplyr::full_join(fit_panom_sp) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Temp. anomaly (°C)",
                                     "Precip. anomaly (mm)"))) |> 
  dplyr::left_join(sp_key)

xpos <- com_fit |> 
  dplyr::group_by(cov_name) |> 
  dplyr::filter(cov == min(cov)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(xpos = ifelse(cov == 0, 0.05, cov + (0.05*abs(cov)) )) |> 
  dplyr::select(cov_name, cov = xpos) |> 
  tibble::add_column(mean = 1.05)

dict <- c(
  area = "Reserve area (km2)", 
  open = "Percent open habitat",
  ratio = "Perimeter-to-area ratio",
  tanom_br = "Temp. anomaly (°C)",
  panom_br = "Precip. anomaly (mm)"
)

com_stars <- stars |> 
  dplyr::filter( var %in% c("area", "open", "ratio", "tanom_br", "panom_br")) |>
  dplyr::mutate( cov_name = dict[var]) |> 
  # tibble::add_column(cov_name = c("Precip. anomaly (mm)",
  # "Reserve area (km2)", "Percent open habitat", "Perimeter-to-area ratio")) |> 
  dplyr::left_join(xpos) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Temp. anomaly (°C)",
                                     "Precip. anomaly (mm)")))

# marginal effects of temp and precip - assemblage level, with 
# species mean in the background
( panel_a <- ggplot() + 
    scale_y_continuous(limits = c(0, 1.1), 
                       breaks = c(0, 0.5, 1)) +
    coord_cartesian(clip = "off") +
    geom_line(data = filter(sp_fit,
                            cov_name %in% 
                              c("Temp. anomaly (°C)",
                                "Precip. anomaly (mm)")),
              aes(x = cov, 
                  y = mean,
                  group = factor(code)),
              alpha = 0.2,
              linewidth = 0.4,
              color = MetBrewer::MetPalettes$Hiroshige[[1]][10]) +
    facet_wrap(~cov_name, scales = "free_x", nrow = 1)  +
    geom_ribbon( data = filter(com_fit,
                               cov_name %in% 
                                 c("Temp. anomaly (°C)",
                                   "Precip. anomaly (mm)")),
                 aes(x = cov, ymin = l95, ymax = u95),
                 color = NA,
                 fill = "#a62517",
                 alpha = 0.3) +
    geom_ribbon( data = filter(com_fit,
                               cov_name %in% 
                                 c("Temp. anomaly (°C)",
                                   "Precip. anomaly (mm)")),
                 aes(x = cov, ymin = l68, ymax = u68),
                 color = NA,
                 fill = "#a62517",
                 alpha = 0.4) +
    geom_line(data = filter(com_fit,
                            cov_name %in% 
                              c("Temp. anomaly (°C)",
                                "Precip. anomaly (mm)")), 
              aes(x = cov, y = mean),
              color = "#a62517",
              linewidth = 1) +
    theme_minimal() + 
    scale_color_viridis_d() +
    labs(x = "Covariate value", 
         y = "Occurrence probability") +
    ggplot2::theme(axis.line = element_line( linewidth = 0.15, color = "black"),
                   panel.grid = element_line(linewidth = 0.1, color = "gray80"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 8, color = "black"),
                   strip.text = element_text(size = 8, color = "black", face = "bold"), 
                   axis.text = element_text(size = 7, color = "black"),
                   panel.background = element_rect(color = NA, fill = "white"), 
                   plot.background = element_rect(color = NA, fill = "white"),
                   plot.title = element_text(color = "black", size = 10),
                   legend.position = "none"))

# species-level precip marginal effects with species following
# 'just right' hypothesis highlighted
( panom_p <- ggplot() +
    geom_line(data = filter(sp_fit,
                            cov_name == "Precip. anomaly (mm)" &
                              !(code %in% c("weme", "lesp", "savs"))),
              aes(x = cov, y = mean, group = code),
              alpha = 0.2,
              linewidth = 0.4,
              color = MetBrewer::MetPalettes$Hiroshige[[1]][10]) +
    geom_line(data = filter(sp_fit,
                            cov_name == "Precip. anomaly (mm)" &
                              (code %in% c("weme", "lesp", "savs"))),
              aes(x = cov, y = mean, group = code),
              alpha = 1,
              color = MetBrewer::MetPalettes$Isfahan1[[1]][3],
              linewidth = 1.3) +
    theme_minimal() + 
    scale_color_viridis_d() +
    labs(x = "Covariate value", 
         y = "Occurrence probability") +
    ggplot2::theme(axis.line = element_line( linewidth = 0.15, color = "black"),
                   panel.grid = element_line(linewidth = 0.1, color = "gray80"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 8, color = "black"),
                   strip.text = element_text(size = 8, color = "black", face = "bold"), 
                   axis.text = element_text(size = 7, color = "black"),
                   panel.background = element_rect(color = NA, fill = "white"), 
                   plot.background = element_rect(color = NA, fill = "white"),
                   plot.title = element_text(color = "black", size = 10),
                   legend.position = "none"))

# temperature - no species followed the 'just right' hypothesis
( tanom_p <- ggplot() +
    geom_line(data = filter(sp_fit,
                            cov_name == "Temp. anomaly (°C)"),
              aes(x = cov, y = mean, group = code),
              alpha = 0.2,
              linewidth = 0.4,
              color = MetBrewer::MetPalettes$Hiroshige[[1]][10]) +
   theme_minimal() + 
    scale_color_viridis_d() +
    labs(x = "Covariate value", 
         y = "Occurrence probability") +
    ggplot2::theme(axis.line = element_line( linewidth = 0.15, color = "black"),
                   panel.grid = element_line(linewidth = 0.1, color = "gray80"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 8, color = "black"),
                   strip.text = element_text(size = 8, color = "black", face = "bold"), 
                   axis.text = element_text(size = 7, color = "black"),
                   panel.background = element_rect(color = NA, fill = "white"), 
                   plot.background = element_rect(color = NA, fill = "white"),
                   plot.title = element_text(color = "black", size = 10),
                   legend.position = "none"))


panel_a / ( tanom_p | panom_p) +
  plot_layout(guides = "collect")

setwd(here::here("figures"))
ggsave(
  "figure_s01.png", 
  width = 4, 
  height = 3, 
  units = "in", 
  dpi = 600)

# here's how I evaluated whether species showed support for 'just right' hypothesis
# looped through and plotted marginal effects for each species
# if the curve was a "smile", that is opposite what the JRH predicts
# if the curve was a "frown" (hump shaped), that is what the JRH predicts
# did this for temp and precip and scrolled through each species
codes <- unique(sp_fit$code)
for( i in 1:length(unique(sp_fit$code))) {
  
  dat <- sp_fit |> 
    dplyr::filter(code == codes[i]) |> 
    dplyr::filter(cov_name == "Temp. anomaly (°C)")
  
  p <- ggplot(dat, aes(x = cov, y = mean)) +
    geom_line(linewidth = 3) +
    labs(title = codes[i])
  
  print(p)
  
}