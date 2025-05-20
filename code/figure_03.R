# create figure 3
# visualize community-level and species-level covariate effects and influence on trend
library(here)
library(tidyverse)
library(patchwork)
library(MCMCvis)
library(MetBrewer)
library(reshape2)

load(here::here("results/prairie_bird_trends_global2025-02-28.RData"))

sp_key <- readr::read_csv(here::here("data/sp_key.csv"))

# format data as we did in analysis to get unscaled covariate values
# read in bird data; contains columns for the focal covariates as well
d <- readr::read_csv(here::here("data/mn_prairie_bird_data_clean.csv"))

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
  dplyr::select(site, point, year, area, ratio, open, anom, sp) |> 
  # critical line here - we are narrowing down to point-species-year combiations (combos we estimate z, aka latent occurence state, for)
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

area.sc <- scale(log(z_info$area))
ratio.sc <- scale(z_info$ratio)
open.sc <- scale(z_info$open)
anom.sc <- scale(z_info$anom)
year.sc <- scale(z_info$year)

param_key <- tibble(
  param = 1:12, 
  name = c("intercept", "yr", "anom", "area", "open", "ratio", 
           "area.ratio", 
           "anom.yr", "area.yr", "open.yr", "ratio.yr", 
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
  dplyr::summarise( l95 = quantile(val, c(0.025)), 
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

p_anom <- tidyr::expand_grid(
  yr = c(2008, 2023),
  anom = seq(from = min(anom.sc), to = max(anom.sc), length.out = 20)) |> 
  dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))

param_key <- tibble(
  param = 1:12, 
  name = c("intercept", "yr", "anom", "area", "open", "ratio", 
           "area.ratio", 
           "anom.yr", "area.yr", "open.yr", "ratio.yr", 
           "area.ratio.yr"))

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
  anom = seq(from = min(anom.sc), to = max(anom.sc), length.out = 10))

fit_open <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr * 0 +                      # 2
                      anom.x * 0 +                    # 3
                      area.x * 0 +                    # 4
                      open.x * open.y +             # 5 
                      ratio.x * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * 0 * 0 +             # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0
    )) |> 
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
  dplyr::mutate(
    psi = plogis(     intercept +                   # 1
                        yr * 0 +                      # 2
                        anom.x * 0 +                  # 3
                        area.x * area.y +             # 4
                        open.x * 0 +                  # 5 
                        ratio.x * 0,                  # 6
                      area.ratio * 0 +              # 7
                        anom.yr * 0 +                 # 8
                        area.yr * 0  +                # 9
                        open.yr * 0 +                 # 10
                        ratio.yr * 0 +                # 11
                        area.ratio.yr * 0             # 12
    )) |> 
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
  dplyr::mutate(
    psi = plogis(       intercept +                   # 1
                          yr * 0 +                 # 2
                          anom.x * 0 +                    # 3
                          area.x * 0 +             # 4
                          open.x * 0 +                    # 5 
                          ratio.x * ratio.y,                    # 6
                        area.ratio * 0 +              # 7
                          anom.yr * 0 +                 # 8
                          area.yr * 0  +    # 9
                          open.yr * 0 +                 # 10
                          ratio.yr * 0 +                # 11
                          area.ratio.yr * 0             # 12
    )) |> 
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

fit_anom <- dplyr::cross_join(
  post, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(         intercept +                   # 1
                            yr * 0 +                 # 2
                            anom.x * anom.y +                    # 3
                            area.x * 0 +             # 4
                            open.x * 0 +                    # 5 
                            ratio.x * 0,                    # 6
                          area.ratio * 0 +              # 7
                            anom.yr * 0 +                 # 8
                            area.yr * 0  +    # 9
                            open.yr * 0 +                 # 10
                            ratio.yr * 0 +                # 11
                            area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(iter, anom = anom.y, psi ) |> 
  dplyr::group_by(anom) |> 
  dplyr::summarise( mean = mean(psi), 
                    l95 = quantile(psi, c(0.025)), 
                    l68 = quantile(psi, c(0.160)), 
                    u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( anom = anom*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center")) |> 
  dplyr::rename(cov = anom) |> 
  tibble::add_column( cov_name = "Precip. anomaly (mm)")

com_fit <- dplyr::full_join(
  fit_area, 
  fit_open) |> 
  dplyr::full_join(fit_ratio) |> 
  dplyr::full_join(fit_anom) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
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
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr * 0 +                      # 2
                      anom.x * 0 +                    # 3
                      area.x * 0 +                    # 4
                      open.x * open.y +             # 5 
                      ratio.x * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * 0 * 0 +             # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0
    )) |> 
  dplyr::select(sp, iter, open = open.y, psi ) |> 
  dplyr::group_by( sp, open) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( open = open*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center")) |> 
  dplyr::rename(cov = open) |> 
  tibble::add_column( cov_name = "Percent open habitat")

fit_area_sp <- dplyr::cross_join(
  post_sp, 
  pdat
) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr * 0 +                      # 2
                      anom.x * 0 +                    # 3
                      area.x * area.y +                    # 4
                      open.x * 0 +             # 5 
                      ratio.x * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * 0 * 0 +             # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0
    )) |> 
  dplyr::select(sp, iter, area = area.y, psi ) |> 
  dplyr::group_by( sp, area) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( area = area*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center")) |> 
  dplyr::rename(cov = area) |> 
  dplyr::mutate(cov = exp(cov)) |> 
  tibble::add_column( cov_name = "Reserve area (km2)")

fit_ratio_sp <- dplyr::cross_join(
  post_sp, 
  pdat
) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr * 0 +                      # 2
                      anom.x * 0 +                    # 3
                      area.x *0 +                    # 4
                      open.x * 0 +             # 5 
                      ratio.x * ratio.y,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * 0 * 0 +             # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0
    )) |> 
  dplyr::select(sp, iter, ratio = ratio.y, psi ) |> 
  dplyr::group_by( sp, ratio) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( ratio = ratio*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center")) |> 
  dplyr::rename(cov = ratio) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio")

fit_anom_sp <- dplyr::cross_join(
  post_sp, 
  pdat
) |> 
  dplyr::mutate(
    psi = plogis(     intercept +                   # 1
                        yr * 0 +                      # 2
                        anom.x * anom.y +                    # 3
                        area.x *0 +                    # 4
                        open.x * 0 +             # 5 
                        ratio.x * 0,                    # 6
                      area.ratio * 0 +              # 7
                        anom.yr * 0 +                 # 8
                        area.yr * 0 +                 # 9
                        open.yr * 0 * 0 +             # 10
                        ratio.yr * 0 +                # 11
                        area.ratio.yr * 0
    )) |> 
  dplyr::select(sp, iter, anom = anom.y, psi ) |> 
  dplyr::group_by( sp, anom) |> 
  dplyr::summarise( mean = mean(psi)) |> 
  dplyr::mutate( anom = anom*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center")) |> 
  dplyr::rename(cov = anom) |> 
  tibble::add_column( cov_name = "Precip. anomaly (mm)")

sp_fit <- dplyr::full_join(fit_area_sp, fit_open_sp) |> 
  dplyr::full_join(fit_ratio_sp) |> 
  dplyr::full_join(fit_anom_sp) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Precip. anomaly (mm)")))

x_pos <- com_fit |>
  dplyr::group_by(cov_name) |>
  dplyr::summarise(min = min(cov), range = max(cov) - min(cov)) |>
  dplyr::mutate(min_plus = min + 0.05 * range) |>
  dplyr::select(cov_name, cov = min_plus)

com_stars <- stars |> 
  dplyr::filter( var %in% c("area", "open", "ratio", "anom")) |> 
  tibble::add_column(mean = 1.05) |> 
  tibble::add_column(cov_name = c("Precip. anomaly (mm)", "Reserve area (km2)", "Percent open habitat", "Perimeter-to-area ratio")) |> 
  dplyr::left_join(x_pos) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Precip. anomaly (mm)")))

( panel_a <- ggplot() + 
    scale_y_continuous(limits = c(0, 1.1), 
                       breaks = c(0, 0.5, 1)) +
    geom_line(data = sp_fit, aes(x = cov,
                                 y = mean,
                                 group = factor(sp)),
              alpha = 0.3,
              linewidth = 0.5,
              color = MetBrewer::MetPalettes$Hiroshige[[1]][10]) +
    facet_wrap(~cov_name, scales = "free_x")  +
    geom_ribbon( data = com_fit, aes(x = cov, ymin = l95, ymax = u95),
                 color = NA,
                 fill = "#a62517",
                 alpha = 0.3) +
    geom_ribbon( data = com_fit, aes(x = cov, ymin = l68, ymax = u68),
                 color = NA,
                 fill = "#a62517",
                 alpha = 0.4) +
    geom_line(data = com_fit, aes(x = cov, y = mean),
              color = "#a62517",
              linewidth = 2) + 
    geom_text( data = com_stars, aes(x = cov, y = mean, label = lab)) +
    theme_minimal() + 
    scale_color_viridis_d() +
    labs(x = "Covariate value", 
         y = "Occurrence probability",
         title = "(a)") +
    ggplot2::theme(axis.line = element_line( linewidth = 0.15, color = "black"),
                   panel.grid = element_line(linewidth = 0.1, color = "gray80"),
                   axis.title = element_text(size = 8, color = "black"),
                   strip.text = element_text(size = 8, color = "black"), 
                   axis.text = element_text(size = 7, color = "black"),
                   panel.background = element_rect(color = NA, fill = "white"), 
                   plot.background = element_rect(color = NA, fill = "white"),
                   plot.title = element_text(color = "black", size = 10),
                   legend.position = "none"))


fit_open <- dplyr::cross_join(
  post, 
  p_open) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr.x * yr.y +                 # 2
                      anom * 0 +                    # 3
                      area * 0 +                    # 4
                      open.x * open.y +              # 5 
                      ratio * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * yr.y * open.y +    # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0
    )) |> 
  dplyr::select(iter, yr = yr.y, open = open.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(open) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Percent open habitat")|> 
  dplyr::rename(cov = open) |> 
  dplyr::mutate(cov = cov*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center"))

post_sp <- MCMCvis::MCMCpstr( out, params = "beta", type = "chains")[[1]] |>
  reshape2::melt(c("sp", "param", "iter")) |> 
  tibble::as_tibble() |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

fit_open_sp <- dplyr::cross_join(
  post_sp, 
  p_open) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr.x * yr.y +                 # 2
                      anom * 0 +                    # 3
                      area * 0 +                    # 4
                      open.x * open.y +              # 5 
                      ratio * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * yr.y * open.y +    # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0)) |> 
  dplyr::select(sp, iter, yr = yr.y, open = open.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(sp, open) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Percent open habitat")|> 
  dplyr::rename(cov = open) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(cov = cov*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center"))

fit_area <- dplyr::cross_join(
  post, 
  p_area) |> 
  dplyr::mutate(
    psi = plogis(     intercept +                   # 1
                        yr.x * yr.y +                 # 2
                        anom * 0 +                    # 3
                        area.x * area.y +             # 4
                        open * 0 +                    # 5 
                        ratio * 0,                    # 6
                      area.ratio * 0 +              # 7
                        anom.yr * 0 +                 # 8
                        area.yr * yr.y * area.y  +    # 9
                        open.yr * 0 +                 # 10
                        ratio.yr * 0 +                # 11
                        area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(iter, yr = yr.y, area = area.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(area) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Reserve area (km2)")|> 
  dplyr::rename(cov = area) |> 
  dplyr::mutate(cov = cov*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center")) |> 
  dplyr::mutate(cov = exp(cov ))

fit_area_sp <- dplyr::cross_join(
  post_sp, 
  p_area) |> 
  dplyr::mutate(
    psi = plogis(     intercept +                   # 1
                        yr.x * yr.y +                 # 2
                        anom * 0 +                    # 3
                        area.x * area.y +             # 4
                        open * 0 +                    # 5 
                        ratio * 0,                    # 6
                      area.ratio * 0 +              # 7
                        anom.yr * 0 +                 # 8
                        area.yr * yr.y * area.y  +    # 9
                        open.yr * 0 +                 # 10
                        ratio.yr * 0 +                # 11
                        area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(sp, iter, yr = yr.y, area = area.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(sp, area) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Reserve area (km2)")|> 
  dplyr::rename(cov = area) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(cov = cov*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center")) |> 
  dplyr::mutate(cov = exp(cov ))

fit_ratio <- dplyr::cross_join(
  post, 
  p_ratio) |> 
  dplyr::mutate(
    psi = plogis(       intercept +                   # 1
                          yr.x * yr.y +                 # 2
                          anom * 0 +                    # 3
                          area * 0 +             # 4
                          open * 0 +                    # 5 
                          ratio.x * ratio.y,                    # 6
                        area.ratio * 0 +              # 7
                          anom.yr * 0 +                 # 8
                          area.yr * 0  +    # 9
                          open.yr * 0 +                 # 10
                          ratio.yr * ratio.y*yr.y +                # 11
                          area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(iter, yr = yr.y, ratio = ratio.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(ratio) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio")|> 
  dplyr::rename(cov = ratio) |> 
  dplyr::mutate(cov = cov*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center"))

fit_ratio_sp <- dplyr::cross_join(
  post_sp, 
  p_ratio) |> 
  dplyr::mutate(
    psi = plogis(       intercept +                   # 1
                          yr.x * yr.y +                 # 2
                          anom * 0 +                    # 3
                          area * 0 +             # 4
                          open * 0 +                    # 5 
                          ratio.x * ratio.y,                    # 6
                        area.ratio * 0 +              # 7
                          anom.yr * 0 +                 # 8
                          area.yr * 0  +    # 9
                          open.yr * 0 +                 # 10
                          ratio.yr * ratio.y*yr.y +                # 11
                          area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(sp, iter, yr = yr.y, ratio = ratio.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(sp, ratio) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio")|> 
  dplyr::rename(cov = ratio) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(cov = cov*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center"))

fit_anom <- dplyr::cross_join(
  post, 
  p_anom) |> 
  dplyr::mutate(
    psi = plogis(         intercept +                   # 1
                            yr.x * yr.y +                 # 2
                            anom.x * anom.y +                    # 3
                            area * 0 +             # 4
                            open * 0 +                    # 5 
                            ratio * 0,                    # 6
                          area.ratio * 0 +              # 7
                            anom.yr * yr.y *anom.y +                 # 8
                            area.yr * 0  +    # 9
                            open.yr * 0 +                 # 10
                            ratio.yr * 0 +                # 11
                            area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(iter, yr = yr.y, anom = anom.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(anom) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |>  
  tibble::add_column( cov_name = "Precip. anomaly (mm)") |> 
  dplyr::rename(cov = anom) |> 
  dplyr::mutate(cov = cov*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center"))

fit_anom_sp <- dplyr::cross_join(
  post_sp, 
  p_anom) |> 
  dplyr::mutate(
    psi = plogis(         intercept +                   # 1
                            yr.x * yr.y +                 # 2
                            anom.x * anom.y +                    # 3
                            area * 0 +             # 4
                            open * 0 +                    # 5 
                            ratio * 0,                    # 6
                          area.ratio * 0 +              # 7
                            anom.yr * yr.y *anom.y +                 # 8
                            area.yr * 0  +    # 9
                            open.yr * 0 +                 # 10
                            ratio.yr * 0 +                # 11
                            area.ratio.yr * 0             # 12
    )) |> 
  dplyr::select(sp, iter, yr = yr.y, anom = anom.y, psi ) |> 
  dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
  tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(sp, anom) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |>  
  tibble::add_column( cov_name = "Precip. anomaly (mm)") |> 
  dplyr::rename(cov = anom) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(cov = cov*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center"))

mu_betas <- dplyr::full_join(
  fit_area, 
  fit_open) |> 
  dplyr::full_join(fit_ratio) |> 
  dplyr::full_join(fit_anom) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Precip. anomaly (mm)")))

sp_fits <- dplyr::full_join(
  fit_area_sp, 
  fit_open_sp) |> 
  dplyr::full_join(fit_ratio_sp) |> 
  dplyr::full_join(fit_anom_sp) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Precip. anomaly (mm)")))

focal_stars <- stars |> 
  dplyr::filter( var %in% c("area.yr", "open.yr", "ratio.yr", "anom.yr")) |> 
  tibble::add_column( cov = c(-100, 0.2, 0.1, 0.005), 
                      test = 1.05) |> 
  tibble::add_column(cov_name = c("Precip. anomaly (mm)", "Reserve area (km2)", "Percent open habitat", "Perimeter-to-area ratio")) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Reserve area (km2)", 
                                     "Percent open habitat", 
                                     "Perimeter-to-area ratio", 
                                     "Precip. anomaly (mm)")))

( panel_b <- ggplot() +
    geom_hline(yintercept = 0.5, color = "black", linetype = "dashed",
               linewidth = 0.2) +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = c(0, 0.5, 1)) +
    geom_line(data = sp_fits, 
              aes(x = cov, 
                  y = test, 
                  group = factor(sp)), 
              color = MetBrewer::MetPalettes$Hiroshige[[1]][10], 
              alpha = 0.3,
              linewidth = 0.5) +
    geom_line( data = mu_betas, aes(x = cov, y = test), 
               linewidth = 2,
               color = "#a62517") +
    facet_wrap(~cov_name, scales = "free_x") +
    geom_text( data = focal_stars, aes(x = cov, y = test, label = lab)) +
    theme_minimal() + 
    labs(x = "Covariate value", 
         title = "(b)",
         y = "Probability of decline in occurrence from 2008 to 2023") +
    ggplot2::theme(axis.line = element_line( linewidth = 0.15, color = "black"),
                   panel.grid = element_line( linewidth = 0.1, color = "gray80"),
                   axis.title = element_text(size = 8, color = "black"),
                   strip.text = element_text(size = 8, color = "black"), 
                   axis.text = element_text(size = 7, color = "black"),
                   panel.background = element_rect(color = NA, fill = "white"), 
                   plot.background = element_rect(color = NA, fill = "white"),
                   plot.title = element_text(size = 10, color = "black")))

panel_a | panel_b 

setwd(here::here("figures"))
ggsave(
  "figure_03.png", 
  width = 6, 
  height = 3.5, 
  units = "in", 
  dpi = 600)