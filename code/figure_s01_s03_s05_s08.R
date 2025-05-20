# make figures S1, S3, S5, S8
library(here)
library(tidyverse)
library(MCMCvis)
library(MetBrewer)

load(here::here("results/prairie_bird_trends_global2025-02-28.RData"))

key <- readr::read_csv(here::here("data/sp_key.csv"))

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

pdat <- tibble::tibble(
  yr = mean(year.sc)) |> 
  dplyr::cross_join(
    tibble::tibble(
      open = seq(from = min(open.sc), to = max(open.sc), length.out = 10), 
      area = seq(from = min(area.sc), to = max(area.sc), length.out = 10),
      ratio = seq(from = min(ratio.sc), to = max(ratio.sc), length.out = 10), 
      anom = seq(from = min(anom.sc), to = max(anom.sc), length.out = 10))) 

post_sp <- MCMCvis::MCMCpstr( out, params = "beta", type = "chains")[[1]] |>
  reshape2::melt(c("sp", "param", "iter")) |> 
  tibble::as_tibble() |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

stars <- post_sp |> 
  group_by(sp) |> 
  dplyr::select(sp, iter, yr:area.ratio.yr) |> 
  tidyr::pivot_longer(yr:area.ratio.yr, names_to = "var", values_to = "val") |> 
  dplyr::group_by(sp, var) |> 
  dplyr::summarise( l95 = quantile(val, c(0.025)), 
                    u95 = quantile(val, c(0.975)), 
                    l68 = quantile(val, c(0.160)), 
                    u68 = quantile(val, c(0.840))) |> 
  mutate(sig68 = ifelse(l68 < 0 & u68 > 0, 0, 1), 
         sig95 = ifelse(l95 < 0 & u95 > 0, 0, 1)) |> 
  # dplyr::left_join(int_key) |> 
  mutate(lab = ifelse( sig95 == 1, "**", 
                       ifelse(sig95 == 0 & sig68 == 1, "*", ""))) |> 
  dplyr::select(var, lab) |> 
  left_join(key)

stars_area_sp <- stars |> 
  dplyr::ungroup() |> 
  dplyr::filter(var == "area") |> 
  tibble::add_column( mean = 1.1,
                      cov = -0.25) |> 
  dplyr::mutate(code = toupper(code))


fit_area_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(           intercept +                   # 1
                              yr.x * yr.y +                      # 2
                              anom.x * 0 +                    # 3
                              area.x * area.y +                    # 4
                              open.x * 0 +             # 5 
                              ratio.x * 0,                    # 6
                            area.ratio * 0 +              # 7
                              anom.yr * 0 +                 # 8
                              area.yr * area.y * yr.y +                 # 9
                              open.yr * 0 * 0 +             # 10
                              ratio.yr * 0 +                # 11
                              area.ratio.yr * 0)) |> 
  dplyr::select(sp, iter, yr = yr.y, area = area.y, psi ) |> 
  dplyr::group_by( sp, area, yr) |> 
  dplyr::summarise( mean = mean(psi),
                    l95 = quantile(psi, c(0.025)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( area = area*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center"),
                 yr = yr*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::rename(cov = area) |> 
  tibble::add_column( cov_name = "Reserve area") |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code)) |> 
  dplyr::mutate( cov = exp(cov))

ggplot() +
  geom_ribbon( data = fit_area_sp, aes(x = cov, ymin = l95, ymax = u95), color = NA, alpha = 0.3) +
  geom_line( data = fit_area_sp, aes(x = cov, y = mean), linewidth = 1.5) +
  geom_text(data = stars_area_sp, aes(x = cov, y = mean, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  facet_wrap(~code) +
  labs(x = "Reserve area (km2)", 
       y = "Occurrence probability") +
  theme_minimal() +
  theme( panel.grid = element_blank(), 
         panel.background = element_rect(color = NA, fill = "white"), 
         plot.background = element_rect(color = NA, fill = "white"), 
         strip.text = element_text(color = "black", size = 9), 
         axis.title = element_text(color = "black", size = 9), 
         axis.text = element_text(color = "black", size = 8),
         axis.line = element_line(color = "black", linewidth = 0.1))

setwd(here::here('figures'))
ggsave(
  filename = "figure_s01.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

stars_open_sp <- stars |> 
  dplyr::ungroup() |> 
  dplyr::filter(var == "open") |> 
  tibble::add_column( mean = 1.1,
                      cov = 0.05) |> 
  dplyr::mutate(code = toupper(code))

fit_open_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(   intercept +                   # 1
                      yr.x * yr.y +                      # 2
                      anom.x * 0 +                    # 3
                      area.x * 0 +                    # 4
                      open.x * open.y +             # 5 
                      ratio.x * 0,                    # 6
                    area.ratio * 0 +              # 7
                      anom.yr * 0 +                 # 8
                      area.yr * 0 +                 # 9
                      open.yr * yr.y * open.y +             # 10
                      ratio.yr * 0 +                # 11
                      area.ratio.yr * 0)) |> 
  dplyr::select(sp, iter, yr = yr.y, open = open.y, psi ) |> 
  dplyr::group_by( sp, open, yr) |> 
  dplyr::summarise( mean = mean(psi),
                    l95 = quantile(psi, c(0.025)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( open = open*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center"),
                 yr = yr*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::rename(cov = open) |> 
  tibble::add_column( cov_name = "Percent open habitat") |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

ggplot() +
  geom_ribbon( data = fit_open_sp, aes(x = cov, ymin = l95, ymax = u95), color = NA, alpha = 0.3) +
  geom_line( data = fit_open_sp, aes(x = cov, y = mean), linewidth = 1.5) +
  geom_text(data = stars_open_sp, aes(x = cov, y = mean, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  facet_wrap(~code) +
  labs(x = "Proportion open habitat (250 m)", 
       y = "Occurrence probability") +
  theme_minimal() +
  theme( panel.grid = element_blank(), 
         panel.background = element_rect(color = NA, fill = "white"), 
         plot.background = element_rect(color = NA, fill = "white"), 
         strip.text = element_text(color = "black", size = 9), 
         axis.title = element_text(color = "black", size = 9), 
         axis.text = element_text(color = "black", size = 8),
         axis.line = element_line(color = "black", linewidth = 0.1))

setwd(here::here('figures'))
ggsave(
  filename = "figure_s03.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

stars_ratio_sp <- stars |> 
  dplyr::ungroup() |> 
  dplyr::filter(var == "ratio") |> 
  tibble::add_column( mean = 1.1,
                      cov = 0.005) |> 
  dplyr::mutate(code = toupper(code))


fit_ratio_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(       intercept +                   # 1
                          yr.x * yr.y +                      # 2
                          anom.x * 0 +                    # 3
                          area.x *0 +                    # 4
                          open.x * 0 +             # 5 
                          ratio.x * ratio.y,                    # 6
                        area.ratio * 0 +              # 7
                          anom.yr * 0 +                 # 8
                          area.yr * 0 +                 # 9
                          open.yr * 0 * 0 +             # 10
                          ratio.yr * yr.y*ratio.y +                # 11
                          area.ratio.yr * 0
    )) |> 
  dplyr::select(sp, iter, yr = yr.y, ratio = ratio.y, psi ) |> 
  dplyr::group_by( sp, ratio, yr) |> 
  dplyr::summarise( mean = mean(psi),
                    l95 = quantile(psi, c(0.025)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( ratio = ratio*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center"),
                 yr = yr*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::rename(cov = ratio) |> 
  tibble::add_column( cov_name = "Perimeter-to-area ratio") |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

ggplot() +
  geom_ribbon( data = fit_ratio_sp, aes(x = cov, ymin = l95, ymax = u95), color = NA, alpha = 0.3) +
  geom_line( data = fit_ratio_sp, aes(x = cov, y = mean), linewidth = 1.5) +
  geom_text(data = stars_ratio_sp, aes(x = cov, y = mean, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  facet_wrap(~code) +
  labs(x = "Perimeter-to-area ratio", 
       y = "Occurrence probability") +
  theme_minimal() +
  theme( panel.grid = element_blank(), 
         panel.background = element_rect(color = NA, fill = "white"), 
         plot.background = element_rect(color = NA, fill = "white"), 
         strip.text = element_text(color = "black", size = 9), 
         axis.title = element_text(color = "black", size = 9), 
         axis.text = element_text(color = "black", size = 8),
         axis.line = element_line(color = "black", linewidth = 0.1))

setwd(here::here('figures'))
ggsave(
  filename = "figure_s05.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

stars_anom_sp <- stars |> 
  dplyr::ungroup() |> 
  dplyr::filter(var == "anom") |> 
  tibble::add_column( mean = 1.1,
                      cov = -125) |> 
  dplyr::mutate(code = toupper(code))


fit_anom_sp <- dplyr::cross_join(
  post_sp, 
  pdat) |> 
  dplyr::mutate(
    psi = plogis(         intercept +                   # 1
                            yr.x * yr.y +                      # 2
                            anom.x * anom.y +                    # 3
                            area.x *0 +                    # 4
                            open.x * 0 +             # 5 
                            ratio.x * 0,                    # 6
                          area.ratio * 0 +              # 7
                            anom.yr * anom.y*yr.y +                 # 8
                            area.yr * 0 +                 # 9
                            open.yr * 0 * 0 +             # 10
                            ratio.yr * 0 +                # 11
                            area.ratio.yr * 0)) |> 
  dplyr::select(sp, iter, yr = yr.y, anom = anom.y, psi ) |> 
  dplyr::group_by( sp, anom, yr) |> 
  dplyr::summarise( mean = mean(psi),
                    l95 = quantile(psi, c(0.025)), 
                    u95 = quantile(psi, c(0.975))) |> 
  dplyr::mutate( anom = anom*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center"),
                 yr = yr*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::rename(cov = anom) |> 
  tibble::add_column( cov_name = "Precipitation anomaly") |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

ggplot() +
  geom_ribbon( data = fit_anom_sp, aes(x = cov, ymin = l95, ymax = u95), color = NA, alpha = 0.3) +
  geom_line( data = fit_anom_sp, aes(x = cov, y = mean), linewidth = 1.5) +
  geom_text(data = stars_anom_sp, aes(x = cov, y = mean, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  facet_wrap(~code) +
  labs(x = "Precipitation anomaly", 
       y = "Occurrence probability") +
  theme_minimal() +
  theme( panel.grid = element_blank(), 
         panel.background = element_rect(color = NA, fill = "white"), 
         plot.background = element_rect(color = NA, fill = "white"), 
         strip.text = element_text(color = "black", size = 9), 
         axis.title = element_text(color = "black", size = 9), 
         axis.text = element_text(color = "black", size = 8),
         axis.line = element_line(color = "black", linewidth = 0.1))

setwd(here::here('figures'))
ggsave(
  filename = "figure_s08.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)