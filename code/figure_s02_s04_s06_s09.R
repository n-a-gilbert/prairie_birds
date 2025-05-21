# make figures S2, S4, S6, S9
library(here)
library(tidyverse)
library(MCMCvis)
library(MetBrewer)

load(here::here("results/prairie_bird_trends_global2025-02-28.RData"))

key <- readr::read_csv(here::here("data/sp_key.csv"))

# read in bird data; contains columns for the focal covariates as well
d <- readr::read_csv(here::here("data/mn_prairie_bird_data_clean.csv")) |> 
  dplyr::filter(!is.na(site))

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

p.area <- tidyr::expand_grid( area = seq( from = min(area.sc), 
                                          to = max(area.sc),
                                          length.out = 10),
                              year = c( min(year.sc), 
                                        max(year.sc))) |> 
  # dplyr::cross_join( key ) |> 
  dplyr::mutate(area.unscaled = area*attr(area.sc, "scaled:scale") + attr(area.sc, "scaled:center"),
                year.unscaled = year*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::select(#sp, code, common,
    year, area, area.unscaled, year.unscaled) |> 
  dplyr::mutate(area.unscaled = exp( area.unscaled ))

area_p <- post_sp |> 
  dplyr::select(sp, iter, intercept, yr, area, area.yr) |> 
  dplyr::cross_join(
    p.area) |> 
  dplyr::mutate( p = plogis( intercept + yr*year + area.x*area.y + area.yr * year*area.y )) |> 
  dplyr::select(sp, iter, area = area.unscaled, yr = year.unscaled, p) |> 
  tidyr::pivot_wider(names_from = yr, values_from = p) |> 
  dplyr::mutate( diff = `2023` - `2008`) |>
  dplyr::group_by(sp, area) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  # dplyr::summarise( test = ecdf(diff)(0)) |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

stars_area <- stars |> 
  dplyr::filter(var == "area.yr") |> 
  dplyr::ungroup() |> 
  dplyr::mutate(code = toupper(code)) |> 
  tibble::add_column( area = 0.5, 
                      test = 1.1)

ggplot() + 
  geom_line(data = area_p, aes(x = area, y = test), 
            linewidth = 1.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_hline(yintercept = 0.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][1], linetype = "dashed") +
  theme_minimal() +
  facet_wrap(~code) +
  geom_text( data = stars_area, 
             aes(x = area, y = test, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  labs( y = "Probability of decline in occurrence from 2008 to 2023", 
        x = "Reserve area (km2)") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 9, color = "black"), 
        axis.text = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 9, color = "black"),
        panel.background = element_rect(fill = "white", 
                                        color = NA), 
        plot.background = element_rect(fill ="white", 
                                       color = NA))

setwd(here::here("figures"))
ggsave(
  filename = "figure_s02.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

p.open <- tidyr::expand_grid( open = seq( from = min(open.sc), 
                                          to = max(open.sc),
                                          length.out = 50),
                              year = c( min(year.sc), 
                                        max(year.sc))) |> 
  dplyr::mutate(open.unscaled = open*attr(open.sc, "scaled:scale") + attr(open.sc, "scaled:center"),
                year.unscaled = year*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::select(year, open, open.unscaled, year.unscaled)

open_p <-
  post_sp |> 
  dplyr::select(sp, iter, intercept, yr, open, open.yr) |> 
  dplyr::cross_join(
    p.open) |> 
  dplyr::mutate( p = plogis( intercept + yr*year + open.x*open.y + open.yr*year*open.y )) |> 
  dplyr::select(sp, iter, open = open.unscaled, yr = year.unscaled, p) |> 
  tidyr::pivot_wider(names_from = yr, values_from = p) |> 
  dplyr::mutate( diff = `2023` - `2008`) |>
  dplyr::group_by(sp, open) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  # dplyr::summarise( test = ecdf(diff)(0)) |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

stars_open <- stars |> 
  dplyr::filter(var == "open.yr") |> 
  dplyr::ungroup() |> 
  dplyr::mutate(code = toupper(code)) |> 
  tibble::add_column( open = 0.1, 
                      test = 1.1)

ggplot() +
  facet_wrap(~code) + 
  geom_line(data = open_p, aes(x = open, y = test),
            linewidth = 1.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_hline(yintercept = 0.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][1], linetype = "dashed") +
  theme_minimal() +
  geom_text(data = stars_open, aes(x = open, y = test, label = lab)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  labs( y = "Probability of decline in occurrence from 2008 to 2023", 
        x = "Proportion open habitat (250m)") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 9, color = "black"), 
        axis.text = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 9, color = "black"),
        panel.background = element_rect(fill = "white", 
                                        color = NA), 
        plot.background = element_rect(fill ="white", 
                                       color = NA))

setwd(here::here("figures"))
ggsave(
  filename = "figure_s04.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

p.ratio <- tidyr::expand_grid( ratio = seq( from = min(ratio.sc), 
                                            to = max(ratio.sc),
                                            length.out = 10),
                               year = c( min(year.sc), 
                                         max(year.sc))) |> 
  # dplyr::cross_join( key ) |> 
  dplyr::mutate(ratio.unscaled = ratio*attr(ratio.sc, "scaled:scale") + attr(ratio.sc, "scaled:center"),
                year.unscaled = year*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::select(#sp, code, common,
    year, ratio, ratio.unscaled, year.unscaled)

ratio_p <-
  post_sp |> 
  dplyr::select(sp, iter, intercept, yr, ratio, ratio.yr) |> 
  dplyr::cross_join(
    p.ratio) |> 
  dplyr::mutate( p = plogis( intercept + yr*year + ratio.x*ratio.y + ratio.yr*year*ratio.y )) |> 
  dplyr::select(sp, iter, ratio = ratio.unscaled, yr = year.unscaled, p) |> 
  tidyr::pivot_wider(names_from = yr, values_from = p) |> 
  dplyr::mutate( diff = `2023` - `2008`) |>
  dplyr::group_by(sp, ratio) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  # dplyr::summarise( test = ecdf(diff)(0)) |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

stars_ratio <- stars |> 
  dplyr::filter(var == "ratio.yr") |> 
  dplyr::ungroup() |> 
  dplyr::mutate(code = toupper(code)) |> 
  tibble::add_column( ratio = 0.005, 
                      test = 1.1)

ggplot() + 
  geom_line(data = ratio_p, aes(x = ratio, y = test),
            linewidth = 1.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_hline(yintercept = 0.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][1], linetype = "dashed") +
  geom_text(data = stars_ratio, aes(x = ratio, y = test, label = lab)) +
  theme_minimal() +
  facet_wrap(~code) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  labs( y = "Probability of decline in occurrence from 2008 to 2023", 
        x = "Perimeter-to-area ratio") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 9, color = "black"), 
        axis.text = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 9, color = "black"),
        panel.background = element_rect(fill = "white", 
                                        color = NA), 
        plot.background = element_rect(fill ="white", 
                                       color = NA))
setwd(here::here("figures"))
ggsave(
  filename = "figure_s06.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)

p.anom <- tidyr::expand_grid( anom = seq( from = min(anom.sc), 
                                          to = max(anom.sc),
                                          length.out = 10),
                              year = c( min(year.sc), 
                                        max(year.sc))) |> 
  # dplyr::cross_join( key ) |> 
  dplyr::mutate(anom.unscaled = anom*attr(anom.sc, "scaled:scale") + attr(anom.sc, "scaled:center"),
                year.unscaled = year*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center")) |> 
  dplyr::select(#sp, code, common,
    year, anom, anom.unscaled, year.unscaled)

anom_p <-
  post_sp |> 
  dplyr::select(sp, iter, intercept, yr, anom, anom.yr) |> 
  dplyr::cross_join(
    p.anom) |> 
  dplyr::mutate( p = plogis( intercept + yr*year + anom.x*anom.y + anom.yr*year*anom.y )) |> 
  dplyr::select(sp, iter, anom = anom.unscaled, yr = year.unscaled, p) |> 
  tidyr::pivot_wider(names_from = yr, values_from = p) |> 
  dplyr::mutate( diff = `2023` - `2008`) |>
  dplyr::group_by(sp, anom) |> 
  dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
  # dplyr::summarise( test = ecdf(diff)(0)) |> 
  dplyr::left_join(key) |> 
  dplyr::mutate(code = toupper(code))

stars_anom <- stars |> 
  dplyr::filter(var == "anom.yr") |> 
  dplyr::ungroup() |> 
  dplyr::mutate(code = toupper(code)) |> 
  tibble::add_column( anom = -130, 
                      test = 1.1)

ggplot() +
  geom_line(data = anom_p, aes(x = anom, y = test),
            linewidth = 1.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_hline(yintercept = 0.5, color = MetBrewer::MetPalettes$Hiroshige[[1]][1], linetype = "dashed") +
  geom_text(data = stars_anom, aes(x = anom, y = test, label = lab)) +
  theme_minimal() +
  facet_wrap(~code) + 
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  # scale_x_continuous(breaks = c(0.01, 0.04)) +
  labs( y = "Probability of decline in occurrence from 2008 to 2023", 
        x = "Precipitation anomaly (mm)") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.title = element_text(size = 9, color = "black"), 
        axis.text = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 9, color = "black"),
        panel.background = element_rect(fill = "white", 
                                        color = NA), 
        plot.background = element_rect(fill ="white", 
                                       color = NA))
setwd(here::here("figures"))
ggsave(
  filename = "figure_s09.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 600)