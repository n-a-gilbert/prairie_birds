# make figure s07 to visualize interactions between trend, reserve size, reserve shape
library(here)
library(tidyverse)
library(MCMCvis)
library(MetBrewer)
library(ggh4x)

load(here::here("results/prairie_bird_trends_global2025-02-28.RData"))

sp_key <- readr::read_csv(here::here("data/sp_key.csv"))

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

param_key <- tibble::tibble(
  param = 1:12,
  name = c("intercept", "yr", "anom", "area", "open", "ratio",
           "area.ratio",
           "anom.yr", "area.yr", "open.yr", "ratio.yr",
           "area.ratio.yr"))

betas <- MCMCvis::MCMCsummary( out, params = "beta", probs = c(0.025, 0.160, 0.840, 0.975))

beta_df <- tibble::as_tibble(betas, rownames = "param") |> 
  tidyr::separate(param, into = c("species", "param"), sep = ",") |> 
  dplyr::mutate(across(species:param, function(x) readr::parse_number(x))) |> 
  dplyr::select( sp = species, param, mean, l68 = `16%`, u68 = `84%`, l95 = `2.5%`, u95 = `97.5%`, rhat = Rhat, neff = n.eff )

obligates <- c( "Bobolink", "Dickcissel", "Grasshopper Sparrow", "LeConte's Sparrow", "Savannah Sparrow", 
                "Sedge Wren", "Upland Sandpiper", "Western Meadowlark")

pdat <- tidyr::expand_grid(
  yr = seq(from = min(year.sc), to = max(year.sc), length.out = 10),
  # area = seq(from = min(area.sc), to = max(area.sc), length.out = 20),
  area = c(min(area.sc), max(area.sc)),
  ratio = c(min(ratio.sc), max(ratio.sc)))

post_sp <- MCMCvis::MCMCpstr( out, params = "beta", type = "chains")[[1]] |>
  reshape2::melt(c("sp", "param", "iter")) |> 
  tibble::as_tibble() |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

tmp <-
  dplyr::cross_join(
    post_sp, 
    pdat) |> 
  dplyr::mutate(
    psi = plogis(       intercept +                   # 1
                          yr.x * yr.y +                      # 2
                          anom * 0 +                  # 3
                          area.x * area.y +             # 4
                          open * 0 +                  # 5 
                          ratio.x * ratio.y,                  # 6
                        area.ratio * area.y * ratio.y +              # 7
                          anom.yr * 0 +                 # 8
                          area.yr * area.y * yr.y  +                # 9
                          open.yr * 0 +                 # 10
                          ratio.yr * ratio.y * yr.y +                # 11
                          area.ratio.yr * area.y * ratio.y * yr.y             # 12
    )) |> 
  dplyr::select(sp, iter, yr = yr.y, area = area.y, ratio = ratio.y, psi ) |>
  # dplyr::mutate( yr = ifelse(yr < 0, "2008", "2023")) |> 
  # tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
  # dplyr::mutate(diff = `2023` - `2008`) |> 
  dplyr::group_by(sp, area, ratio, yr) |> 
  # dplyr::summarise(p = sum(diff < 0) / n())
  dplyr::summarise( mean = mean(psi),
                    l95 = quantile(psi, c(0.025)),
                    # l68 = quantile(psi, c(0.160)), 
                    # u68 = quantile(psi, c(0.840)), 
                    u95 = quantile(psi, c(0.975)))

p <- tmp |> 
  mutate(yr = round( yr*attr(year.sc, "scaled:scale") + attr(year.sc, "scaled:center"), 0)) |> 
  left_join(sp_key) |> 
  dplyr::mutate( type = ifelse(area < 0 & ratio < 0, "Small, simple", 
                               ifelse(area < 0 & ratio > 0, "Small, complex", 
                                      ifelse(area > 0 & ratio < 0, "Large, simple", 
                                             "Large, complex")))) |> 
  ggplot(aes(x = yr, y = mean, color = type)) +
  # geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = l95, ymax = u95, fill = type), alpha = 0.3, color = NA) +
  geom_line( linewidth = 1) + 
  facet_wrap(~common) +
  scale_color_manual(
    values = MetBrewer::MetPalettes$Hiroshige[[1]][c(1,3,6,8)]) +
  scale_fill_manual(
    values = MetBrewer::MetPalettes$Hiroshige[[1]][c(1,3,6,8)]) +
  theme_minimal() +
  scale_x_continuous(breaks = c(2008, 2015, 2023)) +
  labs(x = "Year", y = "Probability of occurrence", color = "Reserve", fill = "Reserve") +
  
  theme(axis.text.x = element_text(angle = 90, color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text = element_text(color = "black", size = 8), 
        axis.title = element_text(color = "black", size = 8), 
        # legend.title = element_text(color = "black", size = 8),
        legend.position = "bottom",
        axis.line = element_line(color = "black", linewidth = 0.2), 
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.box.margin = margin(-15, 0, 0, 0, unit = "pt"),
        legend.text = element_text(color = "black", size = 7),
        panel.background = element_rect(color = NA, fill = "white"), 
        plot.background = element_rect(color = NA, fill = "white"))  +
  guides(fill = guide_legend(ncol = 4), 
         color = guide_legend(ncol = 4))

setwd(here::here("figures"))
ggsave(
  "figure_s07.png", 
  width = 8, 
  height = 8,
  units = "in", 
  dpi = 600)
