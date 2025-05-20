# create figure 2: estimates of species trends

library(here)
library(tidyverse)
library(MCMCvis)
library(MetBrewer)
library(ggh4x)

load(here::here("results/prairie_bird_trends_global2025-02-28.RData"))

sp_key <- readr::read_csv(here::here("data/sp_key.csv"))

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

sp_trends <- beta_df |>
  dplyr::left_join(param_key) |> 
  dplyr::left_join(sp_key) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  dplyr::select(sp, common, type, param, name, mean:neff) |> 
  dplyr::filter( name == "yr") |> 
  dplyr::arrange(type, mean)  |> 
  dplyr::mutate( common = factor(common, levels = unique(common)))

com <- MCMCvis::MCMCsummary( out, params = "mu_beta", probs = c(0.025, 0.160, 0.840, 0.975)) |> 
  tibble::as_tibble(rownames = "param") |> 
  dplyr::mutate( param = readr::parse_number(param)) |> 
  dplyr::left_join(param_key) |> 
  dplyr::filter(name == "yr") |> 
  dplyr::select( param, name, mean, l68 = `16%`, u68 = `84%`, l95 = `2.5%`, u95 = `97.5%`, rhat = Rhat, neff = n.eff ) |> 
  dplyr::cross_join(
    tibble::tibble( type = c("Obligate", "Facultative"))) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative")))

ggplot2::ggplot() +
  ggplot2::facet_wrap(~type, ncol = 1, scales = "free_y") +
  ggh4x::force_panelsizes(rows = c(1,3)) +
  ggplot2::geom_rect(data = com, 
                     aes( xmin = l95, xmax = u95), 
                     ymin = -Inf, 
                     ymax = Inf,
                     fill = "#a62517",
                     alpha = 0.25) +
  ggplot2::geom_rect(data = com, 
                     aes( xmin = l68, xmax = u68), 
                     ymin = -Inf, 
                     ymax = Inf,
                     fill = "#a62517",
                     alpha = 0.25) +
  ggplot2::geom_vline(xintercept = 0,
                      color = "black",
                      linetype = "dashed",
                      linewidth = 0.25) +
  ggplot2::geom_errorbar(data = sp_trends, aes(y = common, xmin = l68, xmax = u68),
                         width = 0,
                         color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                         linewidth = 1) +
  ggplot2::geom_errorbar(data = sp_trends, aes(y = common, xmin = l95, xmax = u95),
                         color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                         width = 0) +
  ggplot2::geom_point(data = sp_trends, aes(y = common, x = mean),
                      color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                      size = 2.25) +
  ggplot2::theme_minimal() +
  ggplot2::xlab("Occurrence trend") +
  ggplot2::theme(axis.title.y = element_blank(),
                 axis.line = element_line( linewidth = 0.15, color = "black"),
                 panel.grid = element_line( linewidth = 0.1, color = "gray80"),
                 axis.title.x = element_text(size = 10, color = "black"),
                 axis.text.x = element_text(size = 9, color = "black"),
                 axis.text.y = ggtext::element_markdown(color = "black", size = 9),
                 panel.background = element_rect(color = NA, fill = "white"), 
                 plot.background = element_rect(color = NA, fill = "white"),
                 strip.text = element_text(size = 10, color = "black", face = "bold"))  

setwd(here::here("figures"))
ggsave(
  filename = "figure_02.png", 
  width = 3, 
  height = 4.5, 
  units = "in", 
  dpi = 600)

sp_trends |> 
  dplyr::mutate( 
    no_trend = ifelse( l68 < 0 & u68 > 0, 1, 0),
    decline68 = ifelse(u68 < 0, 1, 0),
    decline95 = ifelse(u95 < 0, 1, 0),
    increase68 = ifelse(l68 > 0, 1, 0), 
    increase95 = ifelse(l95 > 0, 1, 0)) |> 
  dplyr::group_by(type) |> 
  dplyr::summarise( 
    no_trend = sum(no_trend),
    decline68 = sum(decline68), 
    decline95 = sum(decline95), 
    increase68 = sum(increase68), 
    increase95 = sum(increase95)) 

sp_trends |> 
  dplyr::mutate( 
    no_trend = ifelse( l68 < 0 & u68 > 0, 1, 0),
    decline68 = ifelse(u68 < 0, 1, 0),
    decline95 = ifelse(u95 < 0, 1, 0),
    increase68 = ifelse(l68 > 0, 1, 0), 
    increase95 = ifelse(l95 > 0, 1, 0)) |> 
  dplyr::summarise( 
    no_trend = sum(no_trend),
    decline68 = sum(decline68), 
    decline95 = sum(decline95), 
    increase68 = sum(increase68), 
    increase95 = sum(increase95)) |> mutate(across(no_trend:increase95, function(x) round(x/33, 2)))