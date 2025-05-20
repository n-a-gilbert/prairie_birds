# figure 4: visualize species-level parameters
# make a heatmap for main effects of covariates on occurrence and influence on trend

library(here)
library(tidyverse)
library(MCMCvis)
library(patchwork)
library(reshape2)
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

# posterior of species-level parameters
post_sp <- MCMCvis::MCMCpstr( out, params = "beta", type = "chains")[[1]] |>
  reshape2::melt(c("sp", "param", "iter")) |> 
  tibble::as_tibble() |> 
  dplyr::left_join(param_key) |> 
  dplyr::select(-param) |> 
  tidyr::pivot_wider(names_from = name, values_from = value)

# classify "significance" of parameters
stars <- post_sp |> 
  dplyr::group_by(sp) |> 
  dplyr::select(sp, iter, yr:area.ratio.yr) |> 
  tidyr::pivot_longer(yr:area.ratio.yr, names_to = "var", values_to = "val") |> 
  dplyr::group_by(sp, var) |> 
  dplyr::summarise( l95 = quantile(val, c(0.025)), 
                    u95 = quantile(val, c(0.975)), 
                    l68 = quantile(val, c(0.160)), 
                    u68 = quantile(val, c(0.840)),
                    mean = mean(val)) |> 
  dplyr::mutate(sig68 = ifelse(l68 < 0 & u68 > 0, 0, 1), 
                sig95 = ifelse(l95 < 0 & u95 > 0, 0, 1)) |> 
  dplyr::mutate(lab = ifelse( sig95 == 1, "**", 
                              ifelse(sig95 == 0 & sig68 == 1, "*", ""))) |> 
  dplyr::select(var, sig68, sig95, mean) |> 
  dplyr::left_join(sp_key)

param_name_key <- tibble::tibble(
  var = c("anom.yr", "area.yr", "open.yr", "ratio.yr"),
  cov_name = c("Precip.",
               "Area",
               "Open",
               "Ratio"))

obligates <- c( "Bobolink", "Dickcissel", "Grasshopper Sparrow", "LeConte's Sparrow", "Savannah Sparrow", 
                "Sedge Wren", "Upland Sandpiper", "Western Meadowlark")

trend_params <- stars |> 
  dplyr::filter(grepl(".yr", var)) |> 
  dplyr::filter(!var == "area.ratio.yr") |> 
  dplyr::mutate( lab = ifelse(sig68 == 0, "No effect", 
                              ifelse( sig95 == 1 & mean > 0, "Buffers decline (95% confidence)",
                                      ifelse(sig95 == 1 & mean < 0, "Accelerates decline (95% confidence)",
                                             ifelse(sig68 == 1 & sig95 == 0 & mean > 0, "Buffers decline (68% confidence)",
                                                    ifelse(sig68 == 1 & sig95 == 0 & mean < 0, "Accelerates decline (68% confidence)", NA)))))) |> 
  dplyr::mutate( lab = 
                   factor(lab, 
                          levels = c(
                            "Accelerates decline (95% confidence)",
                            "Accelerates decline (68% confidence)",
                            "No effect",
                            "Buffers decline (68% confidence)",
                            "Buffers decline (95% confidence)"))) |> 
  dplyr::left_join(param_name_key) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  dplyr::mutate( cov_name = factor(cov_name, 
                                   levels = c(
                                     "Area", "Open", "Ratio", "Precip."
                                   ))) |> 
  dplyr::mutate(common = factor(common, levels = unique(common))) |> 
  dplyr::mutate(common = factor(common, levels = rev(levels(common))))

( panelb <- ggplot( trend_params, aes(x = cov_name, 
                                      y = common, 
                                      fill = lab)) +
    facet_wrap(~type, ncol = 1, scales = "free_y",
               strip.position = "right") +
    ggh4x::force_panelsizes(rows = c(0.2424242,0.7575758)) +
    geom_tile(color = "black",
              linewidth = 0.1) +
    scale_fill_manual(values = c(
      MetBrewer::MetPalettes$Isfahan1[[1]][2],
      MetBrewer::MetPalettes$Isfahan1[[1]][4],
      "gray90",
      MetBrewer::MetPalettes$Isfahan1[[1]][5],
      MetBrewer::MetPalettes$Isfahan1[[1]][7])) +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    ggtitle("(b)   trends") +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = -0.05,
                                     vjust = -20,
                                     color = "black",
                                     size = 9),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(color = "black", size = 10, face = "bold"),
          plot.title = element_text(color = "black", size = 10, face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size = 8),
          legend.key.height = unit(1.5, "pt"),
          legend.justification = "left",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", color = NA), 
          panel.background = element_rect(fill = "white", color = NA),
          legend.margin = margin(-5, 0, 0, 0)) +
    guides(fill = guide_legend(ncol = 1)) )

param_name_key2 <- tibble(
  var = c("anom", "area", "open", "ratio"), 
  var_name = c("Precip.", "Area", "Open", "Ratio"))

occ_params <- stars |> 
  dplyr::filter(!grepl("\\.", var)) |> 
  dplyr::filter(!var=="yr") |> 
  dplyr::mutate( lab = ifelse(sig68 == 0,
                              "No effect",
                              ifelse( sig95 == 1 & mean > 0, "Increases occurrence (95% confidence)",
                                      ifelse(sig95 == 1 & mean < 0, "Decreases occurrence (95% confidence)",
                                             ifelse(sig68 == 1 & sig95 == 0 & mean > 0, "Increases occurrence (68% confidence)",
                                                    ifelse(sig68 == 1 & sig95 == 0 & mean < 0, "Decreases occurrence (68% confidence)", NA)))))) |> 
  dplyr::mutate( lab = 
                   factor(lab, 
                          levels = c(
                            "Decreases occurrence (95% confidence)",
                            "Decreases occurrence (68% confidence)",
                            "No effect",
                            "Increases occurrence (68% confidence)",
                            "Increases occurrence (95% confidence)"))) |> 
  dplyr::left_join(param_name_key2) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  dplyr::mutate( var_name = factor(var_name, 
                                   levels = c(
                                     "Area", "Open", "Ratio", "Precip."
                                   ))) |> 
  dplyr::mutate(common = factor(common, levels = unique(common))) |> 
  dplyr::mutate(common = factor(common, levels = rev(levels(common))))

( panela <- ggplot( occ_params, aes(x = var_name, 
                                    y = common, 
                                    fill = lab)) + 
    facet_wrap(~type, ncol = 1, scales = "free_y",
               strip.position = "right") +
    ggh4x::force_panelsizes(rows = c(0.2424242,0.7575758)) +
    geom_tile(color = "black",
              linewidth = 0.1) +
    scale_fill_manual(values = c(
      MetBrewer::MetPalettes$OKeeffe1[[1]][2],
      MetBrewer::MetPalettes$OKeeffe1[[1]][4],
      "gray90",
      MetBrewer::MetPalettes$OKeeffe1[[1]][8],
      MetBrewer::MetPalettes$OKeeffe1[[1]][10])) +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    ggtitle("(a)   occurrence") +
    theme(
      plot.title = element_text(color = "black", size = 10, hjust = -0.4, face = "bold"),
      axis.text.x = element_text(angle = 45,
                                 hjust = -0.05,
                                 vjust = -20,
                                 color = "black",
                                 size = 9),
      axis.text.y = element_text(color = "black", size = 9),
      axis.title = element_blank(),
      strip.text = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 8),
      legend.key.height = unit(1.5, "pt"),
      legend.justification = "left",
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(-5, 0, 0, -75)) +
    guides(fill = guide_legend(ncol = 1)) )

panela | panelb

setwd(here::here("figures"))
ggsave(
  "figure_04.png",
  width = 6, 
  height = 5.6, 
  units = "in", 
  dpi = 600)