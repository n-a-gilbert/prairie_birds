# figure 4: visualize species-level parameters
# make a heatmap for main effects of covariates on occurrence and influence on trend

library(here)
library(tidyverse)
library(MCMCvis)
library(patchwork)
library(reshape2)
library(MetBrewer)
library(ggh4x)

sp_key <- readr::read_csv(here::here("data/sp_key.csv"))

setwd(here::here("results"))

dat_tab <- tibble(
  habitat = c(   "pgrass_250", "pgrass_250", "pgrass_3000", "pgrass_3000", "ed_250",   "ed_250",  "ed_3000",  "ed_3000"), 
  tmp_season = c("tanom_br",    "tanom_nb",   "tanom_br",   "tanom_nb",    "tanom_br", "tanom_nb", "tanom_br","tanom_nb"),
  prcp_season = c("panom_br",    "panom_nb",   "panom_br",   "panom_nb",    "panom_br", "panom_nb", "panom_br","panom_nb")) |> 
  dplyr::mutate(rowid = dplyr::row_number())

res <- list(list())
for( i in 1:nrow(dat_tab)){
  load(paste0("0", i, "_", dat_tab[[i, 'habitat']], "_",
              dat_tab[[i, 'tmp_season']], "_",
              dat_tab[[i, 'prcp_season']], ".RData"))
  
  param_key <- tibble(
    param = 1:8, 
    name = c("intercept", "yr", "tanom", "panom", "habitat",
             "tanom.yr", "panom.yr", "habitat.yr")) |> 
    tibble::add_column(hab_name = dat_tab[[i, 'habitat']],
                       tmp_season = dat_tab[[i, 'tmp_season']],
                       prcp_season = dat_tab[[i, 'prcp_season']])
  
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
    dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, yr:habitat.yr) |> 
    tidyr::pivot_longer(yr:habitat.yr, names_to = "var", values_to = "val") |> 
    dplyr::group_by(sp, hab_name, tmp_season, prcp_season, var) |> 
    dplyr::summarise( l95 = quantile(val, c(0.025)), 
                      u95 = quantile(val, c(0.975)), 
                      l68 = quantile(val, c(0.160)), 
                      u68 = quantile(val, c(0.840)),
                      mean = mean(val)) |> 
    dplyr::mutate(sig68 = ifelse(l68 < 0 & u68 > 0, 0, 1), 
                  sig95 = ifelse(l95 < 0 & u95 > 0, 0, 1)) |> 
    dplyr::mutate(lab = ifelse( sig95 == 1, "**", 
                                ifelse(sig95 == 0 & sig68 == 1, "*", ""))) |> 
    dplyr::select(sp, hab_name, tmp_season, prcp_season, var, sig68, sig95, mean) |> 
    dplyr::left_join(sp_key)
  
  res[[i]] <- stars
}

all <- bind_rows(res)

obligates <- c( "Bobolink", "Dickcissel", "Grasshopper Sparrow", "LeConte's Sparrow", "Savannah Sparrow", 
                "Sedge Wren", "Upland Sandpiper", "Western Meadowlark")

weather <- all |> 
  dplyr::filter(hab_name == "pgrass_250") |> 
  dplyr::filter(var == "panom" | var == "tanom") |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
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
  # dplyr::left_join(param_name_key2) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  # dplyr::mutate( var_name = factor(var_name, 
  # levels = c(
  # "Area", "Open", "Ratio", "Temp.", "Precip."
  # ))) |> 
  dplyr::mutate(common = factor(common, levels = unique(common))) |> 
  dplyr::mutate(common = factor(common, levels = rev(levels(common))))

name_key <- tibble::tribble(
  ~var, ~var_lab,
  "ed_250", "Edge",
  "ed_3000", "Edge density (3000 m)", 
  "panom", "Precip.", 
  "pgrass_250", "Grassland",
  "pgrass_3000", "Grassland (3000 m)",
  "tanom", "Temp.")

tmp <- all |> 
  dplyr::filter(var == "habitat") |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
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
  # dplyr::left_join(param_name_key2) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  # dplyr::mutate( var_name = factor(var_name, 
  # levels = c(
  # "Area", "Open", "Ratio", "Temp.", "Precip."
  # ))) |> 
  dplyr::mutate(common = factor(common, levels = unique(common))) |> 
  dplyr::mutate(common = factor(common, levels = rev(levels(common)))) |> 
  dplyr::mutate(var = hab_name) |> 
  full_join(weather) |> 
  full_join(name_key) |> 
  dplyr::mutate(var_lab = factor(var_lab,
                                 levels = c(
                                   "Grassland",
                                   "Grassland (3000 m)",
                                   "Edge",
                                   "Edge density (3000 m)",
                                   "Temp.",
                                   "Precip."))) |>
  dplyr::filter(!grepl("3000 m", var_lab))


( panela <- ggplot( tmp, aes(x = var_lab, 
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
    ggtitle("(a)   Occurrence") +
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

trnd <- all |> 
  dplyr::filter(var %in% c("habitat.yr", "panom.yr", "tanom.yr")) |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
  dplyr::mutate(var = ifelse(var == "panom.yr", "panom", 
                             ifelse(var == "tanom.yr", "tanom", hab_name))) |> 
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
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  dplyr::mutate(common = factor(common, levels = unique(common))) |> 
  dplyr::mutate(common = factor(common, levels = rev(levels(common)))) |> 
  dplyr::full_join(name_key) |> 
  dplyr::mutate(var_lab = factor(var_lab,
                                 levels = c(
                                   "Grassland",
                                   "Grassland (3000 m)",
                                   "Edge",
                                   "Edge density (3000 m)",
                                   "Temp.",
                                   "Precip."))) |>
  dplyr::filter(!grepl("3000 m", var_lab))

( panelb <- ggplot( trnd, aes(x = var_lab, 
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
    ggtitle("(b)   Trends") +
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

panela | panelb

setwd(here::here("figures"))
ggsave(
  "figure_04.png",
  width = 6, 
  height = 5.6, 
  units = "in", 
  dpi = 600)

all |> 
  dplyr::filter(var == "habitat") |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
  dplyr::mutate(group = ifelse(common %in% obligates, "obligate", "facultative")) |> 
  group_by(hab_name, group) |> 
  summarise(sig_pos = sum( sig95 == 1 & mean > 0),
            sig_neg = sum( sig95 == 1 & mean < 0)) |> 
  mutate(tot_sig = sig_pos + sig_neg) |> 
  mutate(scale = parse_number(hab_name)) |> 
  mutate(hab_name = str_remove_all(hab_name, "_250")) |> 
  mutate(hab_name = str_remove_all(hab_name, "_3000")) |> 
  dplyr::select(hab_name, scale, group, tot_sig) |> 
  pivot_wider(names_from = scale, values_from = tot_sig)
  
