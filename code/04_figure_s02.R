# create figure 3
# visualize community-level and species-level covariate effects and influence on trend
library(here)
library(tidyverse)
library(patchwork)
library(MCMCvis)
library(MetBrewer)
library(reshape2)

anomaly <- readr::read_csv(here::here("data/weather_anomaly.csv")) |> 
  dplyr::rename(point = siteID)

hab <- readr::read_csv(here::here("data/habitat_vars.csv"))

# read in bird data; contains columns for the focal covariates as well
d <- readr::read_csv(here::here("data/mn_prairie_bird_data_clean.csv")) |> 
  dplyr::filter(!is.na(site)) |> 
  dplyr::left_join(anomaly) |> 
  dplyr::left_join(hab)

sp_to_select <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::group_by(sp) |> 
  dplyr::summarise(prop = sum( n > 0) / sum(!is.na(n))) |> 
  dplyr::filter(prop >= 0.03) |> 
  # omit swallows (bars and tres) and ducks (mall) which are mostly flyovers, and unknown (unk) birds
  dplyr::filter(!sp %in% c("bars", "mall", "tres", "unk")) |> 
  dplyr::pull(sp)

z_info <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::filter(sp %in% sp_to_select) |> 
  dplyr::select(site, point, year, pgrass_250, pgrass_3000, ed_250, ed_3000, tanom_nb, tanom_br, panom_nb, panom_br, sp) |> 
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

year.sc <- scale(z_info$year)

dat_tab <- tibble(
  habitat = c(   "pgrass_250", "pgrass_250", "pgrass_3000", "pgrass_3000", "ed_250",   "ed_250",  "ed_3000",  "ed_3000"), 
  tmp_season = c("tanom_br",    "tanom_nb",   "tanom_br",   "tanom_nb",    "tanom_br", "tanom_nb", "tanom_br","tanom_nb"),
  prcp_season = c("panom_br",    "panom_nb",   "panom_br",   "panom_nb",    "panom_br", "panom_nb", "panom_br","panom_nb")) |> 
  dplyr::mutate(rowid = dplyr::row_number())

star_list <- list(list())
fit_sp_list <- list(list())
fit_trnd_sp_list <- list(list())
fit_com_list <- list(list())
fit_trnd_com_list <- list(list())

setwd(here::here("results"))
for( i in 1:nrow(dat_tab)){

  if(!file.exists(paste0("0", i, "_", dat_tab[[i, 'habitat']], "_",
                       dat_tab[[i, 'tmp_season']], "_",
                       dat_tab[[i, 'prcp_season']], ".RData"))) {
    next
  }
  
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
  
  post_com <- MCMCvis::MCMCpstr( out, params = "mu_beta", type = "chains")[[1]] |>
    reshape2::melt(c("param", "iter")) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(param = readr::parse_number(as.character(param))) |> 
    dplyr::left_join(param_key) |> 
    dplyr::select(-param) |> 
    tidyr::pivot_wider(names_from = name, values_from = value)
  
  # classify "significance" of parameters
  stars <- post_com |> 
    # dplyr::group_by(sp) |> 
    dplyr::select( hab_name, tmp_season, prcp_season, iter, yr:habitat.yr) |> 
    tidyr::pivot_longer(yr:habitat.yr, names_to = "var", values_to = "val") |> 
    dplyr::group_by(hab_name, tmp_season, prcp_season, var) |> 
    dplyr::summarise( l95 = quantile(val, c(0.025)), 
                      u95 = quantile(val, c(0.975)), 
                      l68 = quantile(val, c(0.160)), 
                      u68 = quantile(val, c(0.840)),
                      mean = mean(val)) |> 
    dplyr::mutate(sig68 = ifelse(l68 < 0 & u68 > 0, 0, 1), 
                  sig95 = ifelse(l95 < 0 & u95 > 0, 0, 1)) |> 
    dplyr::mutate(lab = ifelse( sig95 == 1, "**", 
                                ifelse(sig95 == 0 & sig68 == 1, "*", ""))) |> 
    dplyr::select(hab_name, tmp_season, prcp_season, var, sig68, sig95, mean, l68, u68, l95, u95) |> 
    dplyr::mutate(var = str_replace_all(var, "habitat", hab_name)) #|> 
  # dplyr::filter(grepl(hab_name, var))
  
  pdat <- tidyr::expand_grid(
    habitat = seq(from = min(data$x[,5]), to = max(data$x[,5]), length.out = 20))
  
  hab.sc <- z_info |> 
    dplyr::select(dat_tab[[i, 'habitat']]) |> 
    dplyr::pull() |> 
    scale()
  
  fit_com <- dplyr::cross_join(
    dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, habitat), 
    pdat) |> 
    dplyr::mutate(psi = plogis( intercept +  (habitat.x * habitat.y))) |> 
    dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = habitat.y, psi ) |> 
    dplyr::group_by( hab_name, tmp_season, prcp_season, habitat) |> 
    dplyr::summarise( mean = mean(psi), 
                      l95 = quantile(psi, c(0.025)), 
                      l68 = quantile(psi, c(0.160)), 
                      u68 = quantile(psi, c(0.840)), 
                      u95 = quantile(psi, c(0.975))) |> 
    dplyr::mutate( habitat = habitat*attr(hab.sc, "scaled:scale") + attr(hab.sc, "scaled:center")) |>
    tibble::add_column( cov_name = dat_tab[[i, 'habitat']])
  
  fit_sp <- dplyr::cross_join(
    dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, habitat), 
    pdat) |> 
    dplyr::mutate(psi = plogis( intercept +  (habitat.x * habitat.y))) |> 
    dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = habitat.y, psi ) |> 
    dplyr::group_by( sp, hab_name, tmp_season, prcp_season, habitat) |> 
    dplyr::summarise( mean = mean(psi), 
                      l95 = quantile(psi, c(0.025)), 
                      l68 = quantile(psi, c(0.160)), 
                      u68 = quantile(psi, c(0.840)), 
                      u95 = quantile(psi, c(0.975))) |> 
    dplyr::mutate( habitat = habitat*attr(hab.sc, "scaled:scale") + attr(hab.sc, "scaled:center")) |>
    tibble::add_column( cov_name = dat_tab[[i, 'habitat']])
  
  pdat_trnd <- tidyr::expand_grid(
    yr = c(2008, 2023),
    habitat = seq(from = min(data$x[,5]), to = max(data$x[,5]), length.out = 20)) |> 
    dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))
  
  fit_trnd_com <- dplyr::cross_join(
    dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, habitat, yr, habitat.yr), 
    pdat_trnd) |> 
    dplyr::mutate(psi = plogis( intercept +  (habitat.x * habitat.y) + (yr.x * yr.y) + (habitat.yr*yr.y*habitat.y))) |> 
    dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = habitat.y, yr = yr.y, psi ) |> 
    dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
    tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
    dplyr::mutate(diff = `2023` - `2008`) |> 
    dplyr::group_by(hab_name, tmp_season, prcp_season, habitat) |> 
    dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
    dplyr::mutate( habitat = habitat*attr(hab.sc, "scaled:scale") + attr(hab.sc, "scaled:center")) 
  
  fit_trnd_sp <- dplyr::cross_join(
    dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, habitat, yr, habitat.yr), 
    pdat_trnd) |> 
    dplyr::mutate(psi = plogis( intercept +  (habitat.x * habitat.y) + (yr.x * yr.y) + (habitat.yr*yr.y*habitat.y))) |> 
    dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = habitat.y, yr = yr.y, psi ) |> 
    dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
    tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
    dplyr::mutate(diff = `2023` - `2008`) |> 
    dplyr::group_by(sp, hab_name, tmp_season, prcp_season, habitat) |> 
    dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
    dplyr::mutate( habitat = habitat*attr(hab.sc, "scaled:scale") + attr(hab.sc, "scaled:center")) 
  
  if(i == 1 | i == 2){
    tanom.sc <- z_info |> 
      dplyr::select(dat_tab[[i, 'tmp_season']]) |> 
      dplyr::pull() |> 
      scale()
  
    panom.sc <- z_info |> 
      dplyr::select(dat_tab[[i, 'prcp_season']]) |> 
      dplyr::pull() |> 
      scale()
    
    pdat_tanom <- tidyr::expand_grid(
      tanom = seq(from = min(tanom.sc), to = max(tanom.sc), length.out = 20))
    
    pdat_panom <- tidyr::expand_grid(
      panom = seq(from = min(panom.sc), to = max(panom.sc), length.out = 20))
    
    fit_com_tanom <- dplyr::cross_join(
      dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, tanom), 
      pdat_tanom) |> 
      dplyr::mutate(psi = plogis( intercept +  (tanom.x * tanom.y))) |> 
      dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = tanom.y, psi ) |> 
      dplyr::group_by( hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( mean = mean(psi), 
                        l95 = quantile(psi, c(0.025)), 
                        l68 = quantile(psi, c(0.160)), 
                        u68 = quantile(psi, c(0.840)), 
                        u95 = quantile(psi, c(0.975))) |> 
      dplyr::mutate( habitat = habitat*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |>
      mutate(hab_name = "tanom")
    
    fit_com_panom <- dplyr::cross_join(
      dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, panom), 
      pdat_panom) |> 
      dplyr::mutate(psi = plogis( intercept +  (panom.x * panom.y))) |> 
      dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = panom.y, psi ) |> 
      dplyr::group_by( hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( mean = mean(psi), 
                        l95 = quantile(psi, c(0.025)), 
                        l68 = quantile(psi, c(0.160)), 
                        u68 = quantile(psi, c(0.840)), 
                        u95 = quantile(psi, c(0.975))) |> 
      dplyr::mutate( habitat = habitat*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |>
      mutate(hab_name = "panom")
    
    fit_sp_tanom <- dplyr::cross_join(
      dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, tanom), 
      pdat_tanom) |> 
      dplyr::mutate(psi = plogis( intercept +  (tanom.x * tanom.y))) |> 
      dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = tanom.y, psi ) |> 
      dplyr::group_by( sp, hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( mean = mean(psi), 
                        l95 = quantile(psi, c(0.025)), 
                        l68 = quantile(psi, c(0.160)), 
                        u68 = quantile(psi, c(0.840)), 
                        u95 = quantile(psi, c(0.975))) |> 
      dplyr::mutate( habitat = habitat*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |>
      mutate(hab_name = "tanom")
    
    fit_sp_panom <- dplyr::cross_join(
      dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, panom), 
      pdat_panom) |> 
      dplyr::mutate(psi = plogis( intercept +  (panom.x * panom.y))) |> 
      dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = panom.y, psi ) |> 
      dplyr::group_by( sp, hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( mean = mean(psi), 
                        l95 = quantile(psi, c(0.025)), 
                        l68 = quantile(psi, c(0.160)), 
                        u68 = quantile(psi, c(0.840)), 
                        u95 = quantile(psi, c(0.975))) |> 
      dplyr::mutate( habitat = habitat*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |>
      mutate(hab_name = "panom")
    
    pdat_trnd_tanom <- tidyr::expand_grid(
      yr = c(2008, 2023),
      tanom = seq(from = min(tanom.sc), to = max(tanom.sc), length.out = 20)) |> 
      dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))
    
    pdat_trnd_panom <- tidyr::expand_grid(
      yr = c(2008, 2023),
      panom = seq(from = min(panom.sc), to = max(panom.sc), length.out = 20)) |> 
      dplyr::mutate( yr = ( yr - attr(year.sc, "scaled:center") ) / attr(year.sc, "scaled:scale"))
    
    fit_trnd_com_tanom <- dplyr::cross_join(
      dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, tanom, yr, tanom.yr), 
      pdat_trnd_tanom) |> 
      dplyr::mutate(psi = plogis( intercept +  (tanom.x * tanom.y) + (yr.x * yr.y) + (tanom.yr*yr.y*tanom.y))) |> 
      dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = tanom.y, yr = yr.y, psi ) |> 
      dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
      tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
      dplyr::mutate(diff = `2023` - `2008`) |> 
      dplyr::group_by(hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
      dplyr::mutate( habitat = habitat*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |> 
      dplyr::mutate(hab_name = "tanom")
    
    fit_trnd_com_panom <- dplyr::cross_join(
      dplyr::select(post_com, iter, hab_name, tmp_season, prcp_season, intercept, panom, yr, panom.yr), 
      pdat_trnd_panom) |> 
      dplyr::mutate(psi = plogis( intercept +  (panom.x * panom.y) + (yr.x * yr.y) + (panom.yr*yr.y*panom.y))) |> 
      dplyr::select(hab_name, tmp_season, prcp_season, iter, habitat = panom.y, yr = yr.y, psi ) |> 
      dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
      tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
      dplyr::mutate(diff = `2023` - `2008`) |> 
      dplyr::group_by(hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
      dplyr::mutate( habitat = habitat*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |> 
      dplyr::mutate(hab_name = "panom")
    
    fit_trnd_sp_tanom <- dplyr::cross_join(
      dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, tanom, yr, tanom.yr), 
      pdat_trnd_tanom) |> 
      dplyr::mutate(psi = plogis( intercept +  (tanom.x * tanom.y) + (yr.x * yr.y) + (tanom.yr*yr.y*tanom.y))) |> 
      dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = tanom.y, yr = yr.y, psi ) |> 
      dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
      tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
      dplyr::mutate(diff = `2023` - `2008`) |> 
      dplyr::group_by(sp, hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
      dplyr::mutate( habitat = habitat*attr(tanom.sc, "scaled:scale") + attr(tanom.sc, "scaled:center")) |> 
      dplyr::mutate(hab_name = "tanom")
    
    fit_trnd_sp_panom <- dplyr::cross_join(
      dplyr::select(post_sp, sp, iter, hab_name, tmp_season, prcp_season, intercept, panom, yr, panom.yr), 
      pdat_trnd_panom) |> 
      dplyr::mutate(psi = plogis( intercept +  (panom.x * panom.y) + (yr.x * yr.y) + (panom.yr*yr.y*panom.y))) |> 
      dplyr::select(sp, hab_name, tmp_season, prcp_season, iter, habitat = panom.y, yr = yr.y, psi ) |> 
      dplyr::mutate(yr = ifelse(yr < 0, 2008, 2023)) |> 
      tidyr::pivot_wider(names_from = yr, values_from = psi) |> 
      dplyr::mutate(diff = `2023` - `2008`) |> 
      dplyr::group_by(sp, hab_name, tmp_season, prcp_season, habitat) |> 
      dplyr::summarise( test = sum( diff < 0 ) / n()) |> 
      dplyr::mutate( habitat = habitat*attr(panom.sc, "scaled:scale") + attr(panom.sc, "scaled:center")) |> 
      dplyr::mutate(hab_name = "panom")
    
    star_list[[i]] <- stars
    
    fit_sp_list[[i]] <- full_join(fit_sp, fit_sp_panom) |> 
      full_join(fit_sp_tanom)
    
    fit_trnd_sp_list[[i]] <- full_join(fit_trnd_sp, fit_trnd_sp_panom) |> 
      full_join(fit_trnd_sp_tanom)
    
    fit_com_list[[i]] <- full_join(fit_com, fit_com_tanom) |> 
      full_join(fit_com_panom)
    
    fit_trnd_com_list[[i]] <- full_join(fit_trnd_com, fit_trnd_com_panom) |> 
      full_join(fit_trnd_com_tanom)
    
  } else {
    star_list[[i]] <- stars
    fit_sp_list[[i]] <- fit_sp
    fit_trnd_sp_list[[i]] <- fit_trnd_sp
    fit_com_list[[i]] <- fit_com
    fit_trnd_com_list[[i]] <- fit_trnd_com
  }
  print(paste('finished', i, "of", nrow(dat_tab)))
}

name_key <- tibble::tribble(
  ~hab_name, ~var_lab,
  "ed_250", "Edge",
  "ed_3000", "Edge density (3000 m)", 
  "panom", "Precip.", 
  "pgrass_250", "Grassland",
  "pgrass_3000", "Grassland (3000 m)",
  "tanom", "Temp.")

fit_sp <- dplyr::bind_rows(fit_sp_list) |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
  dplyr::filter(!grepl("3000", hab_name)) |>
  dplyr::left_join(name_key) |> 
  dplyr::mutate(var_lab = factor(var_lab, 
                                 levels = c(
                                   "Grassland",
                                   "Edge",
                                   "Temp.",
                                   "Precip.")))

fit_com <- dplyr::bind_rows(fit_com_list) |> 
  dplyr::filter(tmp_season == "tanom_br") |> 
  dplyr::filter(!grepl("3000", hab_name)) |>
  dplyr::left_join(name_key) |> 
  dplyr::mutate(var_lab = factor(var_lab, 
                                 levels = c(
                                   "Grassland",
                                   "Edge",
                                   "Temp.",
                                   "Precip.")))

xpos <- fit_com |> 
  dplyr::group_by(var_lab) |> 
  dplyr::filter(habitat == min(habitat) | habitat == max(habitat)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(xpos = ifelse(habitat == 0, 0.05*max(habitat), habitat + (0.05*abs(habitat)) )) |> 
  dplyr::select(var_lab, habitat = xpos) |> 
  dplyr::group_by(var_lab) |> 
  dplyr::filter(habitat == min(habitat)) |> 
  tibble::add_column(mean = 1.05)

stars <- dplyr::bind_rows(star_list) 

stars |>
  filter(grepl("pgrass", var) | grepl("ed", var)) |>
  filter(!grepl(".yr", var)) |>
  mutate(scale = parse_number(var)) |>
  separate(var, into = c("var", "junk"), sep = "_") |>
  separate(tmp_season, into = c("junk", "season"), sep = "_") |> 
  dplyr::mutate(var = ifelse(var == "pgrass", "proportion grassland", "edge density")) |> 
  dplyr::mutate(season = ifelse(season == "br", "breeding", "nonbreeding")) |> 
  ggplot(aes(x = mean, y = var, color = factor(scale), shape = season)) +
  geom_vline(xintercept = 0, color = "gray30", linetype = "dashed") +
  geom_errorbar(aes(xmin = l95, xmax = u95, color = factor(scale)),
                position = position_dodge(width = 0.4), width = 0, linewidth = 1) +
  geom_point(position = position_dodge(width = 0.4),
             size = 3) +
  theme_classic() +
  labs(x = "assemblage-level effect",
       color = "buffer radius (m)",
       shape = "season of\nweather anomaly") +
  theme(axis.title.y = element_blank())

ggsave(
  filename = here::here("figures/figure_s02.png"),
  width = 6, 
  height = 6,
  units = "in", 
  dpi = 600
)             
