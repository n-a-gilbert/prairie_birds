# create Table S01 with species and sample size info
library(here)
library(tidyverse)
library(sf)
library(readxl)
library(avotrex)
library(flextable)
library(officer)

data(BirdTree_trees)

tax <- BirdTree_tax |> 
  tibble::as_tibble()

key <- readr::read_csv(here::here("data/sp_key.csv"))

key2 <- key |> 
  dplyr::mutate(English = common) |> 
  dplyr::mutate(English = ifelse(English == "Clay-colored Sparrow", "Clay-coloured Sparrow",
                                 ifelse(English == "Gray Catbird", "Grey Catbird",
                                        ifelse(English == "Ring-necked Pheasant", "Common Pheasant",
                                               ifelse(English == "LeConte's Sparrow", "Le Conte's Sparrow",
                                                      English))))) |> 
  dplyr::left_join(tax) |> 
  dplyr::select(sp, code, common, Genus, Species, BLFamilyLatin) |> 
  dplyr::mutate(Genus = ifelse(common == "Wilson's Snipe", "Gallinago", Genus), 
                Species = ifelse(common == "Wilson's Snipe", "delicata", Species), 
                BLFamilyLatin = ifelse(common == "Wilson's Snipe", "Scolopacidae", BLFamilyLatin)) |> 
  dplyr::mutate(scientific = paste( Genus, Species ) ) |> 
  dplyr::select(sp, code, common, scientific, family = BLFamilyLatin)

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

# "y information" - aka the bird data and associated observation covariates
y_info <- d |> 
  dplyr::mutate(sp = tolower(sp)) |> 
  dplyr::filter(sp %in% sp_to_select) |> 
  dplyr::select(site, point, area, ratio, open, anom, year, visit, date, obs, start, sp, n) |> 
  dplyr::group_by(sp) |> 
  # create species index
  dplyr::mutate(sp_y = dplyr::cur_group_id()) |> 
  dplyr::group_by(obs) |>
  # create observer idnex
  dplyr::mutate(obs_y = dplyr::cur_group_id()) |> 
  # link the bird data to the "z state"
  dplyr::full_join(
    z_info |> 
      dplyr::select(site, point, year, sp, z_index)) |> 
  dplyr::rename(z_index_y = z_index) |> 
  dplyr::mutate(y = ifelse(n > 0, 1, 0)) |> 
  dplyr::ungroup() 

obligates <- c( "Bobolink", "Dickcissel", "Grasshopper Sparrow", "LeConte's Sparrow", "Savannah Sparrow", 
                "Sedge Wren", "Upland Sandpiper", "Western Meadowlark")

sample_info <- y_info |> 
  dplyr::group_by(sp) |> 
  dplyr::summarise( prop = round(sum(n > 0) / sum(!is.na(n)), 2)) |> 
  dplyr::full_join(
    y_info |> 
      dplyr::filter(n > 0) |> 
      dplyr::group_by(sp) |> 
      dplyr::summarise( nsite = length(unique(site)))) |> 
  dplyr::rename(code = sp) |> 
  dplyr::left_join(key2) |> 
  dplyr::mutate(code = toupper(code)) |> 
  dplyr::mutate( type = ifelse(common %in% obligates, "Obligate", "Facultative")) |> 
  dplyr::mutate( type = factor(type, levels = c("Obligate", "Facultative"))) |> 
  dplyr::arrange(type, common) |> 
  dplyr::select(Type = type, Common = common, Code = code, Scientific = scientific, Family = family, ProportionSurveys = prop, NSites = nsite)

flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable(sample_info)  
setwd(here::here("figures"))
tmp <- tempfile(fileext = ".docx")
officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)