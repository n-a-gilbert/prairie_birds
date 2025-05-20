# this script formats data & fits multispecies occupancy model using nimble
library(here)
library(tidyverse)
library(nimble)
library(coda)
library(parallel)

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

# create design matrix for covariates for occurrence
x <- z_info |> 
  dplyr::select(open, area, ratio, anom, year) |> 
  tibble::add_column(int = 1) |> 
  dplyr::mutate(anom = as.numeric(scale(anom)),
                open = as.numeric(scale(open)), 
                ratio = as.numeric(scale(ratio)),
                area = as.numeric(scale(log(area))), # log of reserve area
                year = as.numeric(scale(year))) |>
  # create columns for interaction terms
  dplyr::mutate( area.ratio = area * ratio, 
                 anom.yr = anom * year,
                 open.yr = open * year, 
                 ratio.yr = ratio * year, 
                 area.yr = area * year, 
                 area.ratio.yr = area * ratio * year ) |> 
  dplyr::select( int, year, anom, area, open, ratio, 
                 area.ratio, 
                 anom.yr, area.yr, open.yr, ratio.yr, 
                 area.ratio.yr ) |> 
  as.matrix() |> 
  unname()

# create design matrix for covariates for detection
w <- y_info |> 
  tibble::add_column(int = 1) |> 
  dplyr::select(int, date, start) |> 
  dplyr::mutate(date = as.numeric(scale(date)),
                start = as.numeric(scale(start))) |> 
  dplyr::mutate(date2 = date*date, 
                start2 = start*start) |> 
  as.matrix() |> 
  unname() 

# package up "constants" for nimble
# these are used for indexing and loop control
constants <- list(
  n_beta = ncol(x), 
  n_alpha = ncol(w),
  sp_z = z_info$sp_z,
  sp_site_z = z_info$sp_site_z,
  sp_point_z = z_info$sp_point_z,
  sp_y = y_info$sp_y,
  obs_y = y_info$obs_y,
  z_index_y = y_info$z_index_y,
  n_obs = length(unique(y_info$obs_y)),
  n_sp_site_z = length(unique(z_info$sp_site_z)),
  n_sp_point_z = length(unique(z_info$sp_point_z)),
  n_sp = length(unique(y_info$sp_y)),
  n_z = nrow(z_info), 
  n_y = nrow(y_info))

# package up data for nimble model 
data <- list(
  x = x, # design matrix of occurrence covariates
  w = w, # design matrix of detection covariates
  y = y_info$y) # vector of detection-nondetection data

# code for model in nimble
model_code <- nimble::nimbleCode({
  
  sd_epsilon ~ dexp(1) # standard deviation among observers
  sd_gamma ~ dexp(1) # standard deviation among species-site combos
  sd_kappa ~ dexp(1) # standard deviation among species-point combos
  
  # random effect of observer
  for(i in 1:n_obs){
    epsilon[i] ~ dnorm(0, sd = sd_epsilon)
  }
  
  # random effect of species-site combo
  for(i in 1:n_sp_site_z){
    gamma[i] ~ dnorm(0, sd = sd_gamma)
  }
  
  # random effect of species-point combo
  for(i in 1:n_sp_point_z){
    kappa[i] ~ dnorm(0, sd = sd_kappa)
  }
  
  # hyperparameters for occurrence variables
  for(i in 1:n_beta){
    mu_beta[i] ~ dnorm(0, sd = 1) # community mean of covariate effect
    sd_beta[i] ~ dexp(1) # standard deviation among species for each covariate effect
  }
  
  # hyperparameters for detection variables
  for(i in 1:n_alpha){
    mu_alpha[i] ~ dnorm(0, sd = 1) # community mean of covariate effect
    sd_alpha[i] ~ dexp(1) # standard deviation among species for each covariate effect
  }
  
  # species-level parameters
  for(i in 1:n_sp){ # loop through species
    for( j in 1:n_beta){ # loop through occurrence variables
      beta[i,j] ~ dnorm( mu_beta[j], sd_beta[j]) # species level occurrence coefficients
    }
    for(j in 1:n_alpha){ # loop through detection variables
      alpha[i,j] ~ dnorm( mu_alpha[j], sd_alpha[j]) # species level detection coefficients
    }
  }
  
  for( i in 1:n_z){ # loop through "z-states", aka species-point-year combinations
    z[i] ~ dbern( psi[i] ) # latent ocurrence state z, occurrence probability psi
    logit(psi[i]) <- inprod(beta[sp_z[i],1:n_beta], x[i, 1:n_beta]) + gamma[sp_site_z[i]] + kappa[sp_point_z[i]] # occurrence equation
  }
  
  for( i in 1:n_y){ # loop through vector of detection-nondetection data
    logit(p[i]) <- inprod(alpha[sp_y[i],1:n_alpha],w[i, 1:n_alpha]) + epsilon[obs_y[i]] # detection equation
    y[i] ~ dbern( p[i] * z[z_index_y[i]] )
  }
})

# function to generate random initial values based on the data and constants
inits <- function(data, constants){
  
  z_st <- cbind(y = data$y,
                z_index = constants$z_index_y) |>
    tibble::as_tibble() |> 
    dplyr::group_by(z_index) |>
    dplyr::summarise(max_y = max(y)) |> 
    dplyr::pull(max_y)
  
  mu_beta <- rnorm(constants$n_beta, 0, 0.25)
  sd_beta <- rexp(constants$n_beta, 1)
  mu_alpha <- rnorm(constants$n_alpha, 0, 0.25)
  sd_alpha <- rexp(constants$n_alpha, 1)
  
  beta <- matrix(NA, nrow = constants$n_sp, ncol = constants$n_beta)
  alpha <- matrix(NA, nrow = constants$n_sp, ncol = constants$n_alpha)
  for(i in 1:constants$n_sp){
    for(j in 1:constants$n_beta){
      beta[i,j] <- rnorm(1, mu_beta[j], sd_beta[j])
    }
    for(j in 1:constants$n_alpha){
      alpha[i,j] <- rnorm(1, mu_alpha[j], sd_alpha[j])
    }
  }
  
  sd_gamma <- rexp(1)
  sd_epsilon <- rexp(1)
  sd_kappa <- rexp(1)
  
  init_list <- list(
    mu_beta = mu_beta,
    sd_beta = sd_beta, 
    mu_alpha = mu_alpha,
    sd_alpha = sd_alpha,
    sd_gamma = sd_gamma, 
    sd_epsilon = sd_epsilon,
    sd_kappa = sd_kappa,
    epsilon = rnorm(constants$n_obs, 0, sd_epsilon),
    gamma = rnorm(constants$n_sp_site_z, 0, sd_gamma),
    kappa = rnorm(constants$n_sp_point_z, 0, sd_kappa),
    beta = beta,
    alpha = alpha,
    z = z_st
  )
  return(init_list)
}

# parameters to monitor
params <- c(
  "mu_beta", "sd_beta", 
  "mu_alpha", "sd_alpha",
  "sd_gamma", "sd_epsilon",
  "sd_kappa",
  "beta", "alpha")

# I did the data formatting/prep stuff locally and then fit the model on a supercomputer, so I saved this stuff to ssh it over
# took maybe 1-2 days on the supercompueter
# would probably be fine to run it on the average desktop though
# setwd(here::here("data"))
# save(
#   constants, 
#   data, 
#   model_code, 
#   params,
#   inits,
#   file = "mn_prairie_2025_02_26.RData")

nc <- 3 # number of chains
nburn <- 30000 # number of burn-in iterations 
ni <- nburn + 20000 # total iterations
nt <- 10 # thinning interval

# running chains in parallel for efficiency
start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model_code",
                              "inits", 
                              "data", 
                              "constants", 
                              "params", 
                              "nburn", 
                              "ni", 
                              "nt"))

for(j in seq_along(cl)) {
  set.seed(j)
  init <- inits(data, constants)
  parallel::clusterExport(cl[j], "init")
}

out <- parallel::clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model <- nimble::nimbleModel(code = model_code, 
                               name = "model_code", 
                               constants = constants, 
                               data = data, 
                               inits = init)
  
  Cmodel <- nimble::compileNimble(model)
  modelConf <- nimble::configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- nimble::buildMCMC(modelConf)
  CmodelMCMC <- nimble::compileNimble(modelMCMC, project = model)
  out1 <- nimble::runMCMC(CmodelMCMC, 
                          nburnin = nburn, 
                          niter = ni, 
                          thin = nt)
  
  return(coda::as.mcmc(out1))
})
end <- Sys.time()
print(end - start)
parallel::stopCluster(cl)

# save model settings
settings <- list(
  n.chains = nc,
  n.iterations = ni,
  n.burnin = nburn,
  n.thin = nt)

# save code, data, constants, model output, and model settings
setwd(here::here("results"))
save( model_code,
      data,
      constants,
      out, 
      settings,
      file = paste0("prairie_bird_trends_global", Sys.Date(), ".RData"))