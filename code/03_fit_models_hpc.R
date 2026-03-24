library(parallel)
library(nimble)
library(coda)

load("list_of_data_lists3.RData")

for(i in 1:length(dat_list)){
  
  model_code <- dat_list[[i]][[3]]
  data <- dat_list[[i]][[2]]
  constants <- dat_list[[i]][[1]]
  params <- dat_list[[i]][[4]]
  inits <- dat_list[[i]][[5]]
  results_file <- dat_list[[i]][[6]]
  
  nc <- 3 # number of chains
  nburn <- 30000 # number of burn-in iterations 
  ni <- nburn + 20000 # total iterations
  nt <- 10 # thinning interval
  
  # ---- FIX: unique build directory per Slurm job ----
  nimbleOptions(buildDir = paste0("nimble_build_",
                                  Sys.getenv("SLURM_JOB_ID")))
  model_name <- paste0("model_",
                       Sys.getenv("SLURM_JOB_ID"),
                       "_", i)
  
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
                                "nt",
                                "model_name"))
  
  for(j in seq_along(cl)) {
    set.seed(j)
    init <- inits(data, constants)
    parallel::clusterExport(cl[j], "init")
  }
  
  out <- parallel::clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    
    model <- nimble::nimbleModel(code = model_code, 
                                 name = model_name, 
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
  
  save( model_code,
        data,
        constants,
        out, 
        settings,
        file = paste0(results_file, ".RData"))
  
print(paste0("finished ", i, " of ", length(dat_list)))

}
