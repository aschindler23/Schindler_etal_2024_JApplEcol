
library(dplyr)
library(tidyr)
library(stringr)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)
# setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Nimble resources")
source("nimble_restart_functions.R")

# setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/revised analysis 070424/data")

### load data
flock_coord <- read.csv("roost_site_coordinates.csv")
count <- read.csv("count_data.csv")

# scale flock coordinates
flock_coord <- flock_coord %>% 
  mutate(lat_scale = as.numeric(scale(lat)),
         lon_scale = as.numeric(scale(lon)))

### indexing
# number of flocks
nflocks <- length(unique(count$flock_id))

# number of groups/classes
K <- 3

# load trend model
load("flock_trends_linear_070424_results.RData")
samples_mat <- as.matrix(rbind(out[[1]][[1]], out[[2]][[1]], out[[3]][[1]]))
samples_mat <- samples_mat[seq(1, 12000, by = 4),]
mu_lambda <- samples_mat[, str_detect(colnames(samples_mat), "beta_trend")]

### LCA model
# specify model in BUGS language 
LCA_model <- nimbleCode({
  
  # Dirichlet prior for mixing proportions
  w[1:K] ~ ddirich(w_alpha[1:K])
  
  for(k in 1:K){ # loop over classes (K = 3 classes)
    # hyperprior for (equal) mixing proportions
    w_alpha[k] <- 1
    
    # LCA covariate means
    class_mu_lambda[k] ~ dnorm(0, sd = 10)
    class_mu_lat[k] ~ dnorm(0, sd = 10)
    class_mu_lon[k] ~ dnorm(0, sd = 10)
    
    # LCA covariate precision
    class_sig_lambda[k] ~ dunif(0, 10)
    class_sig_lat[k] ~ dunif(0, 10)
    class_sig_lon[k] ~ dunif(0, 10)
  } # k
  
  for(f in 1:nflocks){ # loop over flocks
    # latent class membership of each flock
    Z[f] ~ dcat(w[1:K]) 
    for(s in 1:nsim){
      mu_lambda[s, f] ~ dnorm(class_mu_lambda[Z[f]], sd = class_sig_lambda[Z[f]]) 
    }
    latitude[f] ~ dnorm(class_mu_lat[Z[f]], sd = class_sig_lat[Z[f]])
    longitude[f] ~ dnorm(class_mu_lon[Z[f]], sd = class_sig_lon[Z[f]])
  } # f
})

# set constants
nimble_constants <- list(nflocks = nflocks,
                         nsim = nrow(mu_lambda),
                         K = K)

# set data
nimble_data <- list(mu_lambda = mu_lambda, 
                    latitude = flock_coord$lat_scale,
                    longitude = flock_coord$lon_scale
                    )

# initial values
Z_init <- sample(1:K, nflocks, replace = T)
inits_function <- function(){
  list(class_mu_lambda = rnorm(K, 0, 1),
       class_sig_lambda = runif(K, 0.1, 10),
       class_mu_lat = rnorm(K, 0, 1),
       class_sig_lat = runif(K, 0.1, 10),
       class_mu_lon = rnorm(K, 0, 1),
       class_sig_lon = runif(K, 0.1, 10),
       w = c(rep(1/K, K)),
       Z = Z_init
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("LCA_model", "nimble_data", "nimble_constants", "inits", 
                    "getMCMCstate", "getModelState", "getStateVariableNames",
                    "getWAICstate", "setMCMCstate", "setModelState", "setWAICstate"))

for (j in seq_along(cl)) {
  nimble_inits <- inits_function() 
  clusterExport(cl[j], "nimble_inits")
}

# run MCMC in parallel
start_time <- Sys.time()
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(parallel)
  
  # build model
  model <- nimbleModel(code = LCA_model, constants = nimble_constants,  dat =  nimble_data, inits = nimble_inits)
  
  # initialize remaining parameters
  # model$initializeInfo()
  # model$check()
  # model$calculate()
  
  # configure MCMC
  mcmc_Conf  <- configureMCMC(model)
  mcmc_Conf$addMonitors(c("Z"))
  
  # build MCMC
  modelMCMC  <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  # run model
  CmodelMCMC$run(niter = 74000,
                 nburnin = 58000,
                 thin = 4
  )
  
  # extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)
  
  # extract model state and internal state of MCMC 
  modelState <- getModelState(Cmodel)
  mcmcState <- getMCMCstate(mcmc_Conf, CmodelMCMC)
  
  # extract R's random number seed
  seed <- .Random.seed
  
  return(list(
    samples = samples,
    modelState = modelState,
    mcmcState = mcmcState,
    seed = seed
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
file_heading <- paste0("LCA_K", K, Sys.Date())
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples    <- list(chain1 = out[[1]][[1]], 
                   chain2 = out[[2]][[1]], 
                   chain3 = out[[3]][[1]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# traceplots
MCMCtrace(mcmc_list, type = "trace", iter = 4000,
          filename = paste0(file_heading, "_traceplots.pdf"))

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]
rhat <- as.data.frame(rhat)
colnames(rhat) <- c("Rhat", "Rhat_UCI")
