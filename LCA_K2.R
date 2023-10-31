
library(dplyr)
library(tidyr)
library(stringr)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load flock counts and coordinates
count <- read.csv("flock_counts.csv")
flock_coord <- read.csv("flock_coordinates.csv")

### indexing
# number of flocks
nflocks <- length(unique(count$flock2))

# number of years
nyears <- length(unique(count$year2))

# number of groups/classes
K <- 2

# load trend model results
load("flock_trends_mod.RData")
samples_mat <- as.matrix(rbind(out[[1]], out[[2]], out[[3]]))
samples_mat <- samples_mat[seq(1, 12000, by = 2),]
trend <- samples_mat[, str_detect(colnames(samples_mat), "beta_trend")]

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
    
    # LCA covariate precision
    class_sig_lambda[k] ~ dunif(0, 10)
    class_sig_lat[k] ~ dunif(0, 10)
  } # k
  
  for(f in 1:nflocks){ # loop over flocks
    # latent class membership of each flock
    Z[f] ~ dcat(w[1:K]) 
    for(s in 1:nsim){
      trend[s, f] ~ dnorm(class_mu_lambda[Z[f]], sd = class_sig_lambda[Z[f]]) 
    }
    latitude[f] ~ dnorm(class_mu_lat[Z[f]], sd = class_sig_lat[Z[f]])
  } # f
})

# set constants
nimble_constants <- list(nflocks = nflocks,
                         nsim = nrow(trend),
                         K = K)

# set data
nimble_data <- list(trend = trend, latitude = flock_coord$lat_scale)

# initial values
Z_init <- sample(1:K, nflocks, replace = T)
inits_function <- function(){
  list(class_mu_lambda = rnorm(K, 0, 1),
       class_sig_lambda = runif(K, 0.1, 10),
       class_mu_lat = rnorm(K, 0, 1),
       class_sig_lat = runif(K, 0.1, 10),
       w = c(rep(1/K, K)),
       Z = Z_init
  )
}

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("LCA_model", "nimble_data", "nimble_constants"))

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
  model$initializeInfo()
  model$calculate()
  
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
  
  return(as.matrix(CmodelMCMC$mvSamples))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
save(out, file = "LCA_K2_results.RData")

# convert to MCMC list
samples <- list(chain1 = out[[1]], 
                chain2 = out[[2]], 
                chain3 = out[[3]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# traceplots
MCMCtrace(mcmc_list, type = "trace", filename = "LCA_K2_traceplots.pdf")

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]

sum_stats <- MCMCsummary(mcmc_list)
write.csv(sum_stats, "LCA_K2_summary.csv")

