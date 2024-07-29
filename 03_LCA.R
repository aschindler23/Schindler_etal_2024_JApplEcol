
library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load data
# site coordinates
site_coord <- read.csv("roost_site_coordinates.csv")

# count data
count <- read.csv("count_data.csv")

# trend model results
load("site_trends_results.RData")
samples_mat <- as.matrix(rbind(out[[1]][[1]], out[[2]][[1]], out[[3]][[1]]))

# thin trend model results
samples_mat <- samples_mat[seq(1, 12000, by = 4),]

# site trends
beta_trend <- samples_mat[, str_detect(colnames(samples_mat), "beta_trend")]

### indexing
# number of sites
nsites <- length(unique(count$site_id))

# number of groups/classes
K <- 2
# K <- 3

# number of posterior samples for the trend betas
nsamp <- nrow(beta_trend)

### LCA model
# specify model in BUGS language 
LCA_model <- nimbleCode({
  
  # Dirichlet prior for mixing proportions
  w[1:K] ~ ddirich(w_alpha[1:K])
  
  for(k in 1:K){ # loop over classes (K = 3 classes)
    # hyperprior for (equal) mixing proportions
    w_alpha[k] <- 1
    
    # LCA covariate means
    class_mu_trend[k] ~ dnorm(0, sd = 10)
    class_mu_lat[k] ~ dnorm(0, sd = 10)
    class_mu_lon[k] ~ dnorm(0, sd = 10)
    
    # LCA covariate precision
    class_sig_trend[k] ~ dunif(0, 10)
    class_sig_lat[k] ~ dunif(0, 10)
    class_sig_lon[k] ~ dunif(0, 10)
  } # k
  
  for(i in 1:nsites){ # loop over sites
    # latent class membership of each site
    Z[i] ~ dcat(w[1:K]) 
    for(s in 1:nsamp){ # loop over trend model posterior samples
      beta_trend[s, i] ~ dnorm(class_mu_trend[Z[i]], sd = class_sig_trend[Z[i]]) 
    }
    latitude[i] ~ dnorm(class_mu_lat[Z[i]], sd = class_sig_lat[Z[i]])
    longitude[i] ~ dnorm(class_mu_lon[Z[i]], sd = class_sig_lon[Z[i]])
  } # i
})

# set constants
nimble_constants <- list(nsites = nsites,
                         nsamp = nsamp,
                         K = K)

# set data
nimble_data <- list(beta_trend = beta_trend, 
                    latitude = site_coord$lat_scale,
                    longitude = site_coord$lon_scale
                    )

# initial values
inits_function <- function(){
  list(class_mu_trend = rnorm(K, 0, 1),
       class_sig_trend = runif(K, 0.1, 10),
       class_mu_lat = rnorm(K, 0, 1),
       class_sig_lat = runif(K, 0.1, 10),
       class_mu_lon = rnorm(K, 0, 1),
       class_sig_lon = runif(K, 0.1, 10),
       w = c(rep(1/K, K)),
       Z = sample(1:K, nsites, replace = T)
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("LCA_model", "nimble_data", "nimble_constants", "inits"))

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
  
  return(list(
    samples = samples
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
file_heading <- paste0("LCA_K", K)
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
