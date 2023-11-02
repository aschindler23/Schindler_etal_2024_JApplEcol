
library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load and process data
count <- read.csv("flock_counts.csv")
weather <- read.csv("weather_data.csv")
land <- read.csv("land_cover_data.csv")

# scaled count data for density dependence
count_scaled <- count %>% 
  arrange(year, flock) %>% 
  mutate(count = as.numeric(scale(count))) %>% 
  pivot_wider(names_from = year, values_from = count) %>% 
  # repace NAs with 0 (indexing prevents inclusion in model)
  mutate_all(~replace_na(.,0)) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()

# arrange land cover data
grass <-land %>% 
  arrange(year, flock) %>% 
  select(year, flock, grass_scaled) %>% 
  pivot_wider(names_from = year, values_from = grass_scaled) %>% 
  arrange(flock) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()
cereal <-land %>% 
  arrange(year, flock) %>% 
  select(year, flock, cereal_scaled) %>% 
  pivot_wider(names_from = year, values_from = cereal_scaled) %>% 
  arrange(flock) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()
bog <-land %>% 
  arrange(year, flock) %>% 
  select(year, flock, peat_bog_scaled) %>% 
  pivot_wider(names_from = year, values_from = peat_bog_scaled) %>% 
  arrange(flock) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()
ndvi <-land %>% 
  arrange(year, flock) %>% 
  select(year, flock, ndvi_scaled) %>% 
  pivot_wider(names_from = year, values_from = ndvi_scaled) %>% 
  arrange(flock) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()

# initial values for population size
N_init <- count %>% 
  select(year, flock, count) %>% 
  pivot_wider(names_from = year, values_from = count,) %>% 
  mutate_all(~replace_na(.,0)) %>% 
  column_to_rownames(var = "flock") %>% 
  as.matrix()

### indexing
# number of flocks
nflocks <- length(unique(count$flock))

# number of years
nyears <- length(unique(count$year))

# max year
max_years <- count %>% 
  select(flock, year) %>% 
  group_by(flock) %>% 
  summarise(max_year = max(year)) %>% 
  ungroup() %>% 
  arrange(flock) %>% 
  pull(max_year)

# number of environmental covariates
nenv <- 10

# number of classes
K <- 2


### LCA classifications
# setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/LCA model/trend model approach/LCA only")
Z <- read.csv("LCA_K2_assignments.csv")
Z <- Z %>% 
  mutate(Z = if_else(class == "K1", 1, 2)) %>% 
  select(flock, Z) %>% 
  arrange(flock)

### estimate effects of environmental variables on annual lambda
# specify model in BUGS language 
env_model <- nimbleCode({
  
  #===========================
  # prior specification
  #===========================
  
  for(f in 1:nflocks){ # loop over flocks
    # initial population size
    N[f, 1] <- N1[f]
    log_N[f, 1] <- log(N[f, 1])
    
    # process variance for lambda
    sig_lambda[f] ~ dunif(0, 10)
    
    # flock-specific density dependence
    beta_density[f] ~ dnorm(0, sd = 10)
  } # f
  
  # intercept for lambda
  alpha_lambda ~ dnorm(0, sd = 10)
  
  # effects of environment on lambda
  for(e in 1:nenv){ # loop over environmental variables
    for(k in 1:K){ # loop over classes
      beta_env[e, k] ~ dnorm(0, sd = 10)
    } # e
  } # k
  
  #===========================
  # likelihood specification
  #===========================
  
  # process model for true abundance
  for(f in 1:nflocks){ # loop over flocks
    for(t in 1:(max_years[f] - 1)){ # loop over years
      log_N[f, t + 1] <- log_N[f, t] + log_lambda[f, t]
      log_lambda[f, t] ~ dnorm(mu_lambda[f, t], sd = sig_lambda[f])
      mu_lambda[f, t] <- 
        alpha_lambda + 
        beta_density[f] * count[f, t] +
        beta_env[1, Z[f]] * spring_storm_days[t + 1] + 
        beta_env[2, Z[f]] * spring_precip[t + 1] + 
        beta_env[3, Z[f]] * pre_breed_freeze[t + 1] +
        beta_env[4, Z[f]] * post_breed_precip[t + 1] +
        beta_env[5, Z[f]] * fall_storm_days[t + 1] +
        beta_env[6, Z[f]] * fall_freeze[t + 1] +
        beta_env[7, Z[f]] * grass[f, t + 1] + 
        beta_env[8, Z[f]] * cereal[f, t + 1] + 
        beta_env[9, Z[f]] * bog[f, t + 1] +
        beta_env[10, Z[f]] * ndvi[f, t + 1]
      N[f, t + 1] <- exp(log_N[f, t + 1])
    } # t
  } # f
    
  # observation model for count data
  for(n in 1:ncount){ # loop over count data
    y[n] ~ dpois(N[flock[n], year[n]]) 
  }
})

# set constants
nimble_constants <- list(N1 = N_init[,1],
                         max_years = max_years,
                         nflocks = nflocks,
                         flock = count$flock,
                         year = count$year,
                         ncount = nrow(count),
                         K = K, 
                         Z = Z$Z
                    )

# set data
nimble_data <- list(y = count$count,
                    count = count_scaled,
                    spring_storm_days = weather$spring_storm_days_scaled,
                    spring_precip = weather$spring_precip,
                    pre_breed_freeze = weather$pre_breed_freeze,
                    post_breed_precip = weather$post_breed_precip,
                    fall_storm_days = weather$fall_storm_days,
                    fall_freeze = weather$fall_freeze,
                    grass = grass, 
                    cereal = cereal, 
                    bog = bog,
                    ndvi = ndvi)

# initial values
inits_function <- function(){
  list(N = N_init,
       log_N = log(N_init + 0.001),
       mu_lambda = matrix(rnorm(nflocks * nyears, 0, 0.1), nflocks, nyears),
       beta_env = matrix(rnorm(nenv * K, 0, 1), nenv, K),
       beta_density = rnorm(nflocks, 0, 1),
       alpha_lambda = rnorm(1, 0, 1),
       sig_lambda = runif(nflocks, 0.1, 1),
       log_lambda = matrix(rnorm(nflocks * nyears, 0, 0.1), nflocks, nyears)
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("env_model", "nimble_data", "nimble_constants"))

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
  model <- nimbleModel(code = env_model, constants = nimble_constants,  dat =  nimble_data, inits = nimble_inits)
  
  # initialize remaining parameters
  model$initializeInfo()
  model$calculate()
  
  # configure MCMC
  mcmc_Conf  <- configureMCMC(model)
  mcmc_Conf$addMonitors(c("log_lambda"))
  
  # build MCMC
  modelMCMC  <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  # run model
  CmodelMCMC$run(niter = 440000,
                 nburnin = 400000,
                 thin = 10
  )
  
  return(as.matrix(CmodelMCMC$mvSamples))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
file_heading <- paste0("env_K", K, "_")
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples <- list(chain1 = out[[1]], 
                chain2 = out[[2]], 
                chain3 = out[[3]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# traceplots
MCMCtrace(mcmc_list, type = "trace", filename = paste0(file_heading, "_traceplots.pdf"))

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]

# summary statistics
sum_stats <- MCMCsummary(mcmc_list)
write.csv(sum_stats, paste0(file_heading, "_summary.csv"))
