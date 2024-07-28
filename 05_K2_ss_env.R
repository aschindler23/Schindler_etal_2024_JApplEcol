
library(dplyr)
library(tidyr)
library(stringr)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)
#setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Nimble resources")
source("nimble_restart_functions.R")

### load data
#setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/revised analysis 070424/data")
count <- read.csv("count_data.csv")
weather <- read.csv("weather_data.csv") %>% 
  arrange(year_id)
pland <- read.csv("pland_data.csv")
winter_GDD <- read.csv("winter_GDD_data.csv")

### format data
count_scaled <-  count %>%
  group_by(flock_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, flock_id, count) %>% 
  mutate(count = as.numeric(scale(count, center = T, scale = T))) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y") %>% 
  select(!flock_id) %>% 
  as.matrix()

grass <- pland %>% 
  select(flock_id, year_id, scaled_grass) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_grass,
              names_prefix = "y") %>% 
  arrange(flock_id) %>% 
  select(-flock_id) %>% 
  as.matrix()

cereal <- pland %>% 
  select(flock_id, year_id, scaled_cereal) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_cereal,
              names_prefix = "y") %>% 
  arrange(flock_id) %>% 
  select(-flock_id) %>% 
  as.matrix()

bog <- pland %>% 
  select(flock_id, year_id, scaled_peat_bog) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_peat_bog,
              names_prefix = "y") %>% 
  arrange(flock_id) %>% 
  select(-flock_id) %>% 
  as.matrix()

winter_GDD <- winter_GDD %>% 
  select(flock_id, year_id, scaled_winter_GDD) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_winter_GDD,
              names_prefix = "y") %>% 
  arrange(flock_id) %>% 
  select(-flock_id) %>% 
  as.matrix()

### format initial values
N_init <- count %>%
  group_by(flock_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, flock_id, count) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y") %>% 
  select(!flock_id) %>% 
  as.matrix()
N1 <- N_init[,1]
N_init <- N_init + 0.001

# drop NA counts to speed up MCMC sampling
count <- count %>% 
  select(flock_id, year_id, count) %>% 
  arrange(flock_id, year_id) %>% 
  drop_na(count)

### LCA classifications
# setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/revised analysis 070424")
Z <- read.csv("LCA_K2_2024-07-05_assignments.csv")
Z <- Z %>% 
  mutate(Z = if_else(class == "K1", 1, 2)) %>% 
  select(flock, Z) %>% 
  rename(flock_id = flock) %>% 
  arrange(flock_id)

### indexing
# maximum years
max_years <- count %>% 
  select(flock_id, year_id) %>% 
  group_by(flock_id) %>% 
  summarise(max_year = max(year_id)) %>% 
  ungroup() %>% 
  arrange(flock_id)

### indexing
# number of flocks
nflocks <- length(unique(count$flock_id))

# number of years
nyears <- length(unique(count$year_id))

# number of observations in the count data
nobs <- nrow(count)

# number of environmental covariates
nenv <- 10

# number of classes
K <- 2

### state-space model to quantify environmental drivers of abundance trends
# specify model in BUGS language 
ss_env_model <- nimbleCode({
  
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
  for(k in 1:K){
    alpha_lambda[k] ~ dnorm(0, sd = 10)
  }
  
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
        alpha_lambda[Z[f]] + 
        beta_density[f] * count[f, t] +
        beta_env[1, Z[f]] * spring_storm_days[t + 1] + 
        beta_env[2, Z[f]] * spring_precip[t + 1] + 
        beta_env[3, Z[f]] * pre_breed_freeze[t + 1] +
        beta_env[4, Z[f]] * post_breed_precip[t + 1] +
        beta_env[5, Z[f]] * autumn_storm_days[t + 1] +
        beta_env[6, Z[f]] * autumn_freeze[t + 1] +
        beta_env[7, Z[f]] * grass[f, t + 1] + 
        beta_env[8, Z[f]] * cereal[f, t + 1] + 
        beta_env[9, Z[f]] * bog[f, t + 1] +
        beta_env[10, Z[f]] * winter_GDD[f, t]
      N[f, t + 1] <- exp(log_N[f, t + 1])
    } # t
  } # f
  
  # observation model for count data
  for(i in 1:nobs){ # loop over count data observations
    y[i] ~ dpois(N[flock_id[i], year_id[i]])
  } # i
})

# set constants
nimble_constants <- list(nflocks = nflocks,
                         max_years = max_years$max_year,
                         nenv = nenv,
                         K = K, 
                         Z = Z$Z,
                         flock_id = count$flock_id,
                         year_id = count$year_id,
                         nobs = nobs
                    )

# set data
nimble_data <- list(y = count$count,
                    N1 = N1,
                    count = count_scaled,
                    spring_storm_days = weather$scaled_spring_storm_days,
                    spring_precip = weather$scaled_spring_stage_precip,
                    pre_breed_freeze = weather$scaled_pre_breed_freeze,
                    post_breed_precip = weather$scaled_post_breed_precip,
                    autumn_storm_days = weather$autumn_storm_days,
                    autumn_freeze = weather$scaled_autumn_freeze,
                    grass = grass, 
                    cereal = cereal, 
                    bog = bog,
                    winter_GDD = winter_GDD)

# initial values
inits_function <- function(){
  list(N = N_init,
       log_N = log(N_init),
       mu_lambda = matrix(rnorm(nflocks * nyears, 0, 0.1), nflocks, nyears),
       beta_env = matrix(rnorm(nenv * K, 0, 1), nenv, K),
       beta_density = rnorm(nflocks, 0, 1),
       alpha_lambda = rnorm(K, 0, 1),
       sig_lambda = runif(nflocks, 0.1, 1),
       log_lambda = matrix(rnorm(nflocks * nyears, 0, 0.1), nflocks, nyears)
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("ss_env_model", "nimble_data", "nimble_constants", "inits", 
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
  model <- nimbleModel(code = ss_env_model, constants = nimble_constants,  dat =  nimble_data, inits = nimble_inits)
  
  # initialize remaining parameters
  # model$initializeInfo()
  # model$check()
  # model$calculate()
  
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
file_heading <- paste0("K", K, "_ss_env_", Sys.Date())
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples    <- list(chain1 = out[[1]][[1]], 
                   chain2 = out[[2]][[1]], 
                   chain3 = out[[3]][[1]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# traceplots
MCMCtrace(mcmc_list, type = "trace", filename = paste0(file_heading, "_traceplots.pdf"))

for (i in 1:nc) {
  mcmc_list[[i]] <- mcmc_list[[i]][ , which(is.na(colSums(mcmc_list[[i]])) == F)]
}

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]

# summary statistics
prob_above_below_0 <- function(x){
  cdf <- ecdf(x)
  prob_below_0 <- cdf(0)
  prob_above_0 <- 1 - prob_below_0
  if(median(x) < 0){
    return(prob_below_0)
  } else{
    return(prob_above_0)
  }
}

sum_stats <- MCMCsummary(mcmc_list, func = function(x) prob_above_below_0(x), func_name = "f")
write.csv(sum_stats, paste0(file_heading, "_summary.csv"))
