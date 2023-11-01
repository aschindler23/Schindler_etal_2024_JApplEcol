
library(dplyr)
library(tidyr)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load flock count data
count <- read.csv("flock_counts.csv")

# initial values for population size
N_init <- count %>% 
  select(year2, flock2, count) %>% 
  pivot_wider(names_from = year2, values_from = count, 
              names_prefix = "y") %>% 
  mutate_all(~replace_na(.,0)) %>% 
  select(-flock2)

# initial values for mean population size
alpha_init <- count %>% 
  group_by(flock2) %>% 
  summarise(mean_count = mean(count)) %>% 
  ungroup() %>% 
  mutate(log_mean_count = log(mean_count)) %>% 
  arrange(flock2) %>% 
  pull(log_mean_count)

### indexing
# number of flocks
nflocks <- length(unique(count$flock2))

# number of years
nyears <- length(unique(count$year2))

# max year
max_years <- count %>% 
  select(flock2, year2) %>% 
  group_by(flock2) %>% 
  summarise(max_year = max(year2)) %>% 
  ungroup() %>% 
  arrange(flock2) %>% 
  pull(max_year)

# center study years around 0
year_cov <- c(1:nyears) - nyears/2 - .5

### trend model
# specify model in BUGS language 
trend_model <- nimbleCode({
  #===========================
  # priors
  #===========================
  
  for(f in 1:nflocks){ # loop over flocks
    # flock-year random intercept
    alpha[f] ~ dnorm(0, sd = 10)
    
    # flock-specific trend in abundance
    beta_trend[f] ~ dnorm(0, sd = 10)
  } # f
  
  # process variance
  log_N_tau ~ dgamma(0.01, 0.01)
  
  #===========================
  # likelihood specification
  #===========================
  
  # process model for true abundance
  for(f in 1:nflocks){ # loop over flocks
    for(t in 1:(max_years[f])){ # loop over years
      log_N[f, t] ~ dnorm(log_N_mu[f, t], log_N_tau)
      log_N_mu[f, t] <- alpha[f] + beta_trend[f] * year_cov[t]
      N[f, t] <- exp(log_N[f, t])
    } # t
  } # f
  
  # observation model for count data
  for (n in 1:ncount){
    y[n] ~ dpois(N[flock[n], year[n]]) 
  }
})

# set constants
nimble_constants <- list(max_years = max_years,
                         nflocks = nflocks,
                         flock = count$flock2,
                         year = count$year2,
                         ncount = nrow(count))

# set data
nimble_data <- list(y = count$count, year_cov = year_cov)

# initial values
inits_function <- function(){
  list(
    log_N_tau = 1,
    alpha = alpha_init,
    beta_trend = rnorm(nflocks, 0, 1),
    N = as.matrix(N_init),
    log_N = as.matrix(log(N_init + 0.001)),
    log_N_mu = matrix(rep(alpha_init, nyears), nflocks, nyears)
  )
}

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("trend_model", "nimble_data", "nimble_constants"))

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
  model <- nimbleModel(code = trend_model, constants = nimble_constants,  dat =  nimble_data, inits = nimble_inits)
  
  # initialize remaining parameters
  model$initializeInfo()
  model$calculate()
  
  # configure MCMC
  mcmc_Conf  <- configureMCMC(model)
  
  # build MCMC
  modelMCMC  <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  # run model
  CmodelMCMC$run(niter = 216000,
                 nburnin = 200000,
                 thin = 4
  )
  
  return(as.matrix(CmodelMCMC$mvSamples))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
save(out, file = "flock_trends_mod.RData")

# convert to MCMC list
samples <- list(chain1 = out[[1]], 
                chain2 = out[[2]], 
                chain3 = out[[3]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# traceplots
MCMCtrace(mcmc_list, type = "trace", filename = "flock_trends_mod_traceplots.pdf")

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]

# summary statistics
sum_stats <- MCMCsummary(mcmc_list)
write.csv(sum_stats, "flock_trends_mod_summary.csv")

