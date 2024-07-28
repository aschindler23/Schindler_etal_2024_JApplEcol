
library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)
setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Nimble resources")
source("nimble_restart_functions.R")

setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/revised analysis 070424")

### load and process data
count <- read.csv("data/count_data.csv")

# initial values for population size
N_init <- count %>%
  group_by(flock_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, flock_id, count) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y")
N_init <- N_init[,-1] + 0.001

# initial values for mean population size (log scale)
alpha_init <- count %>% 
  group_by(flock_id) %>% 
  summarise(mean_count = mean(count, na.rm = T)) %>%
  ungroup() %>% 
  mutate(log_mean_count = log(mean_count)) %>% 
  arrange(flock_id)

# drop NA counts to speed up MCMC sampling
count <- count %>%
  drop_na() %>%
  arrange(flock_id, year_id)

### indexing
# maximum years
max_years <- count %>% 
  select(flock_id, year_id) %>% 
  group_by(flock_id) %>% 
  summarise(max_year = max(year_id)) %>% 
  ungroup() %>% 
  arrange(flock_id)

# number of flocks
nflocks <- length(unique(count$flock_id))

# number of years
nyears <- length(unique(count$year_id))

year_cov <- c(1:nyears) - nyears/2 - .5

### trend model
# specify model in BUGS language 
trend_model <- nimbleCode({
  #===========================
  # priors and linear models
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
  for(n in 1:ncount){ # loop over count data
    y[n] ~ dpois(N[flock[n], year[n]]) 
  } # n
})

# set constants
nimble_constants <- list(max_years = max_years$max_year,
                         nflocks = nflocks,
                         flock = count$flock_id,
                         year = count$year_id,
                         ncount = nrow(count))

# set data
nimble_data <- list(y = count$count, year_cov = year_cov)

# initial values
inits_function <- function(){
  list(
    log_N_tau = 1,
    alpha = alpha_init$log_mean_count,
    beta_trend = rnorm(nflocks, 0, 1),
    N = as.matrix(round(N_init)),
    log_N = as.matrix(log(N_init)),
    log_N_mu = as.matrix(log(N_init))
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("trend_model", "nimble_data", "nimble_constants", "inits", 
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
  model <- nimbleModel(code = trend_model, constants = nimble_constants,  dat =  nimble_data, inits = nimble_inits)
  
  # initialize remaining parameters
  # model$initializeInfo()
  # model$check()
  # model$calculate()
  
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
file_heading <- "flock_trends_linear_070424"
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples    <- list(chain1 = out[[1]][[1]], 
                   chain2 = out[[2]][[1]], 
                   chain3 = out[[3]][[1]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

for (i in 1:nc) {
  mcmc_list[[i]] <- mcmc_list[[i]][ , which(is.na(colSums(mcmc_list[[i]])) == F)]
}

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]
nonconverged <- rownames(rhat[which(rhat[,1] > 1.1), ])

# traceplots
MCMCtrace(mcmc_list, type = "trace", iter = 4000,
          filename = paste0(file_heading, "_traceplots.pdf"))

# function to calculate percent of posterior above or below 0
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

# caculate summary statistics
sum_stats <- MCMCsummary(mcmc_list, func = function(x) prob_above_below_0(x), func_name = "f")

# save summary statistics
write.csv(sum_stats, paste0(file_heading, "_summary.csv"))

