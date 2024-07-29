
library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load and process data
# count data
count <- read.csv("count_data.csv")

# initial values for population size
N_init <- count %>%
  group_by(site_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, site_id, count) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y")
N_init <- N_init[,-1] + 0.001

# initial values for mean population size (log scale)
alpha_init <- count %>% 
  group_by(site_id) %>% 
  summarise(mean_count = mean(count, na.rm = T)) %>%
  ungroup() %>% 
  mutate(log_mean_count = log(mean_count)) %>% 
  arrange(site_id)

# drop NA counts to speed up MCMC sampling
count <- count %>%
  drop_na() %>%
  arrange(site_id, year_id)

### indexing
# maximum years
max_years <- count %>% 
  select(site_id, year_id) %>% 
  group_by(site_id) %>% 
  summarise(max_year = max(year_id)) %>% 
  ungroup() %>% 
  arrange(site_id)

# number of sites
nsites <- length(unique(count$site_id))

# number of years
nyears <- length(unique(count$year_id))

year_cov <- c(1:nyears) - nyears/2 - .5

### trend model
# specify model in BUGS language 
trend_model <- nimbleCode({
  #===========================
  # priors and linear models
  #===========================
  
  for(i in 1:nsites){ # loop over sites
    # site-year random intercept
    alpha[i] ~ dnorm(0, sd = 10)
    
    # site-specific trend in abundance
    beta_trend[i] ~ dnorm(0, sd = 10)
  } # i
  
  # process variance
  log_N_tau ~ dgamma(0.01, 0.01)
  
  #===========================
  # likelihood specification
  #===========================
  
  # process model for true abundance
  for(i in 1:nsites){ # loop over sites
    for(t in 1:(max_years[i])){ # loop over years
      log_N[i, t] ~ dnorm(log_N_mu[i, t], log_N_tau)
      log_N_mu[i, t] <- alpha[i] + beta_trend[i] * year_cov[t]
      N[i, t] <- exp(log_N[i, t])
    } # t
  } # i
  
  # observation model for count data
  for(n in 1:ncount){ # loop over count data
    y[n] ~ dpois(N[site[n], year[n]]) 
  } # n
})

# set constants
nimble_constants <- list(max_years = max_years$max_year,
                         nsites = nsites,
                         site = count$site_id,
                         year = count$year_id,
                         ncount = nrow(count))

# set data
nimble_data <- list(y = count$count, year_cov = year_cov)

# initial values
inits_function <- function(){
  list(
    log_N_tau = 1,
    alpha = alpha_init$log_mean_count,
    beta_trend = rnorm(nsites, 0, 1),
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

clusterExport(cl, c("trend_model", "nimble_data", "nimble_constants", "inits"))

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
  
  return(list(
    samples = samples
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
file_heading <- "site_trends"
save(out, file = paste0(file_heading, "_results.RData"))

# convert to MCMC list
samples    <- list(chain1 = out[[1]][[1]], 
                   chain2 = out[[2]][[1]], 
                   chain3 = out[[3]][[1]])

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

for(c in 1:nc) {
  mcmc_list[[c]] <- mcmc_list[[c]][ , which(is.na(colSums(mcmc_list[[c]])) == F)]
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

# calculate summary statistics
sum_stats <- MCMCsummary(mcmc_list, func = function(x) prob_above_below_0(x), func_name = "f")

# save summary statistics
write.csv(sum_stats, paste0(file_heading, "_summary.csv"))

