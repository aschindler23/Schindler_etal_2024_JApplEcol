
library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(MCMCvis)

### load data
# count data
count <- read.csv("count_data.csv")

# weather data
weather <- read.csv("weather_data.csv") %>% 
  arrange(year_id)

# % land cover data
pland <- read.csv("pland_data.csv")

# growing degree days data
winter_GDD <- read.csv("winter_GDD_data.csv")

### format data
# scaled count data (to estimate density dependence)
count_scaled <-  count %>%
  group_by(site_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, site_id, count) %>% 
  mutate(count = as.numeric(scale(count, center = T, scale = T))) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y") %>% 
  select(!site_id) %>% 
  as.matrix()

# % grass cover
grass <- pland %>% 
  select(site_id, year_id, scaled_grass) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_grass,
              names_prefix = "y") %>% 
  arrange(site_id) %>% 
  select(-site_id) %>% 
  as.matrix()

# % cereal cover
cereal <- pland %>% 
  select(site_id, year_id, scaled_cereal) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_cereal,
              names_prefix = "y") %>% 
  arrange(site_id) %>% 
  select(-site_id) %>% 
  as.matrix()

# % bog cover
bog <- pland %>% 
  select(site_id, year_id, scaled_peat_bog) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_peat_bog,
              names_prefix = "y") %>% 
  arrange(site_id) %>% 
  select(-site_id) %>% 
  as.matrix()

# growing degree days
winter_GDD <- winter_GDD %>% 
  select(site_id, year_id, scaled_winter_GDD) %>% 
  pivot_wider(names_from = year_id, values_from = scaled_winter_GDD,
              names_prefix = "y") %>% 
  arrange(site_id) %>% 
  select(-site_id) %>% 
  as.matrix()

### format initial values
# initial values for population size
N_init <- count %>%
  group_by(site_id) %>% 
  mutate(count = round(zoo::na.approx(count, na.rm = F))) %>% 
  mutate(count = zoo::na.locf(count, fromLast = T, na.rm = F)) %>% 
  ungroup() %>% 
  mutate(count = replace_na(count, 0)) %>% 
  select(year_id, site_id, count) %>% 
  pivot_wider(names_from = year_id, values_from = count, 
              names_prefix = "y") %>% 
  select(!site_id) %>% 
  as.matrix()
N1 <- N_init[,1]
N_init <- N_init + 0.001

# drop NA counts to speed up MCMC sampling
count <- count %>% 
  select(site_id, year_id, count) %>% 
  arrange(site_id, year_id) %>% 
  drop_na(count)

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

# number of observations in the count data
nobs <- nrow(count)

# number of environmental covariates
nenv <- 10

# number of classes
K <- 2
# K <- 3

### LCA classifications
Z <- read.csv(paste0("LCA_K", K, "_assignments.csv")) %>% 
  mutate(Z = as.numeric(str_remove(class, "K"))) %>% 
  select(site, Z) %>% 
  rename(site_id = site) %>% 
  arrange(site_id)


### state-space model to quantify environmental drivers of abundance trends
# specify model in BUGS language 
ss_env_model <- nimbleCode({
  
  #===========================
  # prior specification
  #===========================
  
  for(i in 1:nsites){ # loop over sites
    # initial population size
    N[i, 1] <- N1[i]
    log_N[i, 1] <- log(N[i, 1])
    
    # process variance for lambda
    sig_lambda[i] ~ dunif(0, 10)
    
    # site-specific density dependence
    beta_density[i] ~ dnorm(0, sd = 10)
  } # i
  
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
  for(i in 1:nsites){ # loop over sites
    for(t in 2:(max_years[i])){ # loop over years
      log_N[i, t] <- log_N[i, t - 1] + log_lambda[i, t]
      log_lambda[i, t] ~ dnorm(mu_lambda[i, t], sd = sig_lambda[i])
      mu_lambda[i, t] <- 
        alpha_lambda[Z[i]] + 
        beta_density[i] * count[i, t - 1] +
        beta_env[1, Z[i]] * spring_storm_days[t] + 
        beta_env[2, Z[i]] * spring_precip[t] + 
        beta_env[3, Z[i]] * pre_breed_freeze[t] +
        beta_env[4, Z[i]] * post_breed_precip[t] +
        beta_env[5, Z[i]] * autumn_storm_days[t] +
        beta_env[6, Z[i]] * autumn_freeze[t] +
        beta_env[7, Z[i]] * grass[i, t] + 
        beta_env[8, Z[i]] * cereal[i, t] + 
        beta_env[9, Z[i]] * bog[i, t] +
        beta_env[10, Z[i]] * winter_GDD[i, t - 1]
      N[i, t] <- exp(log_N[i, t])
    } # t
  } # i
  
  # observation model for count data
  for(n in 1:nobs){ # loop over count data observations
    y[n] ~ dpois(N[site_id[n], year_id[n]])
  } # n
})

# set constants
nimble_constants <- list(nsites = nsites,
                         max_years = max_years$max_year,
                         nenv = nenv,
                         K = K, 
                         Z = Z$Z,
                         site_id = count$site_id,
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
       mu_lambda = matrix(rnorm(nsites * nyears, 0, 0.1), nsites, nyears),
       beta_env = matrix(rnorm(nenv * K, 0, 1), nenv, K),
       beta_density = rnorm(nsites, 0, 1),
       alpha_lambda = rnorm(K, 0, 1),
       sig_lambda = runif(nsites, 0.1, 1),
       log_lambda = matrix(rnorm(nsites * nyears, 0, 0.1), nsites, nyears)
  )
}

inits <- inits_function()

# Select number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("ss_env_model", "nimble_data", "nimble_constants", "inits"))

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
  
  return(list(
    samples = samples
  ))
})
end_time <- Sys.time()
(run_time <- end_time - start_time)
stopCluster(cl)

# save output
file_heading <- paste0("K", K, "_ss_env")
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
