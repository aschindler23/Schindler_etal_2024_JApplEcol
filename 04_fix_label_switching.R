
library(tidyverse)
library(label.switching)
library(coda)
library(MCMCvis)


### define number of classes
K <- 2
# K <- 3

### load and format results
load(paste0("LCA_K", K, "_results.RData"))

# define number of chains
nc <- 3

# convert to MCMC list
samples <- list()
for(c in 1:nc){
  samples[[c]] <- as.data.frame(out[[c]][[1]])
}

### fix label switching
# combine all chains for parameters indexed by K
samples_mat <- as.matrix(bind_rows(samples))
Z <- samples_mat[, str_detect(colnames(samples_mat), "Z")]
w <- samples_mat[, str_detect(colnames(samples_mat), "w")]
class_mu_trend <- samples_mat[, str_detect(colnames(samples_mat), "class_mu_trend")]
class_sig_trend <- samples_mat[, str_detect(colnames(samples_mat), "class_sig_trend")]
class_mu_lat <- samples_mat[, str_detect(colnames(samples_mat), "class_mu_lat")]
class_sig_lat <- samples_mat[, str_detect(colnames(samples_mat), "class_sig_lat")]
class_mu_lon <- samples_mat[, str_detect(colnames(samples_mat), "class_mu_lon")]
class_sig_lon <- samples_mat[, str_detect(colnames(samples_mat), "class_sig_lon")]


# run function to identify label switching
nMCMC <- nrow(samples_mat)
run <- ecr.iterative.1(z = Z, K = K, maxiter = 100000)
neworder <- run$permutations

# reorder classes to correct labels
w_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_mu_trend_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_sig_trend_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_mu_lat_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_sig_lat_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_mu_lon_reorder <- matrix(0, nrow = nMCMC, ncol = K)
class_sig_lon_reorder <- matrix(0, nrow = nMCMC, ncol = K)
Z_reorder <- matrix(0, nrow = nMCMC, ncol = 59)

for(i in 1:nMCMC) {
  w_reorder[i,] <- w[i,][neworder[i,]]
  class_mu_trend_reorder[i,] <- class_mu_trend[i,][neworder[i,]]
  class_sig_trend_reorder[i,] <- class_sig_trend[i,][neworder[i,]]
  class_mu_lat_reorder[i,] <- class_mu_lat[i,][neworder[i,]]
  class_sig_lat_reorder[i,] <- class_sig_lat[i,][neworder[i,]]
  class_mu_lon_reorder[i,] <- class_mu_lon[i,][neworder[i,]]
  class_sig_lon_reorder[i,] <- class_sig_lon[i,][neworder[i,]]
  for(k in 1:K){
    Z_reorder[i,which(Z[i,] == k)] <- neworder[i, k]  
  }
}

# replace results with corrected class assignments
samples <- list()
for(c in 1:nc){
  samples[[c]] <- out[[c]][[1]]
}

nMCMC_1chain <- nrow(samples_mat) / nc
first <- c(1, nMCMC_1chain + 1, nMCMC_1chain * 2 + 1)
last <- c(nMCMC_1chain, nMCMC_1chain * 2, nMCMC_1chain * 3)
for(c in 1:nc){
  samples[[c]][,str_detect(colnames(samples_mat), "Z")] <- Z_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_mu_trend")] <- class_mu_trend_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_sig_trend")] <- class_sig_trend_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_mu_lat")] <- class_mu_lat_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_sig_lat")] <- class_sig_lat_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_mu_lon")] <- class_mu_lon_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "class_sig_lon")] <- class_sig_lon_reorder[first[c]:last[c],]
  samples[[c]][,str_detect(colnames(samples_mat), "w")] <- w_reorder[first[c]:last[c],]
}

mcmc_list <- as.mcmc.list(lapply(samples, mcmc))

# save results
file_heading <- paste0("LCA_K", K)
save(mcmc_list, file = paste0(file_heading, "_results_corrected.RData"))

# corrected traceplots
MCMCtrace(mcmc_list, type = "trace", iter = 4000,
          filename = paste0(file_heading, "_traceplots_corrected.pdf"))

# remove NA parameters
for(c in 1:nc) {
  mcmc_list[[c]] <- mcmc_list[[c]][ , which(is.na(colSums(mcmc_list[[c]])) == F)]
}

# assess convergence
rhat <- gelman.diag(mcmc_list, multivariate = F)
rhat <- unlist(rhat$psrf)
rhat[which(rhat[,1] > 1.1), ]
rhat <- as.data.frame(rhat)
colnames(rhat) <- c("Rhat", "Rhat_UCI")

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

# calculate and save class assignments
calc_prop_Z <- function(x, r, K, n){
  x <- x[, r]
  prop_Z <- as.data.frame(matrix(NA, 1, K))
  for(k in 1:K){
    prop_Z[k] <- length(x[x == k]) / n
    colnames(prop_Z)[k] <- paste0("K", k)
  }
  prop_Z$site <- r
  prop_Z
}

Z_prop <- lapply(1:dim(Z)[2], calc_prop_Z, x = Z_reorder, K = K, n = nMCMC)
Z_prop <- bind_rows(Z_prop) 
Z_prop$class <- colnames(Z_prop)[max.col(Z_prop[1:K])]
table(Z_prop$class)
write.csv(Z_prop, paste0(file_heading, "_assignments.csv"), row.names = F)
