library(tidyverse)
library(NbClust)

### load trend results
# load results and select median trend estimates
trend_sum <- read.csv("site_trends_summary.csv") %>% 
  select(X, X50.) %>% 
  filter(str_detect(X, "beta_trend")) %>% 
  rename(param = X, median_trend = X50.)

# convert to a matrix
trend_median <- as.matrix(select(trend_sum, median_trend))

### compute 21 different indices
# create vector of index names
indices <- c("ball", "beale", "ch", "cindex", "db", 
             "dindex", "duda", "dunn", "frey", "gamma", 
             "gplus", "hartigan", "hubert", "kl", "mcclain", 
             "ptbiserial", "ratkowsky", "sdbw", "sdindex", "silhouette", 
             "tau")

# create data frame
n_clust_results <- data.frame(method = indices, n_clust = NA)


# calculate indices and fill data frame with results
seq <- 1:length(indices)
seq <- seq[-c(6, 13)] # remove two graphical methods
for(i in 1:seq) {
  n_clust_results[i, 2] <- as.numeric(NbClust(data = trend_median, 
                                              diss = NULL, distance = "euclidean", 
                                              min.nc = 2, max.nc = 10, 
                                              method = "kmeans", 
                                              index = indices[i])$Best.nc[1])
}

# plot graphical methods and interpret
NbClust(data = trend_median, 
        diss = NULL, 
        distance = "euclidean", 
        min.nc = 2, 
        max.nc = 10, 
        method = "kmeans", 
        index = "dindex")
n_clust_results[6, 2] <- 3

NbClust(data = trend_median, 
        diss = NULL, 
        distance = "euclidean", 
        min.nc = 2, 
        max.nc = 10, 
        method = "kmeans", 
        index = "hubert")
n_clust_results[13, 2] <- 4


# generate table
table(n_clust_results$n_clust)
write.csv(n_clust_results, file = "Cluster Validity indices results.csv", row.names = F)
