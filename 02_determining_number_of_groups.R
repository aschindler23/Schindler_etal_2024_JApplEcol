
setwd("C:/Users/gth492/OneDrive - University of Saskatchewan/Flock count analysis/revised analysis 070424")

library(tidyverse)
library(factoextra)
library(NbClust)

# load trend results
trend_sum <- read.csv("flock_trends_linear_070424_summary.csv") %>% 
  select(X, X50.) %>% 
  filter(str_detect(X, "beta_trend")) %>% 
  rename(param = X, median_trend = X50.)
trend_median <- as.matrix(select(trend_sum, median_trend))

# composite of 23 different indices
indices <- c("ball", "beale", "ch", "cindex", "db", "dindex", "duda", "dunn", "frey", "gamma", 
             "gap", "gplus", "hartigan", "hubert", "kl", "mcclain", "pseudot2", "ptbiserial", "ratkowsky", "sdbw",
             "sdindex", "silhouette", "tau")
n.clust.results <- data.frame(method = indices, n.clust = NA)
seq <- 1:length(indices)
seq <- seq[-c(6, 11, 14, 17)] # remove two graphical methods
for (i in seq) {
  n.clust.results[i,2] <- as.numeric(NbClust(data = trend_median, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans", index = indices[i])$Best.nc[1])
}

# plot graphical methods and interpret
NbClust(data = trend_median, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans", index = "hubert")
n.clust.results[14, 2] <- 4
NbClust(data = trend_median, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans", index = "dindex")
n.clust.results[6, 2] <- 3
NbClust(data = trend_median, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans", index = "gap")

# generate table
table(n.clust.results$n.clust)
write.csv(n.clust.results, file = "Cluster Validity indices results.csv", row.names = F)
write.csv(table(n.clust.results$n.clust), file = "Cluster Validity indices results2.csv", row.names = F)
