# compare the EE implementation in ergm to the EE implementation by Stival
# from Stivala's end we will implement his python polblogs example separately to this
# here we will just use ergm EE/SA/RM for the estimation of polblogs

################################################## 
# Run this whenever changes have been made to ergm
installr::uninstall.packages("ergm")
devtools::install_github("youhancheery/ergm")
library(ergm)
##################################################
library(network)
library(tidyr)
library(ggraph)

sd_estimate <- function(net_est) {
  # accepts a network variance-covariance matrix
  return(sqrt(diag(vcov(net_est, sources="estimation"))))
  
}

abs_delta_sd_estimate <- function(net_one, net_two) {
  stat_net_one <- sqrt(diag(vcov(net_one, sources="estimation")))
  stat_net_two <- sqrt(diag(vcov(net_two, sources="estimation")))
  abs_delta <- stat_net_one - stat_net_two
  return(abs_delta)
}

# prepare the polblogs dataset
polblogs <- read.csv("~/math5005/msresearch/data/polblogs_arclist.txt", header = F)
# first 1492 rows are just an index for number of nodes
polblogs <- polblogs[1493:nrow(polblogs),] %>% as.data.frame()
colnames(polblogs) <- "node1"
polblogs <- polblogs %>% separate(node1, c("node1", "node2")) # split on spaces
head(polblogs)
# convert to network
polblogs_net <- as.network(polblogs, matrix.type = "edgelist")
plot(polblogs_net, vertex.cex = 3)

# test stochastic approximation
sa_polblogs <- ergm(polblogs_net ~ edges + triangles, 
                    control = control.ergm(main.method = "Stochastic-Approximation",
                                           MCMC.return.stats = TRUE,
                                           MCMC.runtime.traceplot = TRUE))

# test robbins monro
rm_polblogs <- ergm(polblogs_net ~ edges + triangles,
                    control = control.ergm(main.method = "Robbins-Monro",
                                           MCMC.return.stats = TRUE,
                                           MCMC.runtime.traceplot = TRUE))

# test equilibrium expectation
ee_polblogs <- ergm(polblogs_net ~ edges + triangles, 
                    control = control.ergm(main.method = "EE"))

abs_delta_sd_estimate(sa_polblogs, ee_polblogs)

(sd_estimate(sa_polblogs))
(sd_estimate(ee_polblogs))
