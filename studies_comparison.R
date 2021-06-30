# measurement on the three datasets Pavel's contrastive divergence paper
# states are challenging for a network to estimate

m <- modules::use("math5005/msresearch/utils.R") 

m$reinstall_ergm()
library(ergm)

set.seed(1234)

# load all the different datasets
load("math5005/msresearch/data/supp_datasets.RData")

# create a function that reads any single dataset
# and builds a model off the back of it
# need to time the time taken to reach the result

t <- system.time(test_kap <- ergm(kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE),
                                  control = control.ergm(main.method = "Stochastic-Approximation",
                                                         init = NULL,
                                                         init.method = "zeros",
                                                         MCMC.runtime.traceplot = TRUE,
                                                         MCMC.return.stats = TRUE)))
