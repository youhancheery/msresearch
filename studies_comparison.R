# measurement on the three datasets Pavel's contrastive divergence paper
# states are challenging for a network to estimate

# m$reinstall_ergm()
library(ergm)
library(ergm.count)
library(parallel)

# set.seed(1234)

# Simulation ranges
K.LEVELS <- c(1,2,4,16,128)
M.LEVELS <- c(1,2,4,8)
KM.MAX <- 256
STEPLEN.MARGINS <- c(0.05,0.5,1)
RM.CONFS <- list(list(a0=0.1, c=0.5), list(a0=0.5, c=0.5), list(a0=1, c=0))
REPS <- 5
MPLE <- TRUE
INTERCEPT <- TRUE



# load all the different datasets
load("~/math5005/msresearch/data/supp_datasets.RData")
m <- modules::use("~/math5005/msresearch/utils.R") 

# ecoli
# testing that time_response produces something
output_ecoli2 <- m$time_response(data=ecoli2,
                                 ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE),
                                 seed=1234)

output_ecoli2$`Stochastic-Approximation:zeros`$coef
output_ecoli2$`Stochastic-Approximation:NULL`$coef

ouput_ecoli2_self <- m$time_response(data=ecoli2,
                                     ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE)
                                     + nodemix("self", base = 1),
                                     seed=1234)

# testing that parallel run works
test <- m$run_study(data = ecoli2,
                    ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE))




# time response for each model (3) across each dataset (3) for a total of 9 fits
# kapferer
study_kapferer <- m$time_response(data=kapferer,
                                  formula=kapferer~edges+gwdegree(0.25,fixed=TRUE)+gwesp(0.25, fixed=TRUE)+gwdsp(0.25,fixed=TRUE))
# zach
# TODO: not sure why this isn't working? 
# nonzero deprecated in current version of ergm?
zach_study <- m$time_response(data = zach,
                              formula = zach ~ nonzero+sum+nodefactor("role",base=2)+transitiveweights)

