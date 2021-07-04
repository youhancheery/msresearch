# measurement on the three datasets Pavel's contrastive divergence paper
# states are challenging for a network to estimate

m <- modules::use("~/Documents/math5005/msresearch/utils.R") 

# m$reinstall_ergm()
library(ergm)
library(ergm.count)

set.seed(1234)

# load all the different datasets
load("~/Documents/math5005/msresearch/data/supp_datasets.RData")

# create a function that reads any single dataset
# and builds a model off the back of it
# need to time the time taken to reach the result

# time response for each model (3) across each dataset (3) for a total of 9 fits
# kapferer
kapferer_study <- m$time_response(data=kapferer,
                                  formula=kapferer~edges+gwdegree(0.25,fixed=TRUE)+gwesp(0.25, fixed=TRUE)+gwdsp(0.25,fixed=TRUE))
# ecoli
ecoli2_study <- m$time_response(data=ecoli2,
                                ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE))
ecoli2_self_study <- m$time_response(data=ecoli2,
                                     ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1))
# zach
# TODO: not sure why this isn't working? 
# nonzero deprecated in current version of ergm?
zach_study <- m$time_response(data = zach,
                              formula = zach ~ nonzero+sum+nodefactor("role",base=2)+transitiveweights)



run.study("kap.study", kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE))
run.study("zach.binom.study", zach ~ nonzero+sum+nodefactor("role",base=2)+transitiveweights, reference = ~Binomial(8), response = "contexts", ctrl.cd=control.ergm(CD.maxit=60,CD.trustregion=100000,MCMLE.trustregion=100000))
run.study("zach.cmp.study", zach ~ nonzero+CMP+sum+nodefactor("role",base=2)+transitiveweights, reference  = ~Poisson, response = "contexts", ctrl.cd=control.ergm(CD.maxit=60,CD.trustregion=100000,MCMLE.trustregion=100000))
run.study("kap2.study", kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE))
run.study("ecoli2.study", ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE))
run.study("ecoli2.self.study", ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1))
run.study("laz2.c.study", laz ~ gwesp()+nodecov("seniority")+nodefactor("practice")+nodematch("practice")+nodematch("gender")+nodematch("office"), constraints = ~edges)
run.study("laz2.study", laz ~ edges+gwesp()+nodecov("seniority")+nodefactor("practice")+nodematch("practice")+nodematch("gender")+nodematch("office"))




