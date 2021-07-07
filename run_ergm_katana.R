# test studies_comparison.R on kata
# some modifications applied

# measurement on the three datasets Pavel's contrastive divergence paper
# states are challenging for a network to estimate
m <- modules::use("~/math5005/msresearch/utils.R")
m$reinstall_ergm()

library(ergm)
library(ergm.count)
library(parallel)

set.seed(1234)

REPS <- 5

# load all the different datasets
# load("~/math5005/msresearch/data/supp_datasets.RData")
load("~/msresearch/data/supp_datasets.RData") # for katana



run_stoch_approx <- function(name, data, formula, control, estimate, reps=REPS) {
  # run the stochastic approximation algorithm across different starting points
  final <- list()
  init_method <- NULL
  runtime <- system.time(output <- ergm::ergm(formula=formula, 
                                              control=control,
                                              estimate=estimate,
                                              eval.loglik = FALSE))
  # coef <- stats::coef(output)
  coef <- output$coef
  iter <- output$iterations
  is_fail <- output$failure
  init <- output$coef.init
  final[paste0(name, ":",
               "SA",":",
               init_method)] <- list(list("init"=init,
                                          "timings"=runtime,
                                          "coef"=coef,
                                          "iterations"=iter,
                                          "is_fail"=is_fail)) 
  return(final)
}

out <- run_stoch_approx("ecoli2",
                        ecoli2, 
                        ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE),
                        control.ergm(main.method = "Stochastic-Approximation",
                                     force.main = TRUE,
                                     SA.niterations = 100,
                                     SA.nsubphases = 6,
                                     SA.burnin = 10000,
                                     SA.interval = 5000,
                                     SA.samplesize = 1024,
                                     MCMC.interval = 8),
                        estimate="MLE")















