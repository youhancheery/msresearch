# measurement on the three datasets Pavel's contrastive divergence paper
# states are challenging for a network to estimate

# m$reinstall_ergm()
library(ergm)
library(parallel)

ESTIMATION_METHODS <- c("Stochastic-Approximation", "Robbins-Monro", "EE", "MCMLE")
STARTING_ESTIMATES <- c("zeros", "MPLE", "zeros_and_edges")
N_CORES <- 16
set.seed(1234)

# load all the different datasets
load("~/math5005/msresearch/data/supp_datasets.RData")
# m <- modules::use("~/math5005/msresearch/utils.R") 

# # ecoli
# # testing that time_response produces something
# output_ecoli2 <- m$time_response(data=ecoli2,
#                                  ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE),
#                                  seed=1234)
# 
# output_ecoli2$`Stochastic-Approximation:zeros`$coef
# output_ecoli2$`Stochastic-Approximation:NULL`$coef
# 
# ouput_ecoli2_self <- m$time_response(data=ecoli2,
#                                      ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE)
#                                      + nodemix("self", base = 1),
#                                      seed=1234)
# 
# # testing that parallel run works
# test <- m$run_study(data = ecoli2,
#                     ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE))
# 
# 
# # time response for each model (3) across each dataset (3) for a total of 9 fits
# # kapferer
# study_kapferer <- m$time_response(data=kapferer,
#                                   formula=kapferer~edges+gwdegree(0.25,fixed=TRUE)+gwesp(0.25, fixed=TRUE)+gwdsp(0.25,fixed=TRUE))
# # zach
# # TODO: not sure why this isn't working? 
# # nonzero deprecated in current version of ergm?
# zach_study <- m$time_response(data = zach,
#                               formula = zach ~ nonzero+sum+nodefactor("role",base=2)+transitiveweights)


# put together the output as a list of lists
.construct_output <- function(name, main_method, init_method, ergm_output, runtime) {
  
  final <- list()
  # extract relevant metrics
  # keep the network itself to get traceplots
  coef <- coef(ergm_output) # TODO change to coef() when running on katana
  iter <- ergm_output$iterations
  is_fail <- ergm_output$failure
  init <- ergm_output$coef.init
  final[paste0(name, ":",
               main_method,":",
               init_method)] <- list(list("init"=init,
                                          "timings"=runtime,
                                          "coef"=coef,
                                          "iterations"=iter,
                                          "is_fail"=is_fail))
  print(runtime)
  print(coef)
  return(final)
}

# run ergm across all starting points, depending on the method specified by the user
# TODO add more control arguments to the function arguments as we test them
run_stoch_approx <- function(name, starting_est, data, formula, n_term_formula,
                             SA.nsubphases = 10,
                             main_method="Stochastic-Approximation",
                             estimate="MLE") {
  
  if(!(main_method %in% ESTIMATION_METHODS)) {stop("Main method not accurately input - see ergm docs for valid starting point")}
  if(!(starting_est %in% STARTING_ESTIMATES)) {stop(paste("Starting point not in", STARTING_ESTIMATES))}
  
  if (starting_est == "zeros") {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      SA.niterations = 2000,
      SA.nsubphases = SA.nsubphases,
      SA.burnin = 10000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE,
      MCMC.interval = 20)
    # pretty sure we don't need to include reference for zeros starting point?
    runtime <- system.time(ergm_output <- ergm(formula=formula,
                                               control=control,
                                               estimate=estimate,
                                               eval.loglik = FALSE))
    print(ergm_output$coef)
  } else if(starting_est == "zeros_and_edges") {
  
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # need to run ergm once for the edges which we'll use as our starting point
    first_pass <- ergm(formula = data ~ edges,
                       control = control.ergm(main.method = "MCMLE"),
                       estimate = "MLE",
                       eval.loglik = FALSE)
    # NOTE we assume that the first argument is "edges"
    edge <- coef(first_pass)
    theta0 <- c(edge[[1]], rep(0, n_term_formula - 1)) # everything 0 except edges parameter
    # build control, input NOT null
    control <- control.ergm(
      main.method = main_method,
      init = theta0, # don't need init.method since we supply init
      force.main = TRUE,
      SA.niterations = 2000,
      SA.nsubphases = SA.nsubphases,
      SA.burnin = 10000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE)
    # pretty sure we don't need to include reference for zeros starting point?
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               eval.loglik = FALSE))
  } else {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # MPLE starting point
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      SA.niterations = 20,
      SA.nsubphases = SA.nsubphases,
      SA.burnin = 1000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE)
    # need to provide reference =~Bernoulli for MPLE starting point
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               reference=~Bernoulli,
                                               eval.loglik = FALSE))
  }
  final <- .construct_output(name, main_method, starting_est, ergm_output, runtime)
  save(final, ergm_output, file=paste(paste0("outputs/", name), main_method, starting_est, "RData", sep="."))
}

run_mcmle <- function(name, starting_est, data, formula, n_term_formula,
                      MCMLE.maxit = 1000, # control parameters
                      main_method="MCMLE",
                      estimate="MLE") {
  # similar to run_stoch_approx except parameters in control.ergm are mcmle
  
  if(!(main_method %in% ESTIMATION_METHODS)) {stop("Main method not accurately input - see ergm docs for valid starting point")}
  if(!(starting_est %in% STARTING_ESTIMATES)) {stop(paste("Starting point not in", STARTING_ESTIMATES))}
  
  if (starting_est == "zeros") {

    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # i dont believe for zeros we need to add reference=~Bernoulli
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      MCMLE.termination = "Hotelling",
      MCMLE.density.guard = exp(7),
      MCMLE.maxit = MCMLE.maxit,
      MCMLE.metric = "naive", # start with naive but will try "lognormal" later too
      MCMC.samplesize=400,
      MCMC.runtime.traceplot = TRUE,
      MCMC.interval = 10000)
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               eval.loglik = FALSE))
  } else if(starting_est == "zeros_and_edges") {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # need to run ergm once for the edges which we'll use as our starting point
    first_pass <- ergm(formula = data ~ edges,
                       control = control.ergm(main.method = "MCMLE"),
                       estimate = "MLE",
                       eval.loglik = FALSE)
    # NOTE we assume that the first argument is "edges"
    edge <- coef(first_pass)
    theta0 <- c(edge[[1]], rep(0, n_term_formula - 1)) # everything 0 except edges parameter
    # build control, input NOT null
    control <- control.ergm(
      main.method = main_method,
      init = theta0, # don't need init.method since we supply init
      force.main = TRUE,
      MCMLE.termination = "Hotelling",
      MCMLE.density.guard = exp(7),
      MCMLE.maxit = MCMLE.maxit,
      MCMLE.metric = "Median.Likelihood",
      MCMC.interval = 10000,
      MCMC.samplesize = 400,
      MCMC.runtime.traceplot = TRUE)
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               eval.loglik = FALSE))
  } else {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
   # MPLE starting point
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      MCMLE.termination = "Hotelling",
      MCMLE.density.guard = exp(7),
      MCMLE.metric = "Median.Likelihood",
      MCMLE.maxit = MCMLE.maxit,
      MCMC.runtime.traceplot = TRUE,
      MCMC.samplesize = 400,
      MCMC.interval = 10000)
    # need to provide reference =~Bernoulli for MPLE starting point
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               reference=~Bernoulli,
                                               eval.loglik = FALSE))
  }
  final <- .construct_output(name, main_method, starting_est, ergm_output, runtime)
  save(final, ergm_output, file=paste(paste0("outputs/", name), main_method, starting_est, "RData", sep="."))
}

run_ee <- function(name, starting_est, data, formula, n_term_formula,
                   EE.nsubphases = 10,
                   main_method="EE",
                   estimate="MLE") {
  # similar to run_stoch_approx except parameters in control.ergm are ee
  # note that EE is built on top of SA, so the parameters should be the same
  # with the only replacement being obviously SA -> EE
  if(!(main_method %in% ESTIMATION_METHODS)) {stop("Main method not accurately input - see ergm docs for valid starting point")}
  if(!(starting_est %in% STARTING_ESTIMATES)) {stop(paste("Starting point not in", STARTING_ESTIMATES))}
  print(EE.nsubphases) 
  if (starting_est == "zeros") {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      SA.niterations = 2000,
      SA.nsubphases = EE.nsubphases,
      SA.burnin = 10000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE,
      MCMC.interval = 20)
    # pretty sure we don't need to include reference for zeros starting point?
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               eval.loglik = FALSE))
  } else if(starting_est == "zeros_and_edges") {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # need to run ergm once for the edges which we'll use as our starting point
    first_pass <- ergm(formula = data ~ edges,
                       control = control.ergm(main.method = "MCMLE"),
                       estimate = "MLE",
                       eval.loglik = FALSE)
    # NOTE we assume that the first argument is "edges"
    edge <- coef(first_pass)
    theta0 <- c(edge[[1]], rep(0, n_term_formula - 1)) # everything 0 except edges parameter
    # build control, input NOT null
    control <- control.ergm(
      main.method = main_method,
      init = theta0, # don't need init.method since we supply init
      force.main = TRUE,
      EE.niterations = 2000,
      SA.nsubphases = EE.nsubphases,
      SA.burnin = 10000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE)
    # pretty sure we don't need to include reference for zeros starting point?
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               eval.loglik = FALSE))
  } else {
    print(paste("running", main_method, "with", starting_est, "as the starting point"))
    # MPLE starting point
    control <- control.ergm(
      main.method = main_method,
      init = NULL, #otherwise actual theta0. must now input input.method argument
      init.method = starting_est,
      force.main = TRUE,
      SA.niterations = 2000,
      SA.nsubphases = EE.nsubphases,
      SA.burnin = 10000,
      SA.interval = 500,
      SA.samplesize = 1024,
      MCMC.runtime.traceplot = TRUE,
      MCMC.interval = 20)
    # need to provide reference =~Bernoulli for MPLE starting point
    runtime <- system.time(ergm_output <- ergm(formula=formula, 
                                               control=control,
                                               estimate=estimate,
                                               reference=~Bernoulli,
                                               eval.loglik = FALSE))
  }
  final <- .construct_output(name, main_method, starting_est, ergm_output, runtime)
  print(final)
  # save results on local
  save(final, ergm_output, file=paste(paste0("outputs/", name), main_method, starting_est, "RData", sep="."))
}

SA.NSUBPHASES <- 8
EE.NSUBPHASES <- 8
MCMLE.MAXIT   <- 60

###### ECOLI BASE 
###### STOCHASTIC-APPROXIMATION
#run_stoch_approx("ecoli2", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, SA.nsubphases = SA.NSUBPHASES)
#run_stoch_approx("ecoli2", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, SA.nsubphases = SA.NSUBPHASES)
#run_stoch_approx("ecoli2", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, SA.nsubphases = SA.NSUBPHASES)
###### MCMLE
#run_mcmle("ecoli2", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("ecoli2", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("ecoli2", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, MCMLE.maxit = MCMLE.MAXIT) 
#run_mcmle("ecoli2", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, MCMLE.maxit = MCMLE.MAXIT)
###### EQUILIBRIUM EXPECTATION
#run_ee("ecoli2", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, EE.nsubphases = EE.NSUBPHASES)
#run_ee("ecoli2", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, EE.nsubphases = EE.NSUBPHASES)
#run_ee("ecoli2", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE), 6, EE.nsubphases = EE.NSUBPHASES)

###### ECOLI WITH SELF LOOP
###### STOCHASTIC-APPROXIMATION
#run_stoch_approx("ecoli2.self.study", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, SA.nsubphases = SA.NSUBPHASES)
#run_stoch_approx("ecoli2.self.study", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, SA.nsubphases = SA.NSUBPHASES)
#run_stoch_approx("ecoli2.self.study", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, SA.nsubphases = SA.NSUBPHASES)
###### MCMLE
#run_mcmle("ecoli2.self.study", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("ecoli2.self.study", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("ecoli2.self.study", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, MCMLE.maxit = MCMLE.MAXIT)
###### EQUILIBRIUM EXPECTATION
#run_ee("ecoli2.self.study", "zeros", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, EE.nsubphases = EE.NSUBPHASES)
#run_ee("ecoli2.self.study", "MPLE", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, EE.nsubphases = EE.NSUBPHASES)
#run_ee("ecoli2.self.study", "zeros_and_edges", ecoli2, ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1), 8, EE.nsubphases = EE.NSUBPHASES)

###### KAPFERER BASE
###### STOCHASTIC-APPROXIMATION
#SA.NSUBPHASES = 7
run_stoch_approx("kapferer", "zeros", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, SA.nsubphases = SA.NSUBPHASES)
run_stoch_approx("kapferer", "MPLE", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, SA.nsubphases = SA.NSUBPHASES)
run_stoch_approx("kapferer", "zeros_and_edges", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, SA.nsubphases = SA.NSUBPHASES)
###### MCMLE
#run_mcmle("kapferer", "zeros", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("kapferer", "MPLE", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("kapferer", "zeros_and_edges", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, MCMLE.maxit = MCMLE.MAXIT)
###### EQUILIBRIUM EXPECTATION
run_ee("kapferer", "zeros", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, EE.nsubphases = EE.NSUBPHASES)
run_ee("kapferer", "MPLE", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, EE.nsubphases = EE.NSUBPHASES)
run_ee("kapferer", "zeros_and_edges", kapferer, kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 3, EE.nsubphases = EE.NSUBPHASES)

###### KAPFERER2
###### STOCHASTIC-APPROXIMATION
run_stoch_approx("kapferer2", "zeros", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, SA.nsubphases = SA.NSUBPHASES)
run_stoch_approx("kapferer2", "MPLE", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, SA.nsubphases = SA.NSUBPHASES)
run_stoch_approx("kapferer2", "zeros_and_edges", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, SA.nsubphases = SA.NSUBPHASES)
###### MCMLE
#run_mcmle("kapferer2", "zeros", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("kapferer2", "MPLE", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, MCMLE.maxit = MCMLE.MAXIT)
#run_mcmle("kapferer2", "zeros_and_edges", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, MCMLE.maxit = MCMLE.MAXIT)
###### EQUILIBRIUM EXPECTATION
run_ee("kapferer2", "zeros", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, EE.nsubphases = EE.NSUBPHASES)
run_ee("kapferer2", "MPLE", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, EE.nsubphases = EE.NSUBPHASES)
run_ee("kapferer2", "zeros_and_edges", kapferer, kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE), 4, EE.nsubphases = EE.NSUBPHASES)

