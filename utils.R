# utils file for reusable functions

# Run this whenever changes have been made to ergm
reinstall_ergm <- function() {
  installr::uninstall.packages("ergm")
  devtools::install_github("youhancheery/ergm")
}

# estimate_parameters <- function {}

time_response <- function(data, formula) {

  # different methods we will be running through
  # TODO: understand why robbins-monro isn't running properly
  methods <- c("Stochastic-Approximation", "EE")
  
  final <- list()
  
  # initialisation methods are
  # 1. vector of zeros
  # 2, MPLE
  # 3. vector of zeros except edges density parameters (TODO)
  
  init_methods <- c("zeros", "MPLE")
  for(method in methods) {
    for(init_method in init_methods) {
      print(paste("Running estimation using", method, "starting at", init_method))
      if(init_method == "zeros") {
        # starting point initialised at all zeros
        # TODO need to add MCMC.samplesize to be high enough number
        control = ergm::control.ergm(main.method = method, # estimation method
                                     init=NULL, # always null so init.method passes values
                                     init.method=init_method, # 
                                     MCMC.runtime.traceplot = FALSE,
                                     MCMC.return.stats = TRUE,
                                     MCMLE.steplength = .25,
                                     MCMC.burnin = 100000, 
                                     MCMC.interval = 50000, 
                                     MCMC.samplesize = 2500,
                                     MCMLE.maxit = 100)
        runtime <- system.time(output <- ergm::ergm(formula=formula, 
                                                    control=control,
                                                    eval.loglik = FALSE))
        # coef <- stats::coef(output)
        coef <- output$coef
        iter <- output$iterations
        is_fail <- output$failure
        init <- output$coef.init
      } else {
        # using bernoulli reference we have initial starting point at MPLE
        control = ergm::control.ergm(main.method = method, # estimation method
                                     init=NULL, # initial values
                                     init.method=init_method, # 
                                     MCMC.runtime.traceplot = FALSE,
                                     MCMC.return.stats = TRUE,
                                     MCMLE.steplength = .25,
                                     MCMC.burnin = 100000, 
                                     MCMC.interval = 50000, 
                                     MCMC.samplesize = 2500,
                                     MCMLE.maxit = 100)
        runtime <- system.time(output <- ergm::ergm(formula=formula,
                                                    control=control,
                                                    reference=~Bernoulli,
                                                    eval.loglik = FALSE))
        # coef <- stats::coef(output)
        coef <- output$coef
        iter <- output$iterations
        is_fail <- output$failure
        init <- output$coef.init
      }
      
      final[paste0(method, ":", init_method)] <- list(list("init"=init,
                                                           "timings"=runtime,
                                                           "coef"=coef,
                                                           "iterations"=iter,
                                                           "is_fail"=is_fail)) 
    }
  }
  return(final)
}

mple_base <- function() {
  
}
