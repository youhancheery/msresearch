# utils file for reusable functions

# Run this whenever changes have been made to ergm
reinstall_ergm <- function() {
  installr::uninstall.packages("ergm")
  devtools::install_github("youhancheery/ergm")
}

time.ergm <- function(...,control,estimate,seed=0){
  o <- try({
    set.seed(seed) # For reproducibility.
    t <- system.time(out <- ergm(...,control=control,estimate=estimate,eval.loglik=FALSE))
    coef <- coef(out)
    list(timings=t, coef=coef, iterations=if(out$estimate!="MPLE") out$iterations else 1)
  })
  
  if(inherits(o,"try-error")) o <- list(error=TRUE,error.msg=o,error.when="time.ergm") else o$error=FALSE
  
  c(list(est=estimate,ctrl=control,seed=seed),o)
}

time_response <- function(formula, method, reference=~Bernouli) {
  
  control <- control.ergm()
  
  runtime <- system.time(output <- ergm(control=control, eval.loglik = FALSE))
  
}

run.study <- function(name, formula, response=NULL, reference=~Bernoulli, constraints=~., 
                      ctrl.cd=control.ergm(CD.maxit=60,CD.trustregion=100000, CD.Hummel.miss.sample=10),
                      ctrl.simulate=control.simulate.formula(), nsim=1024, 
                      k.levels=K.LEVELS, m.levels=M.LEVELS, km.max=KM.MAX, 
                      cd.metrics=CD.METRICS, steplen.margins=STEPLEN.MARGINS, 
                      reps=REPS, mple=MPLE, intercept=INTERCEPT, rm.confs=RM.CONFS){
  
  ctrl.cd$MCMC.samplesize <- 1024
  ctrl.cd$parallel <- 0
  
  lapply(seq_len(reps), function(rep){
    # Start the intercept algorithm.
    if(intercept){
      limitmyps(PSLIMIT)
      mcparallel({study.fit(name, formula, response=response, reference=reference, constraints=constraints, estimate="intercept", ctrl.cd=ctrl.cd, ctrl.simulate=ctrl.simulate, nsim=nsim, seed=rep)},detached=TRUE)
    }
    # Start the MPLE.
    if(mple){
      limitmyps(PSLIMIT)
      mcparallel({study.fit(name, formula, response=response, reference=reference, constraints=constraints, estimate="MPLE", ctrl.cd=ctrl.cd, ctrl.simulate=ctrl.simulate, nsim=nsim, seed=rep)},detached=TRUE)
    }
    
    # Start the CD RM with the specified parameters.
    for(rm.conf in rm.confs){
      for(k in k.levels){
        for(m in m.levels){
          if(k*m>km.max) next
          ctrl.cdkm <- ctrl.cd
          ctrl.cdkm$CD.method <- "Robbins-Monro"
          ctrl.cdkm$CD.RM.a0 <- rm.conf$a0
          ctrl.cdkm$CD.RM.c <- rm.conf$c
          ctrl.cdkm$CD.nsteps <- k
          ctrl.cdkm$CD.multiplicity <- m
          ctrl.cdkm$CD.nsteps.obs <- k
          ctrl.cdkm$CD.multiplicity.obs <- m
          limitmyps(PSLIMIT)
          mcparallel({study.fit(name, formula, response=response, reference=reference, constraints=constraints, estimate="CD", ctrl.cd=ctrl.cdkm, ctrl.simulate=ctrl.simulate, nsim=nsim, seed=rep)}, detached=TRUE)
        }
      }
    }
    
    
    # Start the CD MCMLE with the specified parameters.
    for(metric in cd.metrics){
      for(steplen.margin in steplen.margins){
        for(k in k.levels){
          for(m in m.levels){
            if(k*m>km.max) next
            ctrl.cdkm <- ctrl.cd
            ctrl.cdkm$CD.metric <- metric
            ctrl.cdkm$CD.steplength.margin <- steplen.margin
            ctrl.cdkm$CD.nsteps <- k
            ctrl.cdkm$CD.multiplicity <- m
            ctrl.cdkm$CD.nsteps.obs <- k
            ctrl.cdkm$CD.multiplicity.obs <- m
            limitmyps(PSLIMIT)
            mcparallel({study.fit(name, formula, response=response, reference=reference, constraints=constraints, estimate="CD", ctrl.cd=ctrl.cdkm, ctrl.simulate=ctrl.simulate, nsim=nsim, seed=rep)}, detached=TRUE)
          }
        }
      }
    }
  })
}

# Fit an ERGM with specified parameters, perform a simulation to
# obtain the measures of the quality of starting values.
study.fit <- function(name, formula, response, reference, constraints, estimate, ctrl.cd, ctrl.simulate, nsim, seed){
  if(estimate=="intercept"){
    nw <- ergm.getnetwork(formula)
    m<-ergm.getmodel(formula, nw, response=response,role="target")
    cn <- unlist(lapply(m$terms, function(t) NVL(names(t$params), t$coef.names)))
    o <- offset.info.formula(formula, response=response)
    epos <- which(cn=="edges" & !o$theta)
    
    theta0 <- rep(0,length(cn))
    if(is.null(response) && constraints!=~. && length(epos)){
      d <- network.edgecount(nw, na.omit=TRUE)/network.dyadcount(nw, na.omit=TRUE)
      epar <- log(d/(1-d))
      theta0[epos] <- epar    
    }else theta0 <- NULL
    
    fit <- list(coef=theta0, est="intercept")
  }else{
    fit <- time.ergm(formula, response=response, reference=reference, constraints=constraints, estimate=estimate, control=ctrl.cd, verbose=FALSE, seed=seed)
    fit <- list(est=fit$est, timings=fit$timings, coef=fit$coef, k=fit$ctrl$CD.nsteps, m=fit$ctrl$CD.multiplicity, metric=fit$ctrl$CD.metric, steplen.margin=fit$ctrl$CD.steplength.margin, iterations=fit$iterations,seed=fit$seed, RM.a0=fit$ctrl$CD.RM.a0, RM.c=fit$ctrl$CD.RM.c, method=fit$ctrl$CD.method)
  }
  
  if(!is.null(fit$coef)){
    gs <- summary(formula, response=response)
    m <- ergm.getmodel(formula, ergm.getnetwork(formula), response=response, role="target")
    obs.stats <- ergm:::.ergm.esteq(fit$coef, m, rbind(gs))
    
    cat("===== Simulating full =====\n")
    print(fit)
    set.seed(fit$seed)
    sim.stats.long <- simulate(formula, response=response, reference=reference, coef=fit$coef, constraints=constraints, esteq=TRUE, control=ctrl.simulate, nsim=nsim, statsonly=TRUE)
    cat("===== End simulating full =====\n")
    
    cat("===== Simulating short =====\n")
    print(fit)
    set.seed(fit$seed)
    sim.stats.short <- simulate(formula, response=response, reference=reference, coef=fit$coef, constraints=constraints, esteq=TRUE, control=modifyList(ctrl.simulate, list(MCMC.burnin=nsim*8, MCMC.interval=8)), nsim=nsim, statsonly=TRUE)
    cat("===== End simulating short =====\n")
    
    cat("===== Fitting ERGM: IS =====\n")
    fit.IS <- time.ergm(formula, response=response, reference=reference, constraints=constraints, estimate="MLE",
                        control=modifyList(ctrl.cd, list(init=fit$coef, MCMLE.metric="naive")), verbose=FALSE, seed=seed)
    cat("===== End fitting ERGM: IS =====\n")
    
    cat("===== Fitting ERGM: Lognormal =====\n")
    fit.LN <- time.ergm(formula, response=response, reference=reference, constraints=constraints, estimate="MLE",
                        control=modifyList(ctrl.cd, list(init=fit$coef, MCMLE.metric="lognormal")), verbose=FALSE, seed=seed)
    cat("===== End fitting ERGM: Lognormal =====\n")
    
    
    # Store the fit information.
    fit <- c(fit, list(Hummel1Step=ERRVL(try(ergm:::.Hummel.steplength(sim.stats.long, obs.stats)), 0),
                       Hummel1Step.short=ERRVL(try(ergm:::.Hummel.steplength(sim.stats.short, obs.stats)), 0),
                       MLE.IS=fit.IS,
                       MLE.LN=fit.LN))
  }
  
  save(fit, file=paste(name, Sys.time(), Sys.getpid(), round(runif(1)*1e10), "RData", sep="."))
}


run.study("kap.study", kapferer ~ edges + gwdegree(0.25, fixed = TRUE) + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE))
run.study("zach.binom.study", zach ~ nonzero+sum+nodefactor("role",base=2)+transitiveweights, reference = ~Binomial(8), response = "contexts", ctrl.cd=control.ergm(CD.maxit=60,CD.trustregion=100000,MCMLE.trustregion=100000))
run.study("zach.cmp.study", zach ~ nonzero+CMP+sum+nodefactor("role",base=2)+transitiveweights, reference  = ~Poisson, response = "contexts", ctrl.cd=control.ergm(CD.maxit=60,CD.trustregion=100000,MCMLE.trustregion=100000))
run.study("kap2.study", kapferer ~ edges + gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE))
run.study("ecoli2.study", ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE))
run.study("ecoli2.self.study", ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = TRUE) + nodemix("self", base = 1))
run.study("laz2.c.study", laz ~ gwesp()+nodecov("seniority")+nodefactor("practice")+nodematch("practice")+nodematch("gender")+nodematch("office"), constraints = ~edges)
run.study("laz2.study", laz ~ edges+gwesp()+nodecov("seniority")+nodefactor("practice")+nodematch("practice")+nodematch("gender")+nodematch("office"))
