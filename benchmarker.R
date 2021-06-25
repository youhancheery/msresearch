# test the ergm package works with new changes
# 1. load a ergm dataset
# 2. test EE algorithm with ergm dataset

################################################## 
# Run this whenever changes have been made to ergm
installr::uninstall.packages("ergm")
devtools::install_github("youhancheery/ergm")
library(ergm)
##################################################

library(dplyr)
library(ggplot2)

sd_estimate <- function(net_est) {
  # accepts a network variance-covariance matrix
  return(sqrt(diag(vcov(net_est, sources="estimation"))))
  
}

# florentine marriages data
data(flo)
flomarriage <- network(flo,directed=FALSE)
flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)

# plot(flomarriage)
# plot(flomarriage, vertex.cex=flomarriage %v% "wealth" / 20, main="Marriage Ties")

# estimate with stochapprox
sa_gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             eval.loglik=FALSE,
             control=control.ergm(main.method = "Stochastic-Approximation",
                                  MCMLE.maxit=20))

# estimate with equilibrium expectation
ee_gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
             eval.loglik=FALSE,
             control=control.ergm(main.method = "EE"))

# calculate the standard deviation between the estimates
ee_sd <- sd_estimate(ee_gest)
sa_sd <- sd_estimate(sa_gest)

# molecular data
data("molecule")

molecule_sa_est <- ergm(molecule ~ edges + degree(1),
                        control = control.ergm(main.method = "Stochastic-Approximation"))
molecule_ee_est <- ergm(molecule ~ edges + degree(1),
                        control = control.ergm(main.method = "EE"))
sa_mol_sd <- sd_estimate(molecule_sa_est)
ee_mol_sd <- sd_estimate(molecule_ee_est)

# magnolia data
data("faux.magnolia.high")

magnolia_sa_est <- ergm(molecule ~ edges + degree(1:2),
                        control = control.ergm(main.method = "Stochastic-Approximation"))
magnolia_ee_est <- ergm(molecule ~ edges + degree(1:2),
                        control = control.ergm(main.method = "EE"))
(sa_magnolia_sd <- sd_estimate(magnolia_sa_est))
(ee_magnolia_sd <- sd_estimate(magnolia_ee_est))

