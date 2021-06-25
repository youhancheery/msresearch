# compare the EE implementation in ergm to the EE implementation by Stival
# from Stivala's end we will implement his python polblogs example separately to this
# here we will just use ergm EE/SA/RM for the estimation of polblogs

################################################## 
# Run this whenever changes have been made to ergm
installr::uninstall.packages("ergm")
devtools::install_github("youhancheery/ergm")
library(ergm)
##################################################



