
#rm(list = ls())

# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  huge,
  glasso,
  stats,
  MASS,
  gdata,
  Matrix,
  polycor,
  corpcor)

source("functions.R")

set.seed(1221)
param <- .25 #dimensionality penalty in eBIC later
n <- 200 # sample size
d <- 750 # dimensionality
nlam <- 50 # number of lambda value for huge routine


### here only use d for which I have chosen appropriate c, i.e.
### d = c(50, 250, 750, 1500, 3000) all output |E| = 200
data <- generate.data(d = d, n = n) # generate according to Fan et al
data_0 <- data[[1]] # data matrix
Omega <- data[[2]] # true precision matrix


#input: continuous data from which ordinal data will be generated 
# proportion of data that should be converted to ordinal
# default: proportion = .5
# number of levels for all ordinal variables (could also think of a vector of size proportion*d here)
# n_O = 3 default
# allow for a mix of 20-leveled discrete variables (count variables) and 3-leveled ones
data_mixed <- make.ordinal(data = data_0, countvar = T) 

# polychoric / polyserial approach to sample correlation matrix rho
rho <- mixed.omega(data_mixed)
#system.time(rho2 <- mixed.omega(data_mixed))

### glasso with sample correlation rho 
huge.result <- huge(rho,nlambda=nlam,method="glasso",verbose=FALSE)

### threshold selection approach 
Omega_hat <- omega.select(x=huge.result, n=n)
### best model selection with huge and then likelihood with glasso package and eBIC
Omega_hat_2 <- omega.select.drton(x=huge.result, n=n, s = rho,  param = 0)

# true number of edges: 
edgenumber(Omega)
# estimated number of edges with thresholding:
edgenumber(Omega_hat)
# estimated number of edges with 2 stage method:
edgenumber(Omega_hat_2)


### performance 
frobenius1 <- base::norm(Omega_hat - Omega, type = "F")
frobenius2 <- base::norm(Omega_hat_2 - Omega, type = "F")

adj_estimate <- abs(Omega_hat) > 0
adj_estimate_2 <- abs(Omega_hat_2) > 0

tpr_1 <- tpr(truth = Omega, estimate = adj_estimate)
tpr_2 <- tpr(truth = Omega, estimate = adj_estimate_2)

fpr_1 <- fpr(truth = Omega, estimate = adj_estimate)
fpr_2 <- fpr(truth = Omega, estimate = adj_estimate_2)






