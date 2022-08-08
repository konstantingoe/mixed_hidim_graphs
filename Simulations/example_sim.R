rm(list = ls())
source("./Packages/packages.R")
source("./Functions/functions.R")
#source("applied_functions.R")

set.seed(13311)

sim <- 20 # simulation runs
n <- 50 # sample size
d <-  50 # dimensionality 
n_E <- 50
t <- .15 # signal strength
nlam <- 30 #50 # number of tuning parameters for graphical lasso
numCores <- 2


plan(multisession, workers = numCores) ## Run in parallel on Linux cluster
  
nonpara_comparison_1 <- future_lapply(future.seed = T, 1:sim, function(k) 
                        ternary_run(n=n, d=d, nlam=nlam, matexport = F,
                        namevector = c("binary" = T, "ordinal" = T, "poisson" = F),
                        unbalanced = .5, mode = "fan", nonpara = T))
  
plan(sequential)

extract.ternary.results(nonpara_comparison_1)


  