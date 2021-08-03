rm(list = ls())
source("packages.R")
source("functions.R")
source("applied_functions.R")

set.seed(1221)

sim <- 100 # simulation runs
n <- c(200,200,300,600) # sample size include high dimension only on cluster
d <- c(50,250,750, 1500) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- c(200,250,750,1500) # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 50 # number of tuning parameters for graphical lasso
param <- .1


X <- rmvnorm(1000, mean = c(-40,40), sigma = matrix(c(1,-.98,-.98,1), nrow = 2, ncol = 2))
cov(X[,1],X[,2])
cov(X[,1][order(X[,1])],X[,2][order(X[,2])]) 

var(X[,1][order(X[,1])])
var(X[,1])


cov(X[,1][order(X[,1])],X[,2][order(X[,2])]) < var(X[,1][order(X[,1])])
