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

n <- 200
X <- rmvnorm(200, mean = c(0,0), sigma = matrix(c(1,.2,.2,1), nrow = 2, ncol = 2))
#X <- X^3
cov(X[,1],X[,2])
cov(X[,1][order(X[,1])],X[,2][order(X[,2])]) 

lambda <- 2
lambda_sample <- abs(rnorm(1,mean = lambda, sd = 2))
ordinal_poisson <- qpois(pnorm(scale(X[,2])),lambda = lambda_sample)
table(ordinal_poisson)

cov(X[,1],ordinal_poisson)
cov(X[,1][order(X[,1])],ordinal_poisson[order(ordinal_poisson)])
polyserial(X[,1],ordinal_poisson)
lord_nonparanormal(X[,1],ordinal_poisson)
spearman(X[,1],ordinal_poisson)


###check the relations 
testing <- ((sum(X[,1])/(n+1)) + (1 - (sum(X[,2])/(n+1))))/n+2 - (sum(X[,1])^2/n*sum(X[,2])^2/n)
cov(X[,1][order(X[,1])],X[,2][order(X[,2])]) 
n*cov(X[,1],X[,2])



