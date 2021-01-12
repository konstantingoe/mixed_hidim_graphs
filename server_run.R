#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

#### sample from multivariate normal
set.seed(1234)


sim = 100  # simulation runs
n <- c(200,200,600) # sample size include high dimension only on cluster
d <- c(50,250,3000) # dimensionality --> include high dimension (3000) only on cluster 
n_E <- 200 # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 50 # number of tuning parameters for graphical lassols()
plan(multisession, workers = 20) ## Run in parallel on local computer


mixed_result_1 <- future_lapply(future.seed =T, 1:sim, function(k) serverrun(n=n[1], d=d[1], n_E = n_E, latent = F, nlam=nlam, matexport = F))

latent_result_1 <- future_lapply(future.seed =T, 1:sim, function(k) serverrun(n=n[1], d=d[1], n_E = n_E, latent = T, nlam=nlam, matexport = F))             






