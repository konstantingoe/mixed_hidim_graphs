#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

sim <- 50 #100 # simulation runs
n <- c(200,200,300,600) # sample size include high dimension only on cluster
d <- c(50,250,750, 1500) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- c(200,250,750,1500) # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 30 #50 # number of tuning parameters for graphical lasso
countvar <- F
latent <- F
matexport <-  F
first_half <- T


print(paste0("number of available cores: ",detectCores()))
if (detectCores() >= 100){
  numCores <-  100
} else {
  numCores <- detectCores()
}

if (first_half == F){
print("Start with d=50")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

mixed_result_1 <- future_lapply(future.seed = T, 1:sim, function(k) 
                    serverrun(n=n[1], d=d[1], latent = T, nlam=nlam, matexport = F,
                      namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                        unbalanced = .5))
plan(sequential)

table_1 <- results_generator(mixed_result_1)
stargazer(table_1, out = "table_1.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    


print("Start with d=250")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

mixed_result_2 <- future_lapply(future.seed = T, 1:sim, function(k) 
                    serverrun(n=n[2], d=d[2], latent = T, nlam=nlam, matexport = F,
                      namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                        unbalanced = .5))
plan(sequential)

table_2 <- results_generator(mixed_result_2)

stargazer(table_2, out = "table_2.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))


print("Start with d=750")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

mixed_result_3 <- future_lapply(future.seed = T, 1:sim, function(k) 
                    serverrun(n=n[3], d=d[3], latent = T, nlam=nlam, matexport = F,
                      namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                        unbalanced = .5))
plan(sequential)

table_3 <- results_generator(mixed_result_3)
stargazer(table_3, out = "table_3.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))

}
print("Start with d=1500")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

mixed_result_4 <- future_lapply(future.seed = T, 1:sim, function(k) 
                    serverrun(n=n[4], d=d[4], latent = T, nlam=nlam, matexport = F,
                      namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                        unbalanced = .5))
plan(sequential)

table_4 <- results_generator(mixed_result_4)
stargazer(table_4, out = "table_4.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[4],"and d=",d[4],"under",sim, "simulation runs."))




