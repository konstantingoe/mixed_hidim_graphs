#### Extended Comparison out Method and Fan et al 

#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

sim <- 100 # simulation runs
n <- c(200,200,300) # sample size include high dimension only on cluster
d <- c(50,250,750) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- 200 # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 50 # number of tuning parameters for graphical lasso
countvar <- F
latent <- F

cat(paste0("number of available cores: ",detectCores()))

if (detectCores() >= 100){
  numCores <-  100
} else {
  numCores <- detectCores()
}

print("Start with d=50")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

result_1 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[1], d=d[1], n_E = n_E, nlam=nlam, matexport = F))

plan(sequential)

table_1 <- kendalls.results(result_1)
stargazer(table_1, out = "table_kendall_1.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    


print("Start with d=250")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

result_2 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[2], d=d[2], n_E = n_E, nlam=nlam, matexport = F))

plan(sequential)

table_2 <- kendalls.results(result_2)
stargazer(table_2, out = "table_kendall_2.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))


print("Start with d=750")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster


result_3 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[3], d=d[3], n_E = n_E, nlam=nlam, matexport = F))

plan(sequential)

table_3 <- kendalls.results(result_3)
stargazer(table_3, out = "table_kendall_3.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))

