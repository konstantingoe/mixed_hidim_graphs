#### Extended Comparison out Method and Fan et al 

#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

sim <- 100 # simulation runs
n <- c(200,200,300,600) # sample size include high dimension only on cluster
d <- c(50,250,750, 1500) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- c(200,250,750,1500) # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 50 # number of tuning parameters for graphical lasso
countvar <- F
latent <- F
matexport = F

cat(paste0("number of available cores: ",detectCores()))

if (detectCores() >= 100){
  numCores <-  100
} else {
  numCores <- detectCores()
}

run <- F
if (run == T) {
print("Start with d=50")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

result_1 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[1], d=d[1], nlam=nlam, matexport = F))

plan(sequential)

table_1 <- kendalls.results(result_1)
stargazer(table_1, out = "table_kendall_1.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    


print("Start with d=250")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

result_2 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[2], d=d[2], nlam=nlam, matexport = F))

plan(sequential)

table_2 <- kendalls.results(result_2)
stargazer(table_2, out = "table_kendall_2.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))


print("Start with d=750")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster


result_3 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[3], d=d[3], nlam=nlam, matexport = F))

plan(sequential)

table_3 <- kendalls.results(result_3)
stargazer(table_3, out = "table_kendall_3.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))

}

print("Start with d=1500")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

result_4 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[4], d=d[4], mode = "fan", nlam=nlam, matexport = F))

plan(sequential)

table_4 <- kendalls.results(result_4)
stargazer(table_4, out = "table_kendall_4.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[4],"and d=",d[4],"under",sim, "simulation runs."))
