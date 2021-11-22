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
nlam <-30 # 50 # number of tuning parameters for graphical lasso
param <- .1
firstrun <- F

print(paste0("number of available cores: ",availableCores()))
if (availableCores() >= 100){
  numCores <-  100
} else {
  numCores <- availableCores()
}

if (firstrun == T){
##### d = 50 ####

print("Start with d=50, f_j(x) = x")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_1 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[1], d = d[1], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 1))

plan(sequential)

table_1 <- extract.kendall.nonpararesults(nonpara_comparison_1)
stargazer(table_1, out = "table_1_binary.tex", summary = F, title=paste("Mixed binary data structure learning comparison n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    

print("continue with d=50, f_j(x) = x^3")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_2 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[1], d = d[1], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 3))
plan(sequential)

table_2 <- extract.kendall.nonpararesults(nonpara_comparison_2)
stargazer(table_2, out = "table_2_binary.tex", summary = F, title=paste("Mixed binary data structure learning with $f_j(x) = x^3$, n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    

##### d = 250 ####

print("Start with d=250, f_j(x) = x")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_3 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[2], d = d[2], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 1))
plan(sequential)

table_3 <- extract.kendall.nonpararesults(nonpara_comparison_3)
stargazer(table_3, out = "table_3_binary.tex", summary = F, title=paste("Mixed binary data structure learning comparison n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))                    

print("continue with d=250, f_j(x) = x^3")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_4 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[2], d = d[2], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 3))
plan(sequential)

table_4 <- extract.kendall.nonpararesults(nonpara_comparison_4)
stargazer(table_4, out = "table_4_binary.tex", summary = F, title=paste("Mixed binary data structure learning with $f_j(x) = x^3$, n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))                    

##### d = 750 #### 

print("Start with d=750, f_j(x) = x")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_5 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[3], d = d[3], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 1))
plan(sequential)

table_5 <- extract.kendall.nonpararesults(nonpara_comparison_5)
stargazer(table_5, out = "table_5_binary.tex", summary = F, title=paste("Mixed binary data structure learning comparison n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))                    

print("continue with d=750, f_j(x) = x^3")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_6 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[3], d = d[3], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 3))
plan(sequential)

table_6 <- extract.kendall.nonpararesults(nonpara_comparison_6)
stargazer(table_6, out = "table_6_binary.tex", summary = F, title=paste("Mixed binary data structure learning with $f_j(x) = x^3$, n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))                    

##### d = 1500 #### 

print("Start with d=1500, f_j(x) = x")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_7 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[4], d = d[4], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 1))
plan(sequential)

table_7 <- extract.kendall.nonpararesults(nonpara_comparison_7)
stargazer(table_7, out = "table_7_binary.tex", summary = F, title=paste("Mixed binary data structure learning comparison n=",n[4],"and d=",d[4],"under",sim, "simulation runs."))                    
}
print("continue with d=1500, f_j(x) = x^3")

plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

nonpara_comparison_8 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall.nonpara(t=t, n = n[4], d = d[4], 
                                                                                                    nlam=nlam, matexport = F, countvar = F, 
                                                                                                    mode = "fan", param = .1, f_j = 3))
plan(sequential)

table_8 <- extract.kendall.nonpararesults(nonpara_comparison_8)
stargazer(table_8, out = "table_8_binary.tex", summary = F, title=paste("Mixed binary data structure learning with $f_j(x) = x^3$, n=",n[4],"and d=",d[4],"under",sim, "simulation runs."))                    


