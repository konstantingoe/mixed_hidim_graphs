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

kendall_result_1 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[1], d=d[1], n_E = n_E, method = "kendall", nlam=nlam, matexport = F))
poly_result_1 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[1], d=d[1], n_E = n_E, method = "poly", nlam=nlam, matexport = F))             

plan(sequential)

table_1 <- round(as_tibble(rbind(c('kendall' = mean(extract.result(kendall_result_1, which = "F")),'sd_p' = sd(extract.result(kendall_result_1, which = "F")),
                                   'polychoric'= mean(extract.result(poly_result_1, which = "F")), 'sd_l' = sd(extract.result(poly_result_1, which = "F"))),
                                 
                                 c(mean(extract.result(kendall_result_1, which = "FPR")),sd(extract.result(kendall_result_1, which = "FPR")), 
                                   mean(extract.result(poly_result_1, which = "FPR")),sd(extract.result(poly_result_1, which = "FPR"))),
                                 
                                 c(mean(extract.result(kendall_result_1, which = "TPR")),sd(extract.result(kendall_result_1, which = "TPR")), 
                                   mean(extract.result(poly_result_1, which = "TPR")),sd(extract.result(poly_result_1, which = "TPR"))),
                                 
                                 c(mean(extract.result(kendall_result_1, which = "AUC")),sd(extract.result(kendall_result_1, which = "AUC")),
                                   mean(extract.result(poly_result_1, which = "AUC")),sd(extract.result(poly_result_1, which = "AUC"))))),4)

stargazer(table_1, out = "table_kendall_1.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    

print("Start with d=250")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster

kendall_result_2 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[2], d=d[2], n_E = n_E, method = "kendall", nlam=nlam, matexport = F))
poly_result_2 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[2], d=d[2], n_E = n_E, method = "poly", nlam=nlam, matexport = F))             

plan(sequential)

table_2 <- round(as_tibble(rbind(c('polychoric' = mean(extract.result(kendall_result_2, which = "F")),'sd_p' = sd(extract.result(kendall_result_2, which = "F")),
                                   'latent data'= mean(extract.result(poly_result_2, which = "F")), 'sd_l' = sd(extract.result(poly_result_2, which = "F"))),
                                 
                                 c(mean(extract.result(kendall_result_2, which = "FPR")),sd(extract.result(kendall_result_2, which = "FPR")), 
                                   mean(extract.result(poly_result_2, which = "FPR")),sd(extract.result(poly_result_2, which = "FPR"))),
                                 
                                 c(mean(extract.result(kendall_result_2, which = "TPR")),sd(extract.result(kendall_result_2, which = "TPR")), 
                                   mean(extract.result(poly_result_2, which = "TPR")),sd(extract.result(poly_result_2, which = "TPR"))),
                                 
                                 c(mean(extract.result(kendall_result_2, which = "AUC")),sd(extract.result(kendall_result_2, which = "AUC")),
                                   mean(extract.result(poly_result_2, which = "AUC")),sd(extract.result(poly_result_2, which = "AUC"))))),4)

stargazer(table_2, out = "table_kendall_2.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[2],"and d=",d[2],"under",sim, "simulation runs."))


print("Start with d=750")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster


kendall_result_3 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[3], d=d[3], n_E = n_E, method = "kendall", nlam=nlam, matexport = F))
poly_result_3 <- future_lapply(future.seed = T, 1:sim, function(k) serverrun.kendall(n=n[3], d=d[3], n_E = n_E, method = "poly", nlam=nlam, matexport = F))             

plan(sequential)

table_3 <- round(as_tibble(rbind(c('polychoric' = mean(extract.result(kendall_result_3, which = "F")),'sd_p' = sd(extract.result(kendall_result_3, which = "F")),
                                   'latent data'= mean(extract.result(poly_result_3, which = "F")), 'sd_l' = sd(extract.result(poly_result_3, which = "F"))),
                                 
                                 c(mean(extract.result(kendall_result_3, which = "FPR")),sd(extract.result(kendall_result_3, which = "FPR")), 
                                   mean(extract.result(poly_result_3, which = "FPR")),sd(extract.result(poly_result_3, which = "FPR"))),
                                 
                                 c(mean(extract.result(kendall_result_3, which = "TPR")),sd(extract.result(kendall_result_3, which = "TPR")), 
                                   mean(extract.result(poly_result_3, which = "TPR")),sd(extract.result(poly_result_3, which = "TPR"))),
                                 
                                 c(mean(extract.result(kendall_result_3, which = "AUC")),sd(extract.result(kendall_result_3, which = "AUC")),
                                   mean(extract.result(poly_result_3, which = "AUC")),sd(extract.result(poly_result_3, which = "AUC"))))),4)


stargazer(table_3, out = "table_kendall_3.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))


output <- list("Mixed with d = 50" = kendall_result_1, "Mixed with d = 250" = kendall_result_2, "Mixed with d = 1500" = kendall_result_3,
               "Latent with d = 50" = poly_result_1, "Latent with d = 250" = poly_result_2, "Latent with d = 1500" = poly_result_3)

save(output, file = "output_kendall.Rdata")

