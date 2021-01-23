#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

sim <- 10 # simulation runs
n <- c(200,200,600) # sample size include high dimension only on cluster
d <- c(50,250,1500) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- 200 # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 50 # number of tuning parameters for graphical lasso

print("starting")

mixed_result_3 <- lapply(1:sim, function(k) serverrun(n=n[3], d=d[3], n_E = n_E, latent = F, nlam=nlam, matexport = F, countvar = T))

print("mixed_3 done, start with latent_3")
latent_result_3 <- lapply(1:sim, function(k) serverrun(n=n[3], d=d[3], n_E = n_E, latent = T, nlam=nlam, matexport = F))             

table_3 <- round(as_tibble(rbind(c('polychoric' = mean(extract.result(mixed_result_3, which = "F")),'sd_p' = sd(extract.result(mixed_result_3, which = "F")),
                                   'latent data'= mean(extract.result(latent_result_3, which = "F")), 'sd_l' = sd(extract.result(latent_result_3, which = "F"))),
                                 
                                 c(mean(extract.result(mixed_result_3, which = "FPR")),sd(extract.result(mixed_result_3, which = "FPR")), 
                                   mean(extract.result(latent_result_3, which = "FPR")),sd(extract.result(latent_result_3, which = "FPR"))),
                                 
                                 c(mean(extract.result(mixed_result_3, which = "TPR")),sd(extract.result(mixed_result_3, which = "TPR")), 
                                   mean(extract.result(latent_result_3, which = "TPR")),sd(extract.result(latent_result_3, which = "TPR"))),
                                 
                                 c(mean(extract.result(mixed_result_3, which = "AUC")),sd(extract.result(mixed_result_3, which = "AUC")),
                                   mean(extract.result(latent_result_3, which = "AUC")),sd(extract.result(latent_result_3, which = "AUC"))))),4)


stargazer(table_3, out = "table_3.tex", summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[3],"and d=",d[3],"under",sim, "simulation runs."))
