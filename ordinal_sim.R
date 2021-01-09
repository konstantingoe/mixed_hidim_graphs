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
nlam <- 30 # number of tuning parameters for graphical lassols()
plan(multisession, workers = (availableCores() - 10)) ## Run in parallel on local computer
options(future.globals.maxSize= 20000*1024^2)


#### data generation via sparse Omega ####

data <- setNames(lapply(seq_along(n), function(i)
          setNames(future_lapply(future.seed = T, 1:sim, function(k)
            generate.data(t=t,n=n[i],d=d[i], n_E = n_E)),nm=1:sim)),nm=paste("d =",d))

# data according to Multivariate normal 
data_0 <- setNames(lapply(seq_along(n), function(i)
            setNames(lapply(1:sim, function(k) 
              data[[i]][[k]][[1]]),nm=1:sim)),nm=paste("d =",d)) 

# underlying undirected graph via precision matrix
Omega <- setNames(lapply(seq_along(n), function(i)
            setNames(lapply(1:sim, function(k) 
              data[[i]][[k]][[2]]),nm=1:sim)),nm=paste("d =",d)) 

#### choose d_1 Variables to be ordinal and let's give them all 3 categories

data_mixed <- setNames(lapply(seq_along(n), function(i)
                setNames(lapply(1:sim, function(k) 
                  make.ordinal(data=data_0[[i]][[k]])),nm=1:sim)),nm=paste("d =",d)) 

### benchmark against unknown latent data! ####
save(data, file = "data.Rdata")
save(data_0, file = "data_0.Rdata")
save(data_mixed, file = "data_mixed.Rdata")
save(Omega, file = "Omega.Rdata")
save(rho_latent, file = "rho_latent.Rdata")
save(rho_pd, file = "rho_pd.Rdata")

rho_pd <- setNames(lapply(seq_along(n), function(i)
            setNames(future_lapply(future.seed = T, 1:sim, function(k) 
              mixed.omega(data = data_mixed[[i]][[k]])),nm=1:sim)),nm=paste("d =",d)) 

rho_latent <- setNames(lapply(seq_along(n), function(i)
                setNames(future_lapply(future.seed = T, 1:sim, function(k) 
                  mixed.omega(data = data_0[[i]][[k]])),nm=1:sim)),nm=paste("d =",d)) 

rm(data, data_0, data_mixed, Omega, rho_latent)
### perform glasso ####  

results_hat <- setNames(lapply(seq_along(n), function(i)
                setNames(future_lapply(future.seed = T, 1:sim, function(k) 
                  glasso.results(Sigma = rho_pd[[i]][[k]], Omega = Omega[[i]][[k]],
                   nlam=nlam, n=n[i])),nm=1:sim)),nm=paste("d =",d)) 

results_latent <- setNames(lapply(seq_along(n), function(i)
                    setNames(future_lapply(future.seed = T, 1:sim, function(k) 
                      glasso.results(Sigma = rho_latent[[i]][[k]], Omega = Omega[[i]][[k]],
                       nlam=nlam, n=n[i])),nm=1:sim)),nm=paste("d =",d)) 










### result object would have 410 GB... not feasible write function so that result object can be deleted afterwards! 
#sim=2

result <- setNames(lapply(seq_along(n), function(i)
            setNames(future_lapply(future.seed = T, 1:sim, function(k) 
              huge(rho_pd[[i]][[k]],nlambda=nlam,method="glasso",verbose=FALSE)),nm=1:sim)),nm=paste("d =",d)) 
            
result.benchmark <- setNames(lapply(seq_along(n), function(i)
                      setNames(future_lapply(future.seed = T, 1:sim, function(k) 
                        huge(rho_latent[[i]][[k]],nlambda=nlam,method="glasso",verbose=FALSE)),nm=1:sim)),nm=paste("d =",d))  

#### use eBIC proposed by Foygel & Drton (2010)  
Omega_hat <- setNames(lapply(seq_along(n), function(i)
              setNames(future_lapply(future.seed = T, 1:sim, function(k) 
                omega.select(x=result[[i]][[k]], n=n[i])),nm=1:sim)),nm=paste("d =",d)) 

Omega_Z <- setNames(lapply(seq_along(n), function(i)
            setNames(future_lapply(future.seed = T, 1:sim, function(k) 
              omega.select(x=result.benchmark[[i]][[k]], n=n[i])),nm=1:sim)),nm=paste("d =",d))  

#### obtain adjacency matrix ####   
adj_1 <- setNames(lapply(seq_along(n), function(i)
          setNames(lapply(1:sim, function(k) 
            (abs(Omega_hat[[i]][[k]]) > .05)),nm=1:sim)),nm=paste("d =",d))  

adj_z <- setNames(lapply(seq_along(n), function(i)
          setNames(lapply(1:sim, function(k) 
            (abs(Omega_Z[[i]][[k]]) > .05)),nm=1:sim)),nm=paste("d =",d))    

adj_0 <- setNames(lapply(seq_along(n), function(i)
          setNames(lapply(1:sim, function(k) 
            (Omega[[i]][[k]] > 0)),nm=1:sim)),nm=paste("d =",d))   
  
### calculate matrix norm ####
frobenius_hat <- setNames(lapply(seq_along(n), function(i)
                  sapply(1:sim, function(k) 
                    base::norm((Omega_hat[[i]][[k]]-Omega[[i]][[k]]), type = "F"))),nm=paste("d =",d))  

frobenius_z <- setNames(lapply(seq_along(n), function(i)
                  sapply(1:sim, function(k) 
                    base::norm((Omega_Z[[i]][[k]]-Omega[[i]][[k]]), type = "F"))),nm=paste("d =",d))    

##### false and true positive rate ####  
fpr_hat <- setNames(lapply(seq_along(n), function(i)
            sapply(1:sim, function(k) 
              fpr(truth = adj_0[[i]][[k]], estimate = adj_1[[i]][[k]]))),nm=paste("d =",d))  
            
tpr_hat <- setNames(lapply(seq_along(n), function(i)
            sapply(1:sim, function(k) 
              tpr(truth = adj_0[[i]][[k]], estimate = adj_1[[i]][[k]]))),nm=paste("d =",d))  
  
fpr_z <- setNames(lapply(seq_along(n), function(i)
            sapply(1:sim, function(k) 
              fpr(truth = adj_0[[i]][[k]], estimate = adj_z[[i]][[k]]))),nm=paste("d =",d)) 
            
tpr_z <- setNames(lapply(seq_along(n), function(i)
            sapply(1:sim, function(k) 
              tpr(truth = adj_0[[i]][[k]], estimate = adj_z[[i]][[k]]))),nm=paste("d =",d)) 

#### pAUC #####
pAUC_hat <- setNames(lapply(seq_along(n), function(i)
              sapply(1:sim, function(k) 
                pAUC(truth = Omega[[i]][[k]], huge_obj = result[[i]][[k]]))),nm=paste("d =",d)) 
              
pAUC_z <- setNames(lapply(seq_along(n), function(i)
            sapply(1:sim, function(k) 
              pAUC(truth = Omega[[i]][[k]], huge_obj = result.benchmark[[i]][[k]]))),nm=paste("d =",d)) 


#### results in table #####
table <- list()
for (i in seq_along(n)){
table[[i]] <- round(as_tibble(rbind(c('polychoric' = mean(frobenius_hat[[i]]),'sd_p' = sd(frobenius_hat[[i]]),'latent data'= mean(frobenius_z[[i]]), 'sd_l' = sd(frobenius_z[[i]])),
                c(mean(fpr_hat[[i]]),sd(fpr_hat[[i]]), mean(fpr_z[[i]]),sd(fpr_z[[i]])),
                c(mean(tpr_hat[[i]]),sd(tpr_hat[[i]]), mean(tpr_z[[i]]),sd(tpr_z[[i]])),
                c(mean(pAUC_hat[[i]]),sd(pAUC_hat[[i]]), mean(pAUC_z[[i]]),sd(pAUC_z[[i]])))),4)
stargazer(table[[i]], summary = F, title=paste("Mixed data structure learning of the precision matrix with n=",n[i],"and d=",d[i],"under",sim, "simulation runs."))

}


