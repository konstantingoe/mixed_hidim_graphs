rm(list = ls())
source("packages.R")
source("functions.R")
source("applied_functions.R")

set.seed(1221)

sim <- 50 #100 # simulation runs
n <- 200 # sample size include high dimension only on cluster
d <- 50 # dimensionality --> include high dimension (1500) only on cluster 
n_E <- c(200,250,750,1500) # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 30 #50 # number of tuning parameters for graphical lasso

nonparanormal_run <- function(n=n, d=d, nlam=nlam, matexport = F,
                              namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                              unbalanced = .5, low = .05, high = .1, sparsity = .1){
  data <- generate.data(t=t, n = 200, d = 50, mode = "fan")
  data_0 <- data[[1]]
  Omega <- data[[2]]
  data <- NULL
  data_mixed <- make.ordinal.general(data_0, namevector = namevector, unbalanced = unbalanced, low = low, high = high)
  
  rho_latent <- mixed.omega(data_0, verbose = F) #oracle
  data_0 <- NULL
  rho <- mixed.omega(data_mixed, verbose = F) #MLestimator
  rho_nonpara <- mixed.omega.paranormal(data_mixed, verbose = F) # adhoc nonparanormal estimator
  data_mixed <- NULL
  
  results_latent <- glasso.results(Sigma=rho_latent, Omega=Omega, nlam=nlam, n=n, matexport = matexport, param = .1)
  results_ml <- glasso.results(Sigma=rho, Omega=Omega, nlam=nlam, n=n, matexport = matexport,  param = .1)
  results_nonpara <- glasso.results(Sigma=rho_nonpara, Omega=Omega, nlam=nlam, n=n, matexport = matexport,  param = .1)
  
  results <- list("latent"=results_latent, "ML" = results_ml, "nonparanormal" = results_nonpara)
  return(results)
}

nonpara_comparison_1 <- future_lapply(future.seed = T, 1:sim, function(k) nonparanormal_run(n=n[1], d=d[1], nlam=nlam, matexport = F,
                          namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
                          unbalanced = .5))

table_1 <-  cbind(c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[1]][[1]])),
             sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[1]][[1]]))),
            c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[2]][[1]])),
              sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[2]][[1]]))),
            c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[3]][[1]])),
              sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[3]][[1]]))))

for (i in 2:4) {
  table_1 <- rbind(table_1,
                   cbind(c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[1]][[i]])),
                           sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[1]][[i]]))),
                         c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[2]][[i]])),
                           sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[2]][[i]]))),
                         c(mean(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[3]][[i]])),
                           sd(sapply(1:sim, function(k) nonpara_comparison_1[[k]][[3]][[i]])))))    
}

rownames(table_1) <- c("Frobenius", "sd_F", "FPR", "sd_{FPR}", "TPR", "sd_{TPR}", "AUC", "sd_{AUC}")
colnames(table_1) <- c("Oracle", "Polyserial ML", "Polyserial nonparanormal")
stargazer(table_1, out = "table_1.tex", summary = F, title=paste("Mixed data structure learning comparison n=",n[1],"and d=",d[1],"under",sim, "simulation runs."))                    


xnew <- mvrnorm(n=2000, mu = c(1,1), Sigma = matrix(c(1,-.27,-.27,1),nrow=2,ncol=2))
#xnew[,1][190:200] <- xnew[,1][190:200]*5

lambda_sample <- abs(rnorm(1,mean = 5, sd = 2))
xnew_poisson <- qpois(pnorm(scale(xnew[,2])),lambda = lambda_sample)

pbin <- runif(1,.2,.9)
xnew_bernoulli <- qbinom(pnorm(scale(xnew[,2])),size=1,prob = pbin)


lord_nonparanormal(xnew[,1],xnew_bernoulli)
polycor::polyserial(xnew[,1],xnew_bernoulli)

lord_nonparanormal(xnew[,1],xnew_poisson)
pointpolynormal(xnew[,1],xnew_poisson)
### that's a problen... use Lord generlized biserial estimator!


lord_nonparanormal <- function(x, y, maxcor = 0.9999, more_verbose = T){
  x <- if (missing(y)){ 
    x
  } else {cbind(x, y)}
  
  x <- as.data.frame(x)
  
  if (any(is.factor(x) == F)){
    if (more_verbose == T) cat("No factor variable specified. I'm taking the one that has fewer than 20 unique values!")
    factor_id <- sapply(x, function(id) length(unique(id)) < 20)
  } else {
    factor_id <- sapply(x, is.factor)
  }
  
  ### if both categorical perform polychoric correlation 
  
  if (sum(factor_id) == 2){
    lord_estimator <- polycor::polychor(x[,1], x[,2])
  } else {
    
    ### retrieve numeric and discrete variable
    numeric_var <- x[,factor_id == F]
    factor_var <- x[,factor_id == T]
    
    ranky <- rank(numeric_var)
    rankmean <- (length(ranky)+1)/2
    n <- length(factor_var)
    
    cummarg_propotions <- c(0,cumsum(table(factor_var)/n))
    sumindex <- n*cummarg_propotions
    
    s_Y <- sqrt(1/(n)*sum((ranky - rankmean)^2))
    
    a_i <- seq_along(as.numeric(levels(as.factor(factor_var))))
    b <- a_i
    for (i in a_i){
      b[i] <- a_i[i]*sum(ranky[order(ranky)[(1+sumindex[i]):sumindex[i+1]]] - rankmean)
    }
    lambda <- 1/(n*s_Y)*sum(b)
    
    s_X <- sqrt(1/(n)*sum((factor_var - mean(factor_var))^2))
    samplecorr <- spearman(numeric_var, factor_var)
    if (abs(samplecorr) > maxcor) 
      samplecorr <- sign(samplecorr) * maxcor
    
    lord_estimator <- samplecorr*s_X/lambda
  }
  if (lord_estimator < 0){
    numeric_var <- -1*numeric_var
    ranky <- rank(numeric_var)
    s_Y <- sqrt(1/(n)*sum((ranky - rankmean)^2))
    for (i in a_i){
      b[i] <- a_i[i]*sum(ranky[order(ranky)[(1+sumindex[i]):sumindex[i+1]]] - rankmean)
    }
    lambda <- 1/(n*s_Y)*sum(b)
    samplecorr <- spearman(numeric_var, factor_var)
    if (abs(samplecorr) > maxcor) 
      samplecorr <- sign(samplecorr) * maxcor
    
    lord_estimator <- -1*samplecorr*s_X/lambda
  }
  return(lord_estimator)
}

