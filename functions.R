
#### applied function for undirected graph learning
# data:             can be supplied in any form, and discrete variables will be declared as factors internally
# verbose:          supplies a notice that no factors are declared when loading the data in
# nlam:             is the vector length of panalty parameters for the glasso
# thresholding:     is an boolean variable. If False, then the EBIC of Foygel, Drton (2010) is 
#                   used (recommended for high dimensional problems) if True then thresholding will be applied (presumably better for medium and low dimensional problems)
# param:            value of additional penalty term in eBIC, default .1
# required packages: stats, polycor, glasso, huge, corpcor, matrix 

mixed.undir.graph <- function(data = data, verbose = T, nlam = 50, thresholding = F, param = .1){
  if (sum(sapply(data, is.factor)) == 0 & verbose == T){
    cat("Warning, there are no factors in the input data.
        I'm checking you input and declare factors for level(x)<20")
  }
  d <- ncol(data)
  n <- nrow(data)
  
  data <- as.data.frame(data)
  factor_ids <- sapply(data, function(id) length(unique(id)) < 20)
  data[,factor_ids] <- lapply(data[,factor_ids], factor)
  
  rho <- matrix(1,d,d)
  
  if (!requireNamespace("polycor", quietly=TRUE))
    stop("Please install package \"polycor\".")
  
  for(i in 1:(d-1)) {
    for(j in (i+1):d){
      if (is.numeric(data[,i]) & is.numeric(data[,j])){
        rho[i,j] <- rho[j,i] <- cor(data[,i], data[,j], method = "pearson")
      }
      if ((is.factor(data[,i]) & is.numeric(data[,j])) |  (is.numeric(data[,i]) & is.factor(data[,j]))) {
        if (is.factor(data[,j])) {
          rho[i,j] <- rho[j,i] <- polycor::polyserial(data[,i], data[,j])
        } else {
          rho[i,j] <- rho[j,i] <- polycor::polyserial(data[,j], data[,i])
        }
      }
      if (is.factor(data[,i]) & is.factor(data[,j])) {
        rho[i,j] <- rho[j,i] <- polycor::polychor(data[,i], data[,j])
      }
    }
  }    
  
  if (!requireNamespace("stringr", quietly=TRUE))
    stop("Please install package \"stringr\".")
  pair <- rho[lower.tri(rho)]
  if(any(abs(pair) > .9)) 
    sapply(seq_along(rho[lower.tri(rho)][which(abs(pair) > .95)]), function (k) warning(paste0('Correlation of the pair ', 
           str_c(as.character(which(rho[lower.tri(rho)][which(abs(pair) > .9)][k] == rho,
           arr.ind = T)[,1]), collapse = ",")),
           ' is close to boundary. Inverse might be misleading. '))
  
  if (!is.positive.definite(rho)) {
    rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
  } else {
    rho_pd <- rho
  }
  #diag(rho_pd) <- 1
  
  if (!requireNamespace("huge", quietly=TRUE))
    stop("Please install package \"huge\".")
  #now with rho_pd we have the sample correlation matrix 
  huge.result <- huge(rho_pd,nlambda=nlam,method="glasso",verbose=FALSE)
  if (thresholding == T) {
    Omega_hat <- omega.select(x=huge.result, n=n)
  } else if (thresholding == F){
    if (!requireNamespace("glasso", quietly=TRUE))
      stop("Please install package \"glasso\".")
    Omega_hat <- omega.select.drton(x=huge.result, n=n, s = rho_pd)
  }
  number_edges <- edgenumber(Omega_hat)
  max_degree <- max(sapply(1:d, function(k) (sum(abs(Omega_hat[k,]) > 0) -1))) 
  
  adj_estimate <- abs(Omega_hat) > 0
  
  output <- list("Estimated Precision Matrix" = Omega_hat, "Adjecency Matrix" = adj_estimate,
                 "Sample Correlation Matrix" = rho_pd, "Edgenumber" = number_edges, "Max Degree" = max_degree)  
  return(output)
}


#function for euclidean norm
euclid_norm <- function(x) sqrt(sum(x^2))

# function for number of edges

edgenumber <- function(Precision=Precision, cut=0){
  sum((abs(Precision) > (0 + cut))[lower.tri((abs(Precision) > (0 + cut)))])
}

  
#function for Omega selection
### potentially set param = 1 (.5)
omega.select <- function(x=x, param = param, n=n){
  stopifnot((class(x)=="huge"))
  d=dim(x$data)[1]
  nlambda <- length(x$lambda)
  cut <- seq(from=0, to=.1, by = .001)
  cutwhich <- rep(0, length(cut))
  for (c in 1:length(cut)){
    eBIC <- rep(0,nlambda)
    for (ind in 1:nlambda) {
      At <- x$icov[[ind]]
      edge <- edgenumber(At,cut=cut[c])
      eBIC[ind] <- -n*x$loglik[ind] + edge*log(n) + 4*edge*param*log(d)
    }
    cutwhich[c] <- which.min(eBIC)
    cutmaxdiff <- which.max(diff(cutwhich))
  }
  At <- x$icov[[cutwhich[cutmaxdiff+1]]]
  #diag(At) <- 1
  At.final <- cov2cor(At)
  
  return(At.final) 
}


omega.select.drton <- function(x=x, param = param, n=n, s = s){
  stopifnot((class(x)=="huge"))
  d=dim(x$data)[1]
  nlambda <- length(x$lambda)
  eBIC <- rep(0,nlambda)
  for (ind in 1:nlambda) {
    huge_path <- x$path[[ind]]
    edge <- edgenumber(huge_path)
    huge_path[upper.tri(huge_path, diag = T)] <- 1
    zero_mat <- which(huge_path == 0, arr.ind = T)
    loglik <- suppressWarnings(glasso::glasso(s = s, rho = 0, nobs = n, zero = zero_mat)$loglik)
    eBIC[ind] <- -2*loglik + edge*log(n) + 4* edge * param * log(d)
  }  
  
  Omega_hat <- x$icov[[which.min(eBIC)]]
  Omega_hat.standardized <- cov2cor(Omega_hat)
  
  return(Omega_hat.standardized)
}

#### write data generating function ####
#### inputs:
#1. signal strength t (make sure to choose s.t. Omega is p.d.)
#2. sample size n
#3. number of vertices (variables) d
#4. sparcity parameter c: controls number of edges i.e. non-zero entries in Omega:
#                         for this simulation use choose.c function
#                         potentially look rather specify number of edges!

generate.data <- function(mode = mode, t=.15, n = 200, d = 50, sparsity = sparsity){
  
    if (mode == "fan"){
      if (d==50){
        c <- -.2798173
      } else if (d==250) {
        c <- -0.03072694
      } else if (d==750) {
        c <- -0.01697741
      } else if (d==1500) {
        c <- -0.01227969
      } 
    #c <- choose.c(n_E=n_E, d=d)
    ##### Initiate the precision matrix \Omega = \Sigma^{-1}
    Omega <- diag(nrow = d,ncol = d)
    
    for(i in 1:(d-1)) {
      for(j in (i+1):d){
        Omega[i,j] <- Omega[j,i] <- t*rbinom(1,1,prob = (1/sqrt(2*pi))*exp(euclid_norm((runif(2, min = 0, max = 1) - runif(2, min = 0, max = 1)))/(2*c)))
        #Omega[i,j] <- Omega[j,i] <- runif(1,min = t, max = 1.3*t)*rbinom(1,1,prob = (1/sqrt(2*pi))*exp(euclid_norm((runif(2, min = 0, max = 1) - runif(2, min = 0, max = 1)))/(2*c)))
      }
    }  
    diag(Omega) <- 1
    #### check number of edges:
  
  if (!is.positive.definite(Omega)) {
    Omega <- as.matrix(nearPD(Omega, corr = T, keepDiag = T)$mat)
  } 
  #diag(Omega) <- 1
  
  } else if (mode == "er"){
    if (!requireNamespace("genscore", quietly=TRUE))
      stop("Please install package \"genscore\".")
    Omega <- cov_cons(mode="er", p=d, seed=NULL, spars=sparsity, eig=0.1)
  }
  edge_number <- edgenumber(Precision = Omega)
  
  #### retrieve Sigma ####
  ### potential alternative here: 
  Sigma <- chol2inv(chol(Omega)) #is simply a lot faster!175
  #Sigma <- solve(Omega)
  Sigma_corr <- cov2cor(Sigma)
  #### Done: Now we can simulate from multivariate normal ####
  data_0 <- mvrnorm(n=n,mu=rep(0,d), Sigma = Sigma_corr)
  return(list("data" = data_0,
              "Omega" = Omega,
              "n_E" = edge_number))
}


### also write function here:
#input: continuous data from which ordinal data will be generated 
# proportion of data that should be converted to ordinal
# number of levels for all ordinal variables (could also think of a vector of size d/2 here)
# allow for a mix of 20-leveled discrete variables (count variables) and 3-leveled ones

make.ordinal <- function(data = data, proportion = .5, n_O = 3, countvar = F, p_count = proportion*1/3){
  d <-  ncol(data)
  d_1 <- d*proportion
  ordinal <- data[,1:d_1]
  
  
  d_12 <- floor(d*p_count)
  if (countvar == F) {
    d_12 <- 0
  }
  d_11 <- d_1 - d_12
  #### generate thresholds via quantiles
  gamma <- list()
  ### reconstruct the binary case in Fan et.al. (2017)
  
  for (col in 1:d_11) { 
    if (n_O == 3){
      gamma[[col]] <- c(-Inf,
                        runif(1,quantile(ordinal[,col])[2],quantile(ordinal[,col])[3]), 
                        runif(1,quantile(ordinal[,col])[3],quantile(ordinal[,col])[4]),
                        Inf) # ordinal case
    } else if (n_O == 2) {
      gamma[[col]] <- c(-Inf,
                        runif(1,-1,1), 
                        Inf) # ordinal case
    } else {
      stop("need to define automatic quantile generation here")
    }
    ordinal[,col] <- cut(ordinal[,col], breaks = gamma[[col]], labels = F, right=T, include.lowest=T)
  }
  
  if (countvar == T) {
    for (col in ((d_11+1):d_1)) {
      gamma[[col]] <- c(-Inf,
                        seq(from = quantile(ordinal[,col])[2]*2, to = quantile(ordinal[,col])[4]*2, length.out=19),
                        Inf)
       
      ordinal[,col] <- cut(ordinal[,col], breaks = gamma[[col]], labels = F, right=T, include.lowest=T)
    }
  }
  data_mixed <- as.data.frame(cbind(ordinal,data[,(d_1+1):d]))
  
  ### declare ordinal variables as factors!
  for (f in 1:d_1) {
    data_mixed[,f] <- factor(data_mixed[,f],ordered = T) 
  }
  return(data_mixed)
}

### write polychoric function ####
### input: data frame with mixed variables where ordinal variables are denoted as factors
### output: positive definite sample correlation matrix \hat{Sigma}
mixed.omega <- function(data=data, verbose = T){
    if (sum(sapply(data, is.factor)) == 0 & verbose == T){
      cat("Warning, there are no factors in the input data.
          Did you declare ordinal variables as factors?")
    }
    d <- ncol(data)
    rho <- matrix(1,d,d)
    
    ### retrieve maxlevels for polychoric:
    maxlevel <- max(sapply(1:ncol(data), function(i) length(levels(data[,i])))) +1
    
    for(i in 1:(d-1)) {
      for(j in (i+1):d){
        if (is.numeric(data[,i]) & is.numeric(data[,j])){
          rho[i,j] <- rho[j,i] <- cor(data[,i], data[,j], method = "pearson")
        }
        if ((is.factor(data[,i]) & is.numeric(data[,j])) |  (is.numeric(data[,i]) & is.factor(data[,j]))) {
          if (is.factor(data[,j])) {
            rho[i,j] <- rho[j,i] <- polycor::polyserial(data[,i], data[,j])
          } else {
            rho[i,j] <- rho[j,i] <- polycor::polyserial(data[,j], data[,i])
          }
        }
        if (is.factor(data[,i]) & is.factor(data[,j])) {
          rho[i,j] <- rho[j,i] <- polycor::polychor(data[,i], data[,j])
          #rho[i,j] <- rho[j,i] <- suppressMessages(psych::polychoric(data[,c(i,j)], max.cat = maxlevel, progress = F)[[1]][1,2])
        }
      }
    }    
    
    if (!requireNamespace("stringr", quietly=TRUE))
      stop("Please install package \"stringr\".")
    pair <- rho[lower.tri(rho)]
    if(any(abs(pair) > .9)) 
      sapply(seq_along(rho[lower.tri(rho)][which(abs(pair) > .95)]), function (k) warning(paste0('Correlation of the pair ', 
             str_c(as.character(which(rho[lower.tri(rho)][which(abs(pair) > .9)][k] == rho,
             arr.ind = T)[,1]), collapse = ",")),
             ' is close to boundary. Inverse might be misleading. '))    
    if (!is.positive.definite(rho)) {
      rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
    } else {
      rho_pd <- rho
    }
    #diag(rho_pd) <- 1
  return(rho_pd)
}

#### Mixed Tau for Fan et.al. Method:

mixed.omega.kendall <- function(data = data, verbose = T){
  
  if (sum(sapply(data, is.factor)) == 0 & verbose == T){
    cat("Warning, there are no factors in the input data.
        Did you declare ordinal variables as factors?")
  }
  d <- ncol(data)
  n <- nrow(data)
  cor1 <- hatR <- tau <- na <- nb <- nc <- nd <- matrix(1,d,d)
  
  for(i in 1:(d-1)) {
    for(j in (i+1):d){
      if (is.numeric(data[,i]) & is.numeric(data[,j])){
        ### Fan et.al.
        cor1[j,i] <- cor1[i,j] <- cor(data[,i],data[,j],method="kendall")
        hatR[i,j] <- hatR[j,i] <- sin(cor1[i,j]*pi/2)
        ###
      }
      if ((is.factor(data[,i]) & is.numeric(data[,j])) |  (is.numeric(data[,i]) & is.factor(data[,j]))) {
        if (is.factor(data[,j])) {
          ### Fan et.al.
          cor1[j,i] <- cor1[i,j] <- cor(data[,i], as.numeric(data[,j]), method="kendall")
          hatR[i,j] <- hatR[j,i] <- 2^(1/2)*sin(cor1[i,j]*pi/2)
          ###
        } else {
          ### Fan et.al.
          cor1[j,i] <- cor1[i,j] <- cor(data[,j],as.numeric(data[,i]),method="kendall")
          hatR[i,j] <- hatR[j,i] <- 2^(1/2)*sin(cor1[i,j]*pi/2)
          ###
        }
      }
      if (is.factor(data[,i]) & is.factor(data[,j])) {
        ### Fan et.al. ### population kendall's tau
        na[i,j] <- na[i,j] <- sum((data[,i] == 1)*(data[,j] == 1))
        nb[i,j] <- nb[i,j] <- sum((data[,i] == 1)*(data[,j] == 0))
        nc[i,j] <- nc[i,j] <- sum((data[,i] == 0)*(data[,j] == 1))
        nd[i,j] <- nd[i,j] <- sum((data[,i] == 0)*(data[,j] == 0))
        
        tau[i,j] <- tau[j,i] <- 2*(na[i,j]*nd[i,j] - nb[i,j]*nc[i,j])/(n*(n-1))
        #tau[i,j] <- tau[j,i] <- cor(as.numeric(data[,i]),as.numeric(data[,j]),method="kendall")
        hatR[i,j] <- hatR[j,i] <- sin(pi*tau[i,j])
        ###
      }
    }
  }    
  
  if (!is.positive.definite(hatR)) {
    hatR_pd <- as.matrix(nearPD(hatR, corr = T, keepDiag = T)$mat)
  } else {
    hatR_pd <- hatR
  }
  return(hatR_pd)
}


#### functions for tpr and fpr

tpr <- function(truth=truth, estimate=estimate){
  n_E <- edgenumber(truth)
  sum((truth != 0 & estimate != 0)[lower.tri((truth != 0 & estimate != 0))]) / n_E
}

fpr <- function(truth=truth, estimate=estimate){
  n_E <- edgenumber(truth)
  d = dim(truth)[1]
  
  sum((truth == 0 & estimate != 0)[lower.tri((truth == 0 & estimate != 0))]) / (d*(d-1)/2 - n_E) 
}


#### input truth: True Precision Matrix
####       huge_obj: huge object e.g. glasso output
AUC <- function(truth = truth, huge_obj = huge_obj){
  stopifnot((class(huge_obj)=="huge"))
  #d = dim(truth)[1]
  #n_E <- sum(truth[lower.tri(truth)] != 0)
  ### here we need adjacency matrices!
  adj_estimate <- lapply(1:nlam, function(q) abs(huge_obj$icov[[q]]) > 0)

  fpr_lambda <- sapply(1:nlam, function(l) fpr(truth = truth,estimate = adj_estimate[[l]]))
  tpr_lambda <- sapply(1:nlam, function(l) tpr(truth = truth,estimate = adj_estimate[[l]]))
  
  fp <- c(fpr_lambda,1)
  tp <- c(tpr_lambda,1)
  dfp <- c(diff(fp),0)
  dtp <- c(diff(tp),0)
  
  AUC <- sum(tp * dfp) + sum(dtp * dfp)/2
  return(AUC)
}



#### need new big function because result object is getting too large
### Sigma: Sample correlation matrix
### Omega: True Precision matrix
### nlam: positive integer indicating number of tuning parameters
### n: sample size


glasso.results <- function(Sigma = Sigma, Omega = Omega, nlam = nlam, n=n, matexport = matexport, thresholding = F, param = param){
  ### this takes a while
  huge.result <- huge(Sigma,nlambda=nlam,method="glasso",verbose=FALSE)
  if (thresholding == T) {
    Omega_hat <- omega.select(x=huge.result, n=n, param = param)
  } else if (thresholding == F){
    Omega_hat <- omega.select.drton(x=huge.result, n=n, s = Sigma, param = param)
  }
  
  frobenius <- base::norm(Omega_hat - Omega, type = "F")
  AUC_hat <- AUC(truth = Omega, huge_obj = huge.result)
  
  adj_estimate <- abs(Omega_hat) > 0
  
  tpr <- tpr(truth = Omega, estimate = adj_estimate)
  fpr <- fpr(truth = Omega, estimate = adj_estimate)
  
  huge.result <- adj_estimate <-  NULL
  
  if (matexport == T){
    output <- list("True Precision Matrix" = Omega, "Estimated Precision Matrix" = Omega_hat, "Frobenius norm" = frobenius,
                 "FPR" = fpr, "TPR" = tpr, "AUC" = AUC_hat)
  } else if (matexport == F) {
    Omega_hat <- Omega <-  NULL
    output <- list("Frobenius norm" = frobenius,
                   "FPR" = fpr, "TPR" = tpr, "AUC" = AUC_hat)
  }
  return(output)
}



#### boil it down even further... we don't really need the data object.

serverrun <- function(mode = "fan", n=n, d=d, sparsity = .1, nlam=nlam, param = .1, matexport = F, 
                      namevector = c("binary" = T, "ordinal" = T, "poisson" = T), 
                      unbalanced = unbalanced, low = .05, high = .1, latent = latent){
  
  data <- generate.data(t=t, n = n, d = d, mode = mode)
  data_0 <- data[[1]]
  Omega <- data[[2]]
  
  data <- NULL
  
  if (latent == T) {
    rho_latent <- mixed.omega(data_0)
    results_latent <- glasso.results(Sigma=rho_latent, Omega=Omega, nlam=nlam, n=n, matexport = matexport, param = param)
  }
  data_mixed <- make.ordinal.general(data_0, namevector = namevector, unbalanced = unbalanced, low = low, high = high)
  data_0 <- NULL
  ### learn sample correlation matrix   
  rho <- mixed.omega(data_mixed)
  data_mixed <- NULL
  
  results <- glasso.results(Sigma=rho, Omega=Omega, nlam=nlam, n=n, matexport = matexport,  param = param)
  
  if (latent == T) {
    results <- list("latent"=results_latent, "mixed" = results)
  }
  return(results)
}

serverrun.kendall <- function(t=.15, n = n, d = d, latent = F, nlam=nlam, matexport = F, countvar = T, mode = mode){
  
  data <- generate.data(t=t, n = n, d = d, mode = mode)
  data_0 <- data[[1]]
  Omega <- data[[2]]
  
  data <- NULL
  
  data_mixed <- make.ordinal(data_0, countvar = countvar, n_O = 2)
  data_0 <- NULL
  ### learn sample correlation matrix   
  rho_kendall <- mixed.omega.kendall(data_mixed)
  rho_poly <- mixed.omega(data_mixed)
  
  data_mixed <- NULL
  
  results_kendall <- glasso.results(Sigma=rho_kendall, Omega=Omega, nlam=nlam, n=n, matexport = matexport)
  results_poly <- glasso.results(Sigma=rho_poly, Omega=Omega, nlam=nlam, n=n, matexport = matexport)
  results <- list("kendall"=results_kendall, "poly" = results_poly)
  return(results)
}



#### extract results

extract.result <- function(results=results, which = c("F", "TPR", "FPR", "AUC")){
  sim <- length(results)
  if (which == "F"){
    extract <- sapply(1:sim, function(k) results[[k]][["Frobenius norm"]])
  } else if (which == "TPR") {
    extract <- sapply(1:sim, function(k) results[[k]][["TPR"]])
  } else if (which == "FPR") {
    extract <- sapply(1:sim, function(k) results[[k]][["FPR"]])  
  } else if (which == "AUC") {
    extract <- sapply(1:sim, function(k) results[[k]][["AUC"]])
  }
  return(extract)
}


mixed.omega_fast <- function(data=data, verbose = F){
  if (sum(sapply(data, is.factor)) == 0 & verbose == T){
    cat("Warning, there are no factors in the input data.
          Did you declare ordinal variables as factors?")
  }
  d <- ncol(data)
  rho <- matrix(1,d,d)
  
  for (j in 1:(d-1)){
    rho[(j+1):d,j] <- sapply(j:(d-1), function(i) 
      mixed_corr(data[,j], data[,i+1]))
  }
  
  upperTriangle(rho) <- lowerTriangle(rho, byrow=TRUE)
  
  
  if (!is.positive.definite(rho)) {
    rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
  } else {
    rho_pd <- rho
  }
  #diag(rho_pd) <- 1
  return(rho_pd)
}



mixed_corr <- function(vec1, vec2){
  
  if (is.numeric(vec1) & is.numeric(vec2)){
    corr <- cor(vec1,vec2, method = "pearson") 
  }  
  if ((is.factor(vec1) & is.numeric(vec2)) |  (is.numeric(vec1) & is.factor(vec2))) {
    if (is.factor(vec2)) {
      corr <- polyserial(vec1, vec2)
    } else {
      corr <- polyserial(vec2, vec1)
    }
  }  
  if (is.factor(vec1) & is.factor(vec2)) {
    corr <- polychor(vec1, vec2)
  }
  return(corr)
}


kendalls.results <- function(result_object = result_object){
  table <- rbind(
    c("kendall" = mean(sapply(1:sim, function(k) result_object[[k]][[1]][["Frobenius norm"]])),"kendall_sd" = sd(sapply(1:sim, function(k) result_object[[k]][[1]][["Frobenius norm"]])),
      "poly" = mean(sapply(1:sim, function(k) result_object[[k]][[2]][["Frobenius norm"]])),"poly_sd" = sd(sapply(1:sim, function(k) result_object[[k]][[2]][["Frobenius norm"]]))),
    
    c(mean(sapply(1:sim, function(k) result_object[[k]][[1]][["FPR"]])),sd(sapply(1:sim, function(k) result_object[[k]][[1]][["FPR"]])),
      mean(sapply(1:sim, function(k) result_object[[k]][[2]][["FPR"]])),sd(sapply(1:sim, function(k) result_object[[k]][[2]][["FPR"]]))),
    
    c(mean(sapply(1:sim, function(k) result_object[[k]][[1]][["TPR"]])),sd(sapply(1:sim, function(k) result_object[[k]][[1]][["TPR"]])),
      mean(sapply(1:sim, function(k) result_object[[k]][[2]][["TPR"]])),sd(sapply(1:sim, function(k) result_object[[k]][[2]][["TPR"]]))),
    
    c(mean(sapply(1:sim, function(k) result_object[[k]][[1]][["AUC"]])),sd(sapply(1:sim, function(k) result_object[[k]][[1]][["AUC"]])),
      mean(sapply(1:sim, function(k) result_object[[k]][[2]][["AUC"]])),sd(sapply(1:sim, function(k) result_object[[k]][[2]][["AUC"]])))
  )
  return(table)
}

make.ordinal.general <- function(data = data, proportion = .5, namevector = c("binary" = T, "ordinal" = T, "poisson" = T), unbalanced = .2, low = .05, high =.1, lambda = 6, num_breaks = round(runif(1,3,10))){
  d <-  ncol(data)
  n <- nrow(data)
  d_1 <- floor(d*proportion)
  ordinal <- data[,1:d_1]
  
  #### now we can do with ordinal what we want 
  
  if (namevector[1] == T & namevector[2] == T & namevector[3] == T) {
    p_devide <- diff(c(0,floor(d_1/sum(namevector)), d_1))
    ordinal_binary <- ordinal[,1:p_devide[1]]
    ordinal_ordinal <- ordinal[,(p_devide[1]+1):p_devide[2]]
    ordinal_poisson <- ordinal[,(p_devide[2]+1):d_1]
  }    
  if (namevector[1] == T & namevector[2] == T & namevector[3] == F) {
    p_devide <- diff(c(0,floor(d_1/sum(namevector))))
    ordinal_binary <- ordinal[,1:p_devide[1]]
    ordinal_ordinal <- ordinal[,(p_devide[1]+1):d_1]
    ordinal_poisson <- NULL
  }  
  if (namevector[1] == T & namevector[2] == F & namevector[3] == T) {
    p_devide <- diff(c(0,floor(d_1/sum(namevector))))
    ordinal_binary <- ordinal[,1:p_devide[1]]
    ordinal_ordinal <- NULL
    ordinal_poisson <- ordinal[,(p_devide[1]+1):d_1]
  }  
  if (namevector[1] == F & namevector[2] == T & namevector[3] == T) {
    p_devide <- diff(c(0,floor(d_1/sum(namevector))))
    ordinal_binary <- NULL
    ordinal_ordinal <- ordinal[,1:p_devide[1]]
    ordinal_poisson <- ordinal[,(p_devide[1]+1):d_1]
  }
  if (namevector[1] == F & namevector[2] == F & namevector[3] == T) {
    ordinal_binary <- NULL
    ordinal_ordinal <- NULL
    ordinal_poisson <- ordinal[,1:d_1]
  }
  if (namevector[1] == F & namevector[2] == T & namevector[3] == F) {
    ordinal_binary <- NULL
    ordinal_ordinal <- ordinal[,1:d_1]
    ordinal_poisson <- NULL
  }
  if (namevector[1] == T & namevector[2] == F & namevector[3] == F) {
    ordinal_binary <- ordinal[,1:d_1]
    ordinal_ordinal <- NULL
    ordinal_poisson <- NULL
  }
  
  if (namevector["binary"] == T){
    #### Xbinary
    #### split binary so as to control fraction of unbalanced binary data
    pbin1 <- runif(floor(ncol(ordinal_binary)*(1-unbalanced)),.4,.6)
    pbin2 <- runif((ncol(ordinal_binary) - floor(ncol(ordinal_binary)*(1-unbalanced))),low,high)
    pbin <- c(pbin1, pbin2)
    for(i in 1:ncol(ordinal_binary)){
      ordinal_binary[,i] <- qbinom(pnorm(scale(ordinal_binary[,i])),size=1,prob = pbin[i])
    }
  }
  if (namevector["ordinal"] == T){
    cum.mat <- list()
    for (k in 1:ncol(ordinal_ordinal)){
      breaks <- runif(num_breaks)
      sum_breaks <- sum(breaks)
      breaks_norm <- sort(breaks)/sum_breaks
      cum.mat[[k]]<- c(0,cumsum(breaks_norm))
    }
    for(i in 1:ncol(ordinal_ordinal)){
      u <- pnorm(scale(ordinal_ordinal[,i]))
      ordinal_ordinal[,i] <- cut(u, breaks = cum.mat[[i]], include.lowest = T, ordered_result = T, labels = 1:(length(cum.mat[[i]])-1))
    }
  }
  if (namevector["poisson"] == T){
    #### Poisson Variables from Threshold generation 
    lambda_sample <- abs(rnorm(ncol(ordinal_poisson),mean = lambda, sd = 2))
    for (i in 1:ncol(ordinal_poisson)) {
      ordinal_poisson[,i] <- qpois(pnorm(scale(ordinal_poisson[,i])),lambda = lambda_sample[i])
    }
  } 
  
  ordinal <- cbind(ordinal_binary, ordinal_ordinal, ordinal_poisson)
  data_mixed <- as.data.frame(cbind(ordinal,data[,(d_1+1):d]))
  ### declare ordinal variables as factors!
  for (f in 1:d_1) {
    data_mixed[,f] <- factor(data_mixed[,f],ordered = T) 
  }
  return(data_mixed)
}

### what we want is basically frobenius norm as a function of n and success probability:
unbalanced.run <- function(mode = "fan", n=n, d=d, sparsity = .1, nlam=nlam, matexport = F, namevector = c("binary" = T, "ordinal" = F, "poisson" = F), 
                           unbalanced = 0, low = .05, high = .1){
  data <- generate.data(mode = mode, t=t, n = n, d = d, sparsity = sparsity)
  data_0 <- data[[1]]
  Omega <- data[[2]]
  data <- NULL
  ### here at this point generate a range of make.ordinals!
  data_mixed <- make.ordinal.general(data_0, namevector = namevector, unbalanced = unbalanced, low = low, high = high)
  data_0 <- NULL
  rho_poly <- mixed.omega(data_mixed)
  data_mixed <- NULL
  results_poly <- glasso.results(Sigma=rho_poly, Omega=Omega, nlam=nlam, n=n, matexport = matexport)

  return(results_poly)
}


results_generator <- function(results_object){
    table <- round(as_tibble(rbind(
              c("latent" = mean(sapply(1:sim, function(k) results_object[[k]]$latent[["Frobenius norm"]])),
                "sd_latent" = sd(sapply(1:sim, function(k) results_object[[k]]$latent[["Frobenius norm"]])),
                "mixed" = mean(sapply(1:sim, function(k) results_object[[k]]$mixed[["Frobenius norm"]])),
                "sd_mixed" = sd(sapply(1:sim, function(k) results_object[[k]]$mixed[["Frobenius norm"]]))),
              
              c(mean(sapply(1:sim, function(k) results_object[[k]]$latent[["FPR"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$latent[["FPR"]])),
                mean(sapply(1:sim, function(k) results_object[[k]]$mixed[["FPR"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$mixed[["FPR"]]))),
              
              c(mean(sapply(1:sim, function(k) results_object[[k]]$latent[["TPR"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$latent[["TPR"]])),
                mean(sapply(1:sim, function(k) results_object[[k]]$mixed[["TPR"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$mixed[["TPR"]]))),
              
              c(mean(sapply(1:sim, function(k) results_object[[k]]$latent[["AUC"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$latent[["AUC"]])),
                mean(sapply(1:sim, function(k) results_object[[k]]$mixed[["AUC"]])),
                sd(sapply(1:sim, function(k) results_object[[k]]$mixed[["AUC"]])))
            )),4)
  return(table)
}
