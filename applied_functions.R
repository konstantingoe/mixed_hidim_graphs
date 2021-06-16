
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
        I'm checking your input and declare factors for level(x)<20")
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
  
  
  if (!requireNamespace("corpcor", quietly=TRUE))
    stop("Please install package \"corpcor\".")
  if (!requireNamespace("Matrix", quietly=TRUE))
    stop("Please install package \"Matrix\".")
  
  ### eigenvalue decomposition and truncate eigenvalues at 0 

  if (!is.positive.definite(rho)) {
    initial_mat_singular <- T
    rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
  } else {
    initial_mat_singular <- F
    rho_pd <- rho
  }
  #diag(rho_pd) <- 1
  
  
  if (!requireNamespace("huge", quietly=TRUE))
    stop("Please install package \"huge\".")
  
  #now with rho_pd we have the sample correlation matrix 
  huge.result <- huge(rho_pd,nlambda=nlam,method="glasso",verbose=FALSE)
  if (thresholding == T) {
    Omega_hat <- omega.select(x=huge.result, n=n, param = param)
  } else if (thresholding == F){
    if (!requireNamespace("glasso", quietly=TRUE))
      stop("Please install package \"glasso\".")
    Omega_hat <- omega.select.drton(x=huge.result, n=n, s = rho_pd, param = param)
  }
  number_edges <- edgenumber(Omega_hat)
  max_degree <- max(sapply(1:d, function(k) (sum(abs(Omega_hat[k,]) > 0) -1))) 
  
  adj_estimate <- abs(Omega_hat) > 0
  
  output <- list("Estimated Precision Matrix" = Omega_hat, "Adjecency Matrix" = adj_estimate,
                 "Sample Correlation Matrix" = rho_pd, "Edgenumber" = number_edges, "Max Degree" = max_degree, "initial_mat_singular" = initial_mat_singular)  
  return(output)
}


### simply calculates the number of edges in a graph
edgenumber <- function(Precision=Precision, cut=0){
  sum((abs(Precision) > (0 + cut))[lower.tri((abs(Precision) > (0 + cut)))])
}


#function for Omega selection
### potentially set param = 1 (.5)
omega.select <- function(x=x, param = param, n=n){
  stopifnot((class(x)=="huge"))
  d=dim(x$data)[1]
  nlam <- length(x$lambda)
  cut <- seq(from=0, to=.1, by = .001)
  cutwhich <- rep(0, length(cut))
  for (c in 1:length(cut)){
    eBIC <- rep(0,nlam)
    for (ind in 1:nlam) {
      At <- x$icov[[ind]]
      edge <- edgenumber(At,cut=cut[c])
      eBIC[ind] <- -n*x$loglik[ind] + edge*log(n) + 4*edge*param*log(d)
    }
    cutwhich[c] <- which.min(eBIC)
    cutmaxdiff <- which.max(diff(cutwhich))
  }
  At <- x$icov[[cutwhich[cutmaxdiff+1]]]
  #diag(At) <- 1
  At.final <- -1*cov2cor(At)
  diag(At.final) <- -1*diag(At.final)
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
  Omega_hat.standardized <- -1*cov2cor(Omega_hat)
  diag(Omega_hat.standardized) <- -1*diag(Omega_hat.standardized)
  return(Omega_hat.standardized)
}



