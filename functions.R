
#function for euclidean norm
euclid_norm <- function(x) sqrt(sum(x^2))

#function for Omega selection
omega.select <- function(x=x, param = 0.1, n=n){
  stopifnot((class(x)=="huge"))
  d=dim(x$data)[1]
  eBIC <- rep(0,nlam)
  for (ind in 1:nlam) {
    At <- x$icov[[ind]]
    #At=solve(cov2cor(solve(At)))
    edge <- sum((abs(At) > .05)[lower.tri((abs(At) > .05))])
    eBIC[ind] <- -n*x$loglik[ind] + edge*log(n) + 4*edge*param*log(d)
  }
  
  indlam <- which.min(eBIC)
  At <- x$icov[[indlam]]
  return(At) 
}

#### write data generating function ####
#### inputs:
#1. signal strength t (make sure to choose s.t. Omega is p.d.)
#2. sample size n
#3. number of vertices (variables) d
#4. sparcity parameter c: controls number of edges i.e. non-zero entries in Omega:
#                         for this simulation use choose.c function
#                         potentially look rather specify number of edges!


generate.data <- function(t=.15, n = 200, d = 50, n_E = 200){
  if (d==50){
    c <- -.2798173
  } else if (d==250) {
    c <- -.0285
  } else if (d==3000) {
    c <- -.002
  }
  #c <- choose.c(n_E=n_E, d=d)
  ##### Initiate the precision matrix \Omega = \Sigma^{-1}
  Omega <- diag(nrow = d,ncol = d)
  
  for(i in 1:(d-1)) {
    for(j in (i+1):d){
      Omega[i,j] <- Omega[j,i] <- t*rbinom(1,1,prob = (1/sqrt(2*pi))*exp(euclid_norm((runif(2, min = 0, max = 1) - runif(2, min = 0, max = 1)))/(2*c)))
    }
  }  
  diag(Omega) <- 1
  #### check number of edges:
  
  edgenumber <- sum(Omega[lower.tri(Omega)] != 0)
  # from possible choose(d,2)
  #### retrieve Sigma ####
  Sigma <- solve(Omega)
  Sigma_corr <- cov2cor(Sigma)
  #### Done: Now we can simulate from multivariate normal ####
  data_0 <- mvrnorm(n=n,mu=rep(0,d), Sigma = Sigma_corr)
  return(list("data" = data_0,
              "Omega" = Omega,
              "n_E" = edgenumber))
}

### choose c function to control number of edges
choose.c <- function(n_E=200, d=50){
  ### what does p_ij have to be:
  p <- (n_E) / choose(d,2)*2
  c <- 1/(4*(log(p) + .5*log(2*pi)))
  ### initiate c 
  #c <- seq(from=-1, to=0, by = .01)
  #prob = sapply(seq_along(c), function(i) (1/sqrt(2*pi))*exp(.5/(2*(c[i]))))
  #c_1 <- c[which.min(abs(prob-p))]
  return(c)
}

### also write function here:
#input: continuous data from which ordinal data will be generated 
# proportion of data that should be converted to ordinal
# number of levels for all ordinal variables (could also think of a vector of size d/2 here)

make.ordinal <- function(data = data, proportion = .5, n_O = 3){
  d = ncol(data)
  d_1 <- d*proportion
  ordinal <- data[,1:d_1]
  
  #### generate thresholds via quantiles
  gamma <- list()
  ### reconstruct the binary case in Fan et.al. (2017)
  
  for (col in 1:ncol(ordinal)) { 
    if (n_O == 3){
      gamma[[col]] <- c(-Inf, runif(1,quantile(ordinal[,col])[2],quantile(ordinal[,col])[3]), runif(1,quantile(ordinal[,col])[3],quantile(ordinal[,col])[4]), Inf) # ordinal case
    } else {
      cat("need to define automatic quantile generation here")
    }
    ordinal[,col] <- cut(ordinal[,col], breaks = gamma[[col]], labels = F, right=T, include.lowest=T)
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
mixed.omega <- function(data=data, verbose = F){
    if (sum(sapply(data, is.factor)) == 0 & verbose == T){
      cat("Warning, there are no factors in the input data.
          Did you declare ordinal variables as factors?")
    }
    d <- ncol(data)
    rho <- matrix(1,d,d)
    
    for(i in 1:(d-1)) {
      for(j in (i+1):d){
        if (is.numeric(data[,i]) & is.numeric(data[,j])){
          rho[i,j] <- rho[j,i] <- cor(data[,i], data[,j], method = "pearson")
        }
        if ((is.factor(data[,i]) & is.numeric(data[,j])) |  (is.numeric(data[,i]) & is.factor(data[,j]))) {
          if (is.factor(data[,j])) {
            rho[i,j] <- rho[j,i] <- polyserial(data[,i], data[,j])
          } else {
            rho[i,j] <- rho[j,i] <- polyserial(data[,j], data[,i])
          }
        }
        if (is.factor(data[,i]) & is.factor(data[,j])) {
          rho[i,j] <- rho[j,i] <- polychor(data[,i], data[,j])
        }
      }
    }    
    
    if (!is.positive.definite(rho)) {
      rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
    } else {
      rho_pd <- rho
    }
  return(rho_pd)
}

#### functions for tpr and fpr

tpr <- function(truth=truth, estimate=estimate){
  n_E <- sum(truth[lower.tri(truth)] != 0)
  sum((truth != 0 & estimate != 0)[lower.tri((truth != 0 & estimate != 0))]) / n_E
}

fpr <- function(truth=truth, estimate=estimate){
  n_E <- sum(truth[lower.tri(truth)] != 0)
  d = dim(truth)[1]
  
  sum((truth == 0 & estimate != 0)[lower.tri((truth == 0 & estimate != 0))]) / (d*(d-1)/2 - n_E) 
}

#### Partial area under the curve conditional on FPR --- method used by McClish (1989)####

#### input truth: True Precision Matrix
####       huge_obj: huge object e.g. glasso output
pAUC <- function(truth = truth, huge_obj = huge_obj){
  stopifnot((class(huge_obj)=="huge"))
  #d = dim(truth)[1]
  #n_E <- sum(truth[lower.tri(truth)] != 0)
  fpr_lambda <- sapply(1:nlam, function(l) fpr(truth = truth,estimate = huge_obj$icov[[l]]))
  tpr_lambda <- sapply(1:nlam, function(l) tpr(truth = truth,estimate = huge_obj$icov[[l]]))
  
  d_fpr_lambda <- c(diff(fpr_lambda),0)
  d_tpr_lambda <- c(diff(tpr_lambda),0)
  #use partial area under the curve: pAUC = .5[1 + (AUC - min(FPR))/(max(DPR) - min(FPR))]
  pAUC <- (1 + ((sum(tpr_lambda * d_fpr_lambda) + sum(d_tpr_lambda * d_fpr_lambda)/2) - min(fpr_lambda))/(max(fpr_lambda) - min(fpr_lambda)))/2
  return(pAUC)
}

