rm(list = ls())
source("packages.R")
source("functions.R")
source("applied_functions.R")

set.seed(1221)


#### try out polyserial nonparanormal

X_tilde <- mvrnorm(n=1000,mu=rep(0,2), Sigma = matrix(c(1,.4,.4,1),nrow = 2,ncol = 2))
pbin <- runif(1,.2,.3)
binarize <- qbinom(pnorm(scale(X_tilde[,2])),size=1,prob = pbin)
X_binary <- cbind(X_tilde[,1], (binarize))

X_binary <- as.data.frame(X_binary)
X_binary[,2] <- factor(X_binary[,2])

thresholds <- function(vector){
  cumprop <- as.numeric(qnorm(table(vector)[-1]/length(vector)))
}

testrun <- pointpolynormal(X_binary)

# Input data should be dataframe or matrix where one column is numeric and the other one factor
pointpolynormal <- function(data){
  data <- as.data.frame(data)
  
  if (any(is.factor(data) == F)){
    factor_id <- sapply(data, function(id) length(unique(id)) < 10)
  } else {
    factor_id <- sapply(data, is.factor)
  }
  
  ### retrieve numeric and discrete variable
  numeric_var <- data[,factor_id == F]
  factor_var <- data[,factor_id == T]
  
  
  ### calculate threshold vector and attach infinity bounds for convenience
  cumprop <- c(-Inf, thresholds(factor_var), Inf)
  ### calculate rank correlation
  samplecorr <- spearman(numeric_var, factor_var)
  
  ### calculate gaussian density evaluated at the estimated thresholds
  densityprob <- dnorm(cumprop)
  
  
  ### calculate sample variance 
  
  ### first P(Y = y_j)
  p_hat <- vector(mode = "numeric", (length(cumprop)-1))
  for (j in 2:length(cumprop)){
    p_hat[j-1] <- pnorm(cumprop[j]) - pnorm(cumprop[j-1])
  }  
  
  ### calulate sample mean
  samplmean <- sum(as.numeric(levels(as.factor(factor_var)))*p_hat)
  ### calculate sample variance
  samplevar <- sum(as.numeric(levels(as.factor(factor_var)))^2*p_hat) - samplmean^2
  
  ### calculate weighting scheme in case levels are not consecutive
  consecutive <- vector(mode = "numeric", length(head(seq_along(as.numeric(levels(as.factor(factor_var)))),-1)))
  for (j in head(seq_along(as.numeric(levels(as.factor(factor_var)))),-1)){
    consecutive[j] <- as.numeric(levels(as.factor(factor_var)))[j+1] - as.numeric(levels(as.factor(factor_var)))[j]
  } 
  
  ### calculate adhoc nonparanormal point polyserial estimator
  adhoc_nonpara <- samplecorr*sqrt(samplevar)/(sum(head(densityprob, -1)*(consecutive)))
  return(adhoc_nonpara)
}
