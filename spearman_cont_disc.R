#### simulating spearman's rho ####

rm(list = ls())
source("packages.R")
source("functions.R")
#
set.seed(1234)
sim <- 200
n <- 100000
grid <- seq(-.9,.9,.1)
breaks <- runif(5)
sum_breaks <- sum(breaks)
breaks_norm <- sort(breaks)/sum_breaks
cum.mat <- c(0,cumsum(breaks_norm))
run <- list()

Sigma <- lapply(seq_along(grid), function(i) matrix(c(1,grid[i], grid[i], 1), nrow = 2, ncol = 2, byrow = T))

for (k in 1:sim){
  Z_X <- lapply(seq_along(grid), function(i) mvrnorm(n=n, mu = c(0,0), Sigma = Sigma[[i]]))
  X_x <- lapply(seq_along(grid), function(i) (Z_X[[i]]))
  continuous <- lapply(seq_along(grid), function(i) (X_x[[i]][,2]))
  ordinal <- lapply(seq_along(grid), function(i) 
              cut(pnorm(scale(X_x[[i]][,1])), breaks = cum.mat, 
                  include.lowest = T, ordered_result = T, labels = 1:(length(cum.mat)-1)))
  
  run[[k]] <- sapply(seq_along(grid), function(i) 
    adhoc_lord_sim(continuous[[i]], ordinal[[i]])-grid[i])
}

### now for plotting

dataset <- sapply(1:sim, function(k) run[[k]][1])

for (i in 2:length(grid)){
  dataset <- cbind(dataset, sapply(1:sim, function(k) run[[k]][i]))
}

dataset <- as.data.frame(dataset)
colnames(dataset) <- grid
library(reshape) 
plotdata <- melt(dataset)

linear <- ggplot(plotdata, aes(x=variable, y=value)) + 
          geom_boxplot() +
          geom_hline(yintercept=0, linetype='dotted', col = 'red')+
          xlab(expression(rho["latent"])) + 
          ylab(expression(hat(rho)["polyserial"]-rho["latent"]))  
linear  

breaks <- runif(9)
sum_breaks <- sum(breaks)
breaks_norm <- sort(breaks)/sum_breaks
cum.mat <- c(0,cumsum(breaks_norm))

#vec <- vector(length = 100)
#for (i in 1:100){
a = 0.15; b = 0.5; c =1.5
library(MASS)
#set.seed(2021)
Z = (mvrnorm(n = 1E6, rep(0, 2), matrix(c(1,a,a,1),2,2)))
z_j = Z[,1]
z_k = Z[,2]
funs <- function (t){
  if (t <= b) {res=0}
  else {res=1}
  return (res)
}
u <- 
x_k <- cut(pnorm(scale(z_k)), breaks = cum.mat, include.lowest = T, ordered_result = T, labels = 1:(length(cum.mat)-1))

x_k =  #sapply(z_k, funs)
vec[i] <- adhoc_lord_sim(z_j,x_k)


funs <- function (t){
  if (t <= c) {res=1}
  else {res=2}
  return (res)
}
x_j = sapply(z_j, funs)

cor(z_j,z_k)

polyserial(z_j,x_k)
fan.case.2(z_j,x_k)

adhoc_lord(z_j,x_k)
adhoc_lord_sim(z_j,x_k)

