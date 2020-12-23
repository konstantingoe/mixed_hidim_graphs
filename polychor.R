rm(list = ls())
source("packages.R")
source("functions.R")

#### sample from multivariate normal
set.seed(1234)


sim = 100  # simulation runs
n <- c(200,200,600) # sample size
d <- c(50,250,3000) # dimensionality
d_1 = d/2
n_E <- 200 # sparsity level of the graph 
t <- .15 # signal strength

frobenius1 <- vector(length = sim)
frobenius2 <- vector(length = sim)

fpr_at1 <- vector(length = sim)
tpr_at1 <- vector(length = sim)

fpr_at2 <- vector(length = sim)
tpr_at2 <- vector(length = sim)

edgenumber <- vector(length = sim)

#### begin simulation ####
for (run in 1:sim){
  
  #### data generation via sparse Omega ####
  data <- generate.data(t=t,n=n[1],d=d[1], n_E = n_E)
  
  data_0 <- data$data # data according to Multivariate normal 
  Omega <- data$Omega # underlying undirected graph via precision matrix
  edgenumber[run] <- data$n_E
  
  ordinal <- data_0[,1:d_1[1]]
  
  #### generate thresholds:
  # potentially use quantiles of variables!
  gamma <- list()
  
  ### reconstruct the binary case in Fan et.al. (2017)
  
  for (col in 1:ncol(ordinal)) { 
    #gamma[[col]] <- c(-Inf, a <- runif(1,-1,+1), runif(1,a,+1.5), Inf) # ordinal case
    gamma[[col]] <- c(-Inf, runif(1,-1,+1), Inf) # binary case
    
    ordinal[,col] <- cut(ordinal[,col], breaks = gamma[[col]], labels = F, right=T, include.lowest=T)
  }
  
  data_mixed <- as.data.frame(cbind(ordinal,data_0[,(d_1[1]+1):d[1]]))
  
  ### declare ordinal variables as factors!
  for (k in 1:d_1[1]) {
    data_mixed[,k] <- factor(data_mixed[,k] -1,ordered = F) 
  }
  
  #### compare performance with Fan et.al. and the Copula stuff ####
  #### reconstruct Sigma
  
  ####calculate kendall tau#####
  cor1 <- tau <- hatR <- na <- nb <- nc <- nd <- matrix(1,d[1],d[1])
  rho <- matrix(1,d[1],d[1])
  nlam <- 50

  
  for(i in 1:(d[1]-1)) {
    for(j in (i+1):d[1]){
      if (is.numeric(data_mixed[,i]) & is.numeric(data_mixed[,j])){
        rho[i,j] <- rho[j,i] <- cor(data_mixed[,i], data_mixed[,j], method = "pearson")
        
        ### Fan et.al.
        cor1[j,i] <- cor1[i,j] <- cor(data_mixed[,i],data_mixed[,j],method="kendall")
        hatR[i,j] <- hatR[j,i] <- sin(cor1[i,j]*pi/2)
        ###
      }
      if ((is.factor(data_mixed[,i]) & is.numeric(data_mixed[,j])) |  (is.numeric(data_mixed[,i]) & is.factor(data_mixed[,j]))) {
        if (is.factor(data_mixed[,j])) {
        rho[i,j] <- rho[j,i] <- polyserial(data_mixed[,i], data_mixed[,j])
        
        ### Fan et.al.
        cor1[j,i] <- cor1[i,j] <- cor(data_mixed[,i], as.numeric(data_mixed[,j]), method="kendall")
        hatR[i,j] <- hatR[j,i] <- 2^(1/2)*sin(cor1[i,j]*pi/2)
        ###
        } else {
        rho[i,j] <- rho[j,i] <- polyserial(data_mixed[,j], data_mixed[,i])
        ### Fan et.al.
        cor1[j,i] <- cor1[i,j] <- cor(data_mixed[,j],as.numeric(data_mixed[,i]),method="kendall")
        hatR[i,j] <- hatR[j,i] <- 2^(1/2)*sin(cor1[i,j]*pi/2)
        ###
        }
      }
      if (is.factor(data_mixed[,i]) & is.factor(data_mixed[,j])) {
        rho[i,j] <- rho[j,i] <- polychor(data_mixed[,i], data_mixed[,j])
        
        ### Fan et.al. ### population kendall's tau
        na[i,j] <- na[i,j] <- sum((data_mixed[,i] == 1)*(data_mixed[,j] == 1))
        nb[i,j] <- nb[i,j] <- sum((data_mixed[,i] == 1)*(data_mixed[,j] == 0))
        nc[i,j] <- nc[i,j] <- sum((data_mixed[,i] == 0)*(data_mixed[,j] == 1))
        nd[i,j] <- nd[i,j] <- sum((data_mixed[,i] == 0)*(data_mixed[,j] == 0))
    
        tau[i,j] <- tau[j,i] <- 2*(na[i,j]*nd[i,j] - nb[i,j]*nc[i,j])/(n[1]*(n[1]-1))
        #tau[i,j] <- tau[j,i] <- cor(as.numeric(data_mixed[,i]),as.numeric(data_mixed[,j]),method="kendall")
        hatR[i,j] <- hatR[j,i] <- sin(pi*tau[i,j])
        ###
      }
    }
  }    
  
  if (!is.positive.definite(rho)) {
    rho_pd <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
  } else {
    rho_pd <- rho
  }
  
  if (!is.positive.definite(hatR)) {
    hatR_pd <- as.matrix(nearPD(hatR, corr = T, keepDiag = T)$mat)
  } else {
    hatR_pd <- hatR
  }
  #### implement Lglasso with BIC #####
  result <- huge(rho_pd,nlambda=nlam,method="glasso",verbose=FALSE)
  result_tau <- huge(hatR_pd,nlambda=nlam,method="glasso",verbose=FALSE)
  
  #### huge.select is not available when input is a covariance matrix!
  
  ####function for eBIC
  first <- omega.select(x=result)
  second <- omega.select(x=result_tau)
  
 # adj <- (abs(At)>0.05)
  adj_1 <- (abs(first) > .05)
  adj_2 <- (abs(second) > .05)
  adj_0 <- (Omega > 0) 
  
  frobenius1[run] <-  base::norm((first-Omega), type = "F")
  frobenius2[run] <-  base::norm((second-Omega), type = "F")
  
  fpr_at1[run] <- sum((adj_0 == 0 & adj_1 != 0)[lower.tri((adj_0 == 0 & adj_1 != 0))]) / (d[1]*(d[1]-1)/2 - edgenumber[run])                            
  fpr_at2[run] <- sum((adj_0 == 0 & adj_2 != 0)[lower.tri((adj_0 == 0 & adj_1 != 0))]) / (d[1]*(d[1]-1)/2 - edgenumber[run])                            
  
  tpr_at1[run] <- sum((adj_0 != 0 & adj_1 != 0)[lower.tri((adj_0 != 0 & adj_1 != 0))]) / edgenumber[run]
  tpr_at2[run] <- sum((adj_0 != 0 & adj_2 != 0)[lower.tri((adj_0 != 0 & adj_2 != 0))]) / edgenumber[run]
  
  #### get FPR and TPR 
  fp <- sapply(seq_along(result$lambda), function(i) sum((adj_0 == 0 & result$icov[[i]] != 0)[lower.tri((adj_0 == 0 & result$icov[[i]] != 0))]))
  fpr <- fp/(d[1]*(d[1]-1)/2 - edgenumber[run])
  tp <- sapply(seq_along(result$lambda), function(i) sum((adj_0 != 0 & result$icov[[i]] != 0)[lower.tri((adj_0 != 0 & result$icov[[i]] != 0))]))
  tpr <- tp/edgenumber[run]
  
  roc <- ggplot(as.data.frame(cbind(fpr,tpr)), aes(x=fpr, y=tpr)) +
  xlim(0, 1) +
  ylim(0,1) +  
  geom_line()
}

table_fan <- round(as_tibble(rbind(c('polychoric' = mean(frobenius1),'sd1' = sd(frobenius1),'kendalls'= mean(frobenius2), 'sd2' = sd(frobenius2)),
c(mean(fpr_at1),sd(fpr_at1), mean(fpr_at2),sd(fpr_at2)),
c(mean(tpr_at1),sd(tpr_at1), mean(tpr_at2),sd(tpr_at2)))),4)

stargazer(table_fan, summary = F)









### some nice graph plot
network <- graph_from_adjacency_matrix(adj_0, weighted=T, mode="undirected", diag=F)
# plot
par(bg="grey13", mar=c(0,0,0,0))
plot(network, 
     vertex.size=12,
     vertex.color="#0065bd", 
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent"
)

network2 <- graph_from_adjacency_matrix(adj_1, weighted=T, mode="undirected", diag=F)
par(bg="grey13", mar=c(0,0,0,0))
plot(network, 
     vertex.size=12,
     vertex.color="#0065bd", 
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent"
)

