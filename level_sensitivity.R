#### let us see how results behave w.r.t. class balance ####

#### Extended Comparison out Method and Fan et al 

#### Simulation with ordinal data ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

sim <- 100 # simulation runs
n <- c(200,200,300,600) # sample size include high dimension only on cluster
d <- c(50,250,750, 1500) # dimensionality --> include high dimension (1500) only on cluster 
n_E <- c(200,250,750,1500) # sparsity level of the graph: amount of edges we want to introduce 
t <- .15 # signal strength
nlam <- 30 # number of tuning parameters for graphical lasso
countvar <- F
latent <- F
matexport = F
reps <- 20
unbalanced.grid <- seq(from=0, to=1, by = .1)

cat(paste0("number of available cores: ",detectCores()))

if (detectCores() >= 50){
  numCores <-  50
} else {
  numCores <- detectCores()
}

print("Start with d=750")
plan(multisession, workers = numCores) ## Run in parallel on Linux cluster
  
lowdim.result.fan <- lapply(seq_along(unbalanced.grid), function(k) future_lapply(future.seed = T, 1:reps, function(i) unbalanced.run(mode = "fan", n=n[1], d=d[3], sparsity = .1, nlam=nlam, matexport = F, namevector = c("binary" = T, "ordinal" = F, "poisson" = F), 
                       unbalanced = unbalanced.grid[k], low = .05, high = .1)))

lowdim.result.er <- lapply(seq_along(unbalanced.grid), function(k) future_lapply(future.seed = T, 1:reps, function(i) unbalanced.run(mode = "er", n=n[1], d=d[3], nlam=nlam, matexport = F, namevector = c("binary" = T, "ordinal" = F, "poisson" = F), 
                       unbalanced = unbalanced.grid[k], low = .05, high = .1)))

plan(sequential)

plotgrid <- rep(unbalanced.grid[1], reps)
for (k in 2:length(unbalanced.grid)){
  plotgrid <- c(plotgrid, rep(unbalanced.grid[k], reps))
}

frobenius.fan <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
                 sapply(1:reps, function(j) lowdim.result.fan[[k]][[j]][[1]])))
AUC.fan <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
           sapply(1:reps, function(j) lowdim.result.fan[[k]][[j]][[4]])))

frobenius.er <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
                sapply(1:reps, function(j) lowdim.result.er[[k]][[j]][[1]])))
AUC.er <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
          sapply(1:reps, function(j) lowdim.result.er[[k]][[j]][[4]])))

plotdata <- as.data.frame(cbind(plotgrid, "frobenius_fan" = frobenius.fan, "AUC_fan" = AUC.fan, "frobenius_er" = frobenius.er, "AUC_er" = AUC.er))

p <- ggplot(plotdata, aes(x=factor(plotgrid), y=frobenius_fan)) +
     geom_boxplot() +
     xlab("Unbalanced Proportion") +
     ylab("Frobenius Norm") +
     scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])
q <- ggplot(plotdata, aes(x=factor(plotgrid), y=AUC_fan)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("AUC") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])

r <- ggplot(plotdata, aes(x=factor(plotgrid), y=frobenius_er)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("Frobenius Norm") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])
s <- ggplot(plotdata, aes(x=factor(plotgrid), y=AUC_er)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("AUC") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])

p_q_r_s <- plot_grid(p,q,r,s, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2)

ggsave("fan_gen.pdf", plot = p_q_r_s, width = 20, height = 10, units = "cm")


