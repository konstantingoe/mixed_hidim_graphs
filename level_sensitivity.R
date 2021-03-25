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
nlam <- 50 # number of tuning parameters for graphical lasso
countvar <- F
latent <- F
matexport = F

unbalanced.grid <- seq(from=0, to=1, by = .05)
reps <- 50
lowdim.result <- lapply(seq_along(unbalanced.grid), function(k) lapply(1:reps, function(i) unbalanced.run(mode = "fan", n=n[1], d=d[1], sparsity = .1, nlam=nlam, matexport = F, namevector = c("binary" = T, "ordinal" = F, "poisson" = F), 
                       unbalanced = unbalanced.grid[k], low = .05, high = .1)))

plotgrid <- rep(unbalanced.grid[1], reps)
for (k in 2:length(unbalanced.grid)){
  plotgrid <- c(plotgrid, rep(unbalanced.grid[k], reps))
}

frobenius <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
              sapply(1:reps, function(j) lowdim.result[[k]][[j]][[1]])))
AUC <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
                sapply(1:reps, function(j) lowdim.result[[k]][[j]][[4]])))

plotdata <- as.data.frame(cbind(plotgrid, "frobenius" = frobenius, "AUC" = AUC))

p <- ggplot(plotdata, aes(x=factor(plotgrid), y=frobenius)) +
     geom_boxplot() +
     xlab("Unbalanced Proportion") +
     ylab("Frobenius Norm") +
     scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])
q <- ggplot(plotdata, aes(x=factor(plotgrid), y=AUC)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("Arean Under the Curve") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])

plot_grid(p,q)



lowdim.result <- lapply(seq_along(unbalanced.grid), function(k) lapply(1:reps, function(i) unbalanced.run(mode = "er", n=n[1], d=d[1], nlam=nlam, matexport = F, namevector = c("binary" = T, "ordinal" = F, "poisson" = F), 
                                                                                                          unbalanced = unbalanced.grid[k], low = .05, high = .1)))

plotgrid <- rep(unbalanced.grid[1], reps)
for (k in 2:length(unbalanced.grid)){
  plotgrid <- c(plotgrid, rep(unbalanced.grid[k], reps))
}

frobenius.plot <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
  sapply(1:reps, function(j) lowdim.result[[k]][[j]][[1]])))
AUC.plot <- unlist(lapply(seq_along(unbalanced.grid), function(k) 
  sapply(1:reps, function(j) lowdim.result[[k]][[j]][[4]])))

plotdata <- as.data.frame(cbind(plotgrid, "frobenius" = frobenius.plot, "AUC" = AUC.plot))

p <- ggplot(plotdata, aes(x=factor(plotgrid), y=frobenius)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("Frobenius Norm") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])
q <- ggplot(plotdata, aes(x=factor(plotgrid), y=AUC)) +
  geom_boxplot() +
  xlab("Unbalanced Proportion") +
  ylab("Area Under the Curve") +
  scale_x_discrete(breaks = levels(factor(plotdata$plotgrid))[c(rep(c(T, F),floor(length(unbalanced.grid)/2)),T)])

plot_grid(p,q)


