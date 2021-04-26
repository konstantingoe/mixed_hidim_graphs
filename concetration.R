
library(GGally)
library(network)
library(sna)
library(ggplot2)
# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]


ggnet2(net, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)

ggsave("example_graph.pdf")









### here input a simple adjecency matrix! 
nlam=50
data <- generate.data(t=.15, n = 200, d = 50, mode = "fan")
data_0 <- data[[1]]
Omega <- data[[2]]

rho_latent <- mixed.omega(data_0)
results_latent <- glasso.results(Sigma=rho_latent, Omega=Omega, nlam=50, n=200, matexport = T, param = .25)

data_mixed <- make.ordinal.general(data_0, namevector = c("binary" = T, "ordinal" = T, "poisson" = T), unbalanced = .2, low = .05, high = .1)
### learn sample correlation matrix   
rho <- mixed.omega(data_mixed)

results <- glasso.results(Sigma=rho, Omega=Omega, nlam=50, n=200, matexport = T,  param = .25)

true_adjecency <- abs(Omega) > 0
latent_adjecency <- abs(results_latent$`Estimate d Precision Matrix`) > 0
estimated_adjecency <- abs(results$`Estimated Precision Matrix`) > 0

net_true = network(true_adjecency, directed = FALSE)
net_latent = network(latent_adjecency, directed = FALSE)
net_estimated = network(estimated_adjecency, directed = FALSE)

ggnet2(net_true, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)
ggnet2(net_latent, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)
ggnet2(net_estimated, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)

