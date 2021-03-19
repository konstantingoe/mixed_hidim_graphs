#### Simulation with ordinal data Drton procedure ####

rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

library(igraph)

create_component <- function(graph = graph, number_hubs = number_hubs, hub_node_degree = hub_node_degree) {
  potential_edges <- seq_along(V(graph))
  degree <- floor(mean(degree(g)))
  
  hub_nodes <- sample(potential_edges, number_hubs, replace = F)
  edge_set <- sample(potential_edges[-hub_nodes], (hub_node_degree - degree)*number_hubs, replace = F)
  edge_hub <- base::split(edge_set, ceiling(seq_along(edge_set)/(hub_node_degree - degree)))
  
  merge_vecs <- function(vec1, vec2){
    stopifnot(length(vec1)==length(vec2))
    merge <- c(vec1[1], vec2[1])
    for (i in 2:length(vec1)){
      merge <- c(merge, c(vec1[i], vec2[i]))
    }
    return(merge)
  }
  
  edge_paths <- unlist(lapply(1:number_hubs, function(k) merge_vecs(rep(hub_nodes[k],hub_node_degree - degree),edge_hub[[k]])))
  graph_out <-  add_edges(graph, edge_paths)
  return(graph_out)
}


hub_graph <- create_component(graph=make_lattice(dimvector = c(10,10), 
                              circular = T), number_hubs = 3, hub_node_degree = 20) 

new_ad <- bdiag(as_adjacency_matrix(create_component(graph=make_lattice(dimvector = c(10,10), 
                                                    circular = T), number_hubs = 3, hub_node_degree = 20)),
                as_adjacency_matrix(create_component(graph=make_lattice(dimvector = c(10,10), 
                                                    circular = T), number_hubs = 3, hub_node_degree = 20))
                )

for (i in 3:10){
  new_ad <- bdiag(new_ad, as_adjacency_matrix(create_component(graph=make_lattice(dimvector = c(10,10),circular = T), number_hubs = 3, hub_node_degree = 20)))
}

isSymmetric.matrix(as.matrix(new_ad))
plot(graph_from_adjacency_matrix(new_ad))





#### check results for large sparsity parameter


data <- simulate.mixed.data(N=100, pfrac=c(50,0,0,0,0,0,50), sparsity = .5)
data_mixed <- as.data.frame(data[[1]])

for (i in 51:ncol(data_mixed)) {
  data_mixed[,i] <- factor(data_mixed[,i], ordered = T)
}

Omega <- data[[2]]
print(paste0("We have ",edgenumber(Omega), " edges in the graph"))
data_0 <- data[[3]]

rho_latent <- mixed.omega(data_0)
rho <- mixed.omega(data_mixed)

huge.result <- huge(rho,nlambda=50,method="glasso",verbose=FALSE)
huge.result.latent <- huge(rho_latent,nlambda=50,method="glasso",verbose=FALSE)

AUC_hat <- AUC(truth = Omega, huge_obj = huge.result)
AUC_hat_latent <- AUC(truth = Omega, huge_obj = huge.result.latent)
AUC_hat
AUC_hat_latent







data <- generate.data(n=100, d=100, n_E=choose(100,2)/2)
data_0 <- data[[1]]
Omega <- data[[2]]
data_mixed <- make.ordinal(data_0, countvar = F, n_O = 2)

rho_latent <- mixed.omega(data_0)
rho <- mixed.omega(data_mixed)

huge.result <- huge(rho,nlambda=50,method="glasso",verbose=FALSE)
huge.result.latent <- huge(rho_latent,nlambda=50,method="glasso",verbose=FALSE)

AUC_hat <- AUC(truth = Omega, huge_obj = huge.result)
AUC_hat_latent <- AUC(truth = Omega, huge_obj = huge.result.latent)
AUC_hat
AUC_hat_latent



