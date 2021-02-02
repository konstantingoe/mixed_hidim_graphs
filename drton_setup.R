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
  new_ad <- bdiag(new_ad, as_adjacency_matrix(create_component(graph=make_lattice(dimvector = c(10,10),                                                                                circular = T), number_hubs = 3, hub_node_degree = 20)))
}

isSymmetric.matrix(as.matrix(new_ad))
plot(graph_from_adjacency_matrix(new_ad))

