
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




sum(ell)