rm(list = ls())
source("packages.R")
source("functions.R")
source("applied_functions.R")

library(rio)
data <- import("/Users/kgoebler/Desktop/TUM Promotion/Projekte/Mixed GMs/Data/TCGA/brain.csv")

set.seed(1234)
gene.data <- data[,-1]

factor_ids <- sapply(gene.data, function(id) length(unique(id)) < 10)
gene.data[,factor_ids] <- lapply(gene.data[,factor_ids], factor)

log.transform <- function(data, epsilon = .001){
  epsilon = epsilon
  data <- (log(data + epsilon) - mean(log(data + epsilon)))/sd(log(data + epsilon))
  return(data)
}

contdata <- sapply(gene.data, is.numeric)
test.data <- gene.data
test.data[,contdata] <- sapply(seq_along(contdata), function(k)
                            log.transform(gene.data[,k]))


#length(which(sapply(gene.data, function(id) length(unique(id)) < 10)))/ncol(gene.data)
#roughly 50% "discrete" variables when 20 is the threshold
#still 38% discrete variables when 10 is the threshold

indices <- sample(1:ncol(gene.data),size = 200)
#length(which(sapply(gene.data[,indices], function(id) length(unique(id)) < 20)))/300



testrun <- mixed.nonpara.graph(gene.data[,indices], verbose = F, param = .5)

d.increasing <- seq(from = 50, to = 1000, by = 50)
#d.increasing <- seq(from = 50, to = 300, by = 50)

indices.list <- lapply(seq_along(d.increasing), function(i) sample(1:ncol(gene.data),size = d.increasing[i]))

corrmat <- sapply(1:length(indices.list), function(t) system.time(mixed.omega.paranormal(gene.data[,indices.list[[t]]], verbose = F))[3])

plotdata <- data.frame(dimension = d.increasing, runtime_in_seconds = corrmat)

ggplot(data=plotdata, aes(x=runtime_in_seconds, y=dimension, group=1)) +
  geom_line()+
  geom_point()

