rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

### load data 

library(rio)
library(stringr)

data <- import("/Users/kgoebler/Desktop/Lektorat/Sabrinas Doktorarbeit/Tabelle ParoVeg_Imputation_kurz.xlsx",
               which = 1, na = "NA")

arcsinevars <- colnames(dplyr::select(data, contains("arcsin")))
renamevars <- colnames(dplyr::select(data, contains("...")))
dropvars <- c("id", renamevars, arcsinevars)

data.raw <- dplyr::select(data, -all_of(dropvars)) 


changevars <- c("crp_t01", "crp_t02", "crp_t03", "holo_t02", "ins_t03",
                "sf_t01", "sf_t02", "sf_t03")

data.raw <- data.raw %>% 
  mutate(score = ifelse(
            SCORE == "niedrig",1, ifelse(
            str_detect(SCORE, "mo"), 2, ifelse(
            SCORE == "hoch", 3, NA))),
         crp_t01 = as.numeric(ifelse(crp_t01 == "<0,3", 0, crp_t01)),
         crp_t02 = as.numeric(ifelse(crp_t02 == "<0,3", 0, crp_t02)),
         crp_t03 = as.numeric(ifelse(crp_t03 == "<0,3", 0, crp_t03)),
         ins_t03 = as.numeric(ifelse(ins_t03 == "folgt", NA, ins_t03)),
         sf_t01 = as.numeric(ifelse(sf_t01 == "<20", 0, sf_t01)),
         sf_t02 = as.numeric(ifelse(sf_t02 == "<20", 0, sf_t02)),
         sf_t03 = as.numeric(ifelse(sf_t03 == "<20", 0, sf_t03)),
         )

data.raw <- dplyr::select(data.raw, -all_of(c("SCORE", "holo_t02", "hy")))


data.final <- filter(data.raw, !(is.na(t02) & is.na(t03)))

### imputation: 

library(missForest)

misstoohigh <- which(sapply(1:ncol(data.final), function(i) sum(is.na(data.final[,i]))/nrow(data.final)) > .2)

data.imp.pre <- dplyr::select(data.final, -all_of(c("t01", "t02","t03", colnames(data.final[,misstoohigh]))))
data.imputed <- missForest(data.imp.pre)

data.complete <- cbind(data.final[,c("t01", "t02","t03")], data.imputed$ximp)


#### Descriptive Analysis ####

# decriptive vars

descvars <- c("age", "sexe", "ra", "m", "py", "nu", "bu", "ge", "gr")

c(mean(filter(data.complete, g == 1)[,descvars[1]]), mean(filter(data.complete, g == 2)[,descvars[1]]),
  wilcox.test(filter(data.complete, g == 1)[,descvars[1]], filter(data.complete, g == 2)[,descvars[1]], exact = F)$p.value)

c((mean(filter(data.complete, g == 1)[,descvars[2]]) -1 )*100, mean((filter(data.complete, g == 2)[,descvars[2]]) -1 )*100,
  chisq.test(table(data.complete[,c("g",descvars[2])]))$p.value)







#### Strategy: Estimate Graph for Baseline, T1, and T2
#### and take differential edges that are present in each graph 
### get baseline vars

timeinvariantnames <- c("g", "age", "sexe", "th", "py", "zig", "alk", "mets", "m", 
                        "medzu", "medred", "nu", "fa", "path", "zvp", "muhy1", "muhy2", "muhy3")

baselinenames <- str_subset(colnames(data.complete), "01")

data.baseline <- dplyr::select(data.complete, all_of(c(baselinenames[-which(baselinenames == "t01")], timeinvariantnames)))

graph.baseline <- mixed.undir.graph(data.baseline, nlam = 50)

library(GGally)
library(ggnetwork)
library(network)
### nice even at this stage this is what we want!
net.baseline <-  network(graph.baseline$`Adjecency Matrix`, directed = FALSE)
network.vertex.names(net.baseline) <- colnames(data.baseline)
ggnet2(net.baseline, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)


t2names <- str_subset(colnames(data.complete), "02")
data.t2 <- dplyr::select(data.complete, all_of(c(t2names, timeinvariantnames)))

graph.t2 <- mixed.undir.graph(data.t2, nlam = 50)
### nice even at this stage this is what we want!
net.t2 <-  network(graph.t2$`Adjecency Matrix`, directed = FALSE)
network.vertex.names(net.t2) <- colnames(data.t2)
ggnet2(net.t2, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)





