rm(list = ls())
source("packages.R")
source("functions.R")

set.seed(1221)

#mypath<- "/soep/kgoebler/data"
mydata <- import("/Users/kgoebler/Desktop/SOEP Projekte/Wealth Acrticle/out/tw_cleaned.dta")


#potentially <- set_na(mydata$residence_debt_filter, na =c("Does not apply" = -2), as.tag = T)
# first impute filter information and then based on these impute wealth components:

##### Filter Imputation


factorvars <- c("sex", "bula_sta", "job_training", "selfempl", "superior", "second_empl", "citizen", "famstd",
                "owner", "residence_debt_filter", "other_estate", "other_estate_debt_filter",
                "assets_filter", "building_contract_filter", "life_insure_filter", "business_holdings_filter",
                "vehicles_filter", "tangibles_filter", "consumer_debt_filter", "education_debt_filter",
                "gborn", "kidsu16", "partner", "inheritance_dummy", "emp_status", "labormarket_dummy")

orderedfactorvars <- c("hhgr", "bik", "ggk", "wuma3", "compsize", "overtime", "residence_value_limits",
                       "residence_debt_limits", "other_estate_value_limits", "other_estate_debt_limits", 
                       "assets_limits", "building_contract_limits", "life_insure_limits", "business_holdings_limits",
                       "vehicles_limits", "tangibles_limits", "consumer_debt_limits", "education_debt_limits",
                       "education", "housecond", "nkids")

continousvars <- setdiff(names(mydata),c(factorvars,orderedfactorvars))

mydata$hhgr <- ifelse(mydata$hhgr > 4, 4, mydata$hhgr)

for (i in 1:length(factorvars)){
  mydata[,factorvars[i]] <- factor(mydata[,factorvars[i]], ordered = F)
}

for (i in 1:length(orderedfactorvars)){
  mydata[,orderedfactorvars[i]] <- factor(mydata[,orderedfactorvars[i]], ordered = T)
}

gg_miss_var(mydata, show_pct = TRUE)

data <- mydata

keepvec <- c("age", "sqmtrs", "hhnetto", "wage_gross_m", "inheritance_dummy",
             "sex", "selfempl", "citizen", "famstd",
             "owner", "kidsu16", "partner", "emp_status",
             "hhgr", "bik", "ggk", "wuma3", "overtime", 
             "education", "housecond", "nkids")

data.try <- dplyr::select(data, all_of(keepvec))

for (i in which(sapply(data.try, is.numeric))){
  data.try <- data.try[data.try[,i] > 0,]
}
data.try <- na.omit(data.try)

data.try2 <- data.try

data.try2 <- data.try2 %>% 
  mutate(hhnetto = (log(hhnetto) - mean(log(hhnetto)))/sd(log(hhnetto)),
         wage_gross_m = (log(wage_gross_m) - mean(log(wage_gross_m)))/sd(log(wage_gross_m)))
data.try2 <- dplyr::select(data.try2, -c("overtime", "emp_status", "kidsu16"))  
  
graph.try <- mixed.undir.graph(data.try2, nlam = 50, param = .1)

max(abs(graph.try$`Estimated Precision Matrix`))

net = network(graph.try$`Adjecency Matrix`, directed = FALSE)
network.vertex.names(net) = colnames(data.try2)

ggnet2(net, size = 12, label = TRUE, label.size = 5, label.alpha = 0.75)

which(abs(graph.try$`Estimated Precision Matrix`) > 1, arr.ind = T)

graph.try$`Sample Correlation Matrix`[13,11]

names(data.try2[,c(11,13)])

table(data.try2$kidsu16, data.try2$hhgr)


### summary of things that could go wrong: empty catogies are problematic, always standardize, if possible log-standardize variables
### build in warning message when correlation exceeds abs(.85)
