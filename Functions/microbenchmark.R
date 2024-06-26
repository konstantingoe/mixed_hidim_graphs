rm(list = ls())

# setwd("/home/goebler/linux/mixed_hidim_graphs")
source("./Packages/packages.R")
source("./Functions/simulation_functions.R")

graph_learner <- function(data, param = .1, nlam = 30, method = c("poly", "bridge", "mle"), only_corrmat = FALSE) {
    n <- nrow(data)
    if (method == "poly") {
        sigma <- hatR_poly(data, verbose = F)
    } else if (method == "bridge") {
        sigma <- hatR_feng(data, verbose = F)
    } else if (method == "mle") {
        sigma <- hatR_mle(data, verbose = F)
    }
    if (only_corrmat) {
        return(sigma)
    } else {
        omega_path <- huge::huge(sigma, nlambda = nlam, method = "glasso", verbose = FALSE)
        omega_hat <- ebic_selection(x = omega_path, n = n, s = sigma, param = param)
        adj_estimate <- abs(omega_hat) > 0
        return(adj_estimate)
    }
}

graph_learner_timing <- function(dimension = 50, sim_runs = 100) {
    data <- generate_data(n = 300, d = dimension, g = cubic)
    data_1 <- data[[2]]

    data <- NULL

    data_mixed <- make_ordinal(data_1, general = TRUE)
    test_bm <- microbenchmark::microbenchmark(
        unit = "milliseconds",
        times = sim_runs,
        oracle = {
            graph_learner(data = data_1, method = "poly", only_corrmat = TRUE)
        },
        poly = {
            graph_learner(data = data_mixed, method = "poly", only_corrmat = TRUE)
        },
        bridge = {
            graph_learner(data = data_mixed, method = "bridge", only_corrmat = TRUE)
        }
    )
    return(summary(test_bm))
}


d_seq <- c(50, 250, 750)
final_seq <- rep(d_seq, each = 3)

corr_mat_benchmark <- data.frame(cbind(
    final_seq,
    do.call(rbind, lapply(
        d_seq, function(k) {
            graph_learner_timing(
                dimension = k,
                sim_runs = 5
            )
        }
    ))
))

print(corr_mat_benchmark)

save(corr_mat_benchmark, file = "./simulation_results/corrmat_speed.RData")

load("./simulation_results/corrmat_speed.RData")

results <- corr_mat_benchmark[, c("expr", "final_seq", "mean")]
results$mean <- results$mean / 1000
names(results) <- c("Method", "Dimension", "Time (s)")

knitr::kable(results, "latex", booktabs = T, digits = 4) %>%
    cat(., file = "./paper/High-Dimensional Mixed Graphs EJS/Tables/timing.tex")
