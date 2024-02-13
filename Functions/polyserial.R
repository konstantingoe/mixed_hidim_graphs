rm(list = ls())

setwd("/Users/kgoebler/hume_revisions/mixed_hidim_graphs")
source("./Packages/packages.R")
source("./Functions/functions.R")

# The palette with grey:
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

latent_ordinalize <- function(x, num_levels = 2, unbalanced = FALSE) {
    if (unbalanced) {
        levels <- runif(num_levels)
        sum_levels <- sum(levels)
        levels_norm <- sort(levels) / sum_levels
        cum_levels <- c(0, cumsum(levels_norm))
    } else {
        cum_levels <- seq(0, 1, length.out = num_levels + 1)
    }
    u <- pnorm(scale(x))
    ordinal <- cut(u,
        breaks = cum_levels,
        include.lowest = T,
        ordered_result = T,
        labels = 1:(length(cum_levels) - 1)
    )

    return(factor(ordinal, ordered = TRUE))
}

polyserial_mle <- function(cont_var, disc_var, maxcor = 0.999, nonparanormal = FALSE) {
    x <- cont_var
    if (nonparanormal) {
        x <- f_hat(cont_var)
    }
    y <- disc_var
    f <- function(pars) {
        rho <- pars[1]

        if (abs(rho) > maxcor) {
            rho <- sign(rho) * maxcor
        }

        if (length(pars) == 1) {
            cts <- c(-Inf, cuts, Inf)
        } else {
            cts <- c(-Inf, pars[-1], Inf)
        }

        if (any(diff(cts) < 0)) {
            return(Inf)
        }
        tau <- (matrix(cts, n, s + 1, byrow = TRUE) -
            matrix(rho * z, n, s + 1)) / sqrt(1 - rho^2)
        -sum(log(dnorm(z) * (pnorm(tau[cbind(indices, y + 1)]) -
            pnorm(tau[cbind(indices, y)])))) # eq (20) in Olsson et.al. (1982)
    }

    z <- scale(x) # (x-mu)/sd
    tab <- table(y)

    n <- sum(tab)
    indices <- 1:n
    s <- length(tab)

    cuts <- qnorm(cumsum(tab) / n)[-s]

    y <- as.numeric(as.factor(y))

    rho <- sqrt((n - 1) / n) * sd(y) * cor(x, y) / sum(dnorm(cuts))

    result <- optim(c(rho, cuts), f, hessian = F, method = "BFGS")
    return(result$par[1])
}

f_hat <- function(x) {
    n <- length(x)
    npn_thresh <- 1 / (4 * (n^0.25) * sqrt(pi * log(n)))
    x <- qnorm(pmin(
        pmax(rank(x) / n, npn_thresh),
        1 - npn_thresh
    ))
    x <- x / sd(x)
    # rm(n, npn.thresh)
    # gc()
    return(x)
}

polyserial_closedform_mle <- function(cont_var, disc_var, maxcorr = 0.999, nonparanormal = FALSE) {
    x <- cont_var
    if (nonparanormal) {
        x <- f_hat(cont_var)
    }
    y <- as.numeric(disc_var)
    z <- scale(x)
    n <- length(x)

    cummarg_propotions <- c(0, cumsum(table(y) / n))
    thresholds <- qnorm(cummarg_propotions)
    rho <- sqrt((n - 1) / n) * sd(y) * cor(x, y) / sum(dnorm(thresholds))


    polyserial_closedform <- function(rho) {
        thresholds_star <- function(rho, j_minus_1 = TRUE) {
            if (j_minus_1) {
                return((thresholds[y] - rho * z) / sqrt(1 - rho^2))
            } else {
                return((thresholds[(y + 1)] - rho * z) / sqrt(1 - rho^2))
            }
        }

        first_summand <- ifelse(
            dnorm(thresholds_star(rho, j_minus_1 = FALSE)) != 0,
            dnorm(thresholds_star(rho, j_minus_1 = FALSE)) * (thresholds[(y + 1)] * rho - z),
            0
        )

        second_summand <- ifelse(
            dnorm(thresholds_star(rho, j_minus_1 = TRUE)) != 0,
            dnorm(thresholds_star(rho, j_minus_1 = TRUE)) * (thresholds[y] * rho - z),
            0
        )

        pre_factor <- (1 - rho^2)^-(3 / 2)

        conditional_prob_inv <- (1 / (pnorm(thresholds_star(rho, j_minus_1 = FALSE))
        - pnorm(thresholds_star(rho, j_minus_1 = TRUE))))


        1 / n * sum(pre_factor * conditional_prob_inv * (first_summand - second_summand))
    }

    res <- nleqslv::nleqslv(rho, polyserial_closedform)
    return(res$x)
}

corr_comp <- function(rho, transform = function(x) 5 * x^5) {
    data <- mvtnorm::rmvnorm(1000, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2))
    transformed_data <- transform(data) # .5 * data^5
    x <- transformed_data[, 1]
    y <- latent_ordinalize(transformed_data[, 2])

    adhoc <- adhoc_lord_sim(cont = x, disc = y)
    closed_form <- polyserial_closedform_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)
    obj_func <- polyserial_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)

    return(c("ad_hoc" = adhoc, "foc" = closed_form, "log_likelihood" = obj_func))
}

run_corr_sim <- function(sim, rho = .5) {
    intermdediate <- do.call(rbind, lapply(seq_len(sim), function(s) corr_comp(rho = rho)))
    mean_sims <- colMeans(intermdediate)
    sd_sims <- apply(intermdediate, 2, sd)
    new_names <- paste(names(sd_sims), "sd", sep = "_")
    names(sd_sims) <- new_names
    return(c(mean_sims, sd_sims))
}

corr_seq <- seq(0, .98, length.out = 50)
test_corr <- data.frame(cbind(
    corr_seq,
    do.call(rbind, lapply(
        corr_seq,
        function(k) run_corr_sim(sim = 100, rho = k)
    ))
))


test_long_df <- melt(test_corr, id.vars = "corr_seq", measure.vars = c("ad_hoc", "foc", "log_likelihood"))
interim <- melt(test_corr, id.vars = "corr_seq", measure.vars = c("ad_hoc_sd", "foc_sd", "log_likelihood_sd"))
melt_df <- cbind(test_long_df, interim[, "value"])
names(melt_df) <- c("corr_seq", "variable", "value", "value_sd")

ggplot(melt_df) +
    geom_line(aes(y = corr_seq - value, x = corr_seq, colour = variable)) +
    geom_point(aes(y = corr_seq - value, x = corr_seq, colour = variable)) +
    geom_ribbon(aes(
        ymin = corr_seq - value - value_sd, ymax = corr_seq - value + value_sd,
        x = corr_seq, colour = variable
    ), alpha = .02) +
    ylab(TeX("$\\Delta$ betweem $\\Sigma^{*}_{jk}$ and $\\hat{\\Sigma}_{jk}$")) +
    xlab(TeX("$\\Sigma^{*}_{jk}$")) +
    labs(colour = "Estimation method") +
    scale_colour_manual(values = cbPalette) +
    ylim(-.1, .1) +
    theme_minimal(base_size = 16)

ggplot2::ggsave(filename = "case_2_difference.pdf", path = "./paper/High-Dimensional Mixed Graphs EJS/Figures/")

#### Some microbenchmarks ####


corr_timing <- function(rho, transform = function(x) 5 * x^5, sample_size = 1000) {
    data <- mvtnorm::rmvnorm(sample_size, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2))
    transformed_data <- transform(data) # .5 * data^5
    x <- transformed_data[, 1]
    y <- latent_ordinalize(transformed_data[, 2])
    test_bm <- microbenchmark::microbenchmark(
        unit = "milliseconds",
        ad_hoc = {
            adhoc_lord_sim(cont = x, disc = y)
        },
        foc = {
            polyserial_closedform_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)
        },
        log_likelihood = {
            polyserial_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)
        }
    )

    return(summary(test_bm))
}

n_seq <- seq(from = 50, to = 10000, by = 50)
final_seq <- rep(n_seq, each = 3)

sample_size_corr <- data.frame(cbind(
    final_seq,
    do.call(rbind, lapply(
        n_seq, function(k) {
            corr_timing(
                rho = .87, sample_size = k
            )
        }
    ))
))

microbench_a <- ggplot(sample_size_corr) +
    geom_line(aes(y = median, x = final_seq, colour = expr)) +
    geom_ribbon(aes(
        ymin = lq, ymax = uq,
        x = final_seq, colour = expr
    ), alpha = .2) +
    ylab("milliseconds") +
    xlab("sample size") +
    labs(colour = "Estimation method") +
    scale_colour_manual(values = cbPalette) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")

ggplot2::ggsave(filename = "case_2_speed.pdf", path = "./paper/High-Dimensional Mixed Graphs EJS/Figures/")

corr_seq <- seq(-.89, .98, length.out = 200)
corr_seq_final <- rep(corr_seq, each = 3)
corr_corr <- data.frame(cbind(
    corr_seq_final,
    do.call(rbind, lapply(
        corr_seq, function(k) {
            corr_timing(
                rho = k, sample_size = 1000
            )
        }
    ))
))


microbench_b <- ggplot(corr_corr) +
    geom_line(aes(y = median, x = corr_seq_final, colour = expr)) +
    geom_ribbon(aes(
        ymin = lq, ymax = uq,
        x = corr_seq_final, colour = expr
    ), alpha = .2) +
    ylab(TeX("milliseconds")) +
    xlab(TeX("$\\Sigma^{*}_{jk}$")) +
    labs(colour = "Method") +
    scale_colour_manual(values = cbPalette) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "top")

ggplot2::ggsave(filename = "case_2_speed_corr.pdf", path = "./paper/High-Dimensional Mixed Graphs EJS/Figures/")

ggpubr::ggarrange(microbench_a, microbench_b,
    # labels = c("A", "B"),
    ncol = 2, nrow = 1
)

ggplot2::ggsave(filename = "case_2_speed_comp.pdf", path = "./paper/High-Dimensional Mixed Graphs EJS/Figures/")
