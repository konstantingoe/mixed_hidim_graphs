rm(list = ls())
source("Packages/packages.R")
source("Functions/functions.R")

num_levels <- 4
data <- rmvnorm(1000, c(0, 0), matrix(c(1, -.27, -.27, 1), 2, 2))
transformed_data <- .5 * data^5
cor(transformed_data)

x <- transformed_data[, 1]
y <- transformed_data[, 2]



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

y_disc <- latent_ordinalize(y)

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

cor(x, y) # sample correlation

polyserial(x, y_disc, ML = TRUE, std.err = TRUE)
polyserial(x, y_disc, ML = F, std.err = F)
polyserial_mle(cont_var = x, disc_var = y_disc)


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

system.time(polyserial_closedform_mle(cont_var = x, disc_var = y_disc, nonparanormal = TRUE))
system.time(adhoc_lord_sim(cont = x, disc = y_disc))
system.time(polyserial_mle(cont_var = x, disc_var = y_disc, nonparanormal = TRUE))

adhoc_lord_sim(cont = x, disc = y_disc)
polyserial_closedform_mle(cont_var = x, disc_var = y_disc, nonparanormal = TRUE)
polyserial_mle(cont_var = x, disc_var = y_disc, nonparanormal = TRUE)

corr_comp <- function(rho, transform = function(x) 5 * x^5) {
    data <- mvtnorm::rmvnorm(1000, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2))
    transformed_data <- transform(data) # .5 * data^5
    x <- transformed_data[, 1]
    y <- latent_ordinalize(transformed_data[, 2])

    adhoc <- adhoc_lord_sim(cont = x, disc = y)
    closed_form <- polyserial_closedform_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)
    obj_func <- polyserial_mle(cont_var = x, disc_var = y, nonparanormal = TRUE)

    return(c("adhoc" = adhoc, "closed_form" = closed_form, "obj_func" = obj_func))
}

run_corr_sim <- function(sim, rho = .5) {
    intermdediate <- do.call(rbind, lapply(seq_len(sim), function(s) corr_comp(rho = rho)))
    mean_sims <- colMeans(intermdediate)
    sd_sims <- apply(intermdediate, 2, sd)
    new_names <- paste(names(sd_sims), "sd", sep = "_")
    names(sd_sims) <- new_names
    return(c(mean_sims, sd_sims))
}

run_corr_sim(sim = 100, rho = .99)
corr_seq <- seq(-.98, .98, length.out = 100)

final_corr <- data.frame(cbind(
    corr_seq,
    do.call(rbind, lapply(
        corr_seq,
        function(k) run_corr_sim(sim = 100, rho = k)
    ))
))

final_corr_melt <- melt(final_corr, id.vars = "corr_seq")

ggplot(final_corr_melt, aes(y = corr_seq, x = value, colour = variable)) +
    geom_point() +
    geom_line() +
    geom_abline() +
    geom_errorbar(aes(xmin = corr_seq - value, xmax = corr_seq + value, colour = variable), width = .3)
