rm(list = ls())
source("Packages/packages.R")
source("Functions/functions.R")

num_levels <- 4
data <- rmvnorm(1000, c(0, 0), matrix(c(1, -.27, -.27, 1), 2, 2))
transformed_data <- .5 * data^5
cor(transformed_data)

x <- transformed_data[, 1]
y <- transformed_data[, 2]

levels <- runif(num_levels)
sum_levels <- sum(levels)
levels_norm <- sort(levels) / sum_levels
cum_levels <- c(0, cumsum(levels_norm))
u <- pnorm(scale(y))
ordinal <- cut(u,
    breaks = cum_levels,
    include.lowest = T,
    ordered_result = T,
    labels = 1:(length(cum_levels) - 1)
)

y_disc <- factor(ordinal, ordered = TRUE)

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
polyserial_closedform_mle(cont_var = x, disc_var = y_disc, nonparanormal = FALSE)
