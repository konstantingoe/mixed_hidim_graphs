source("mixed_hidim_graphs/Packages/packages.R")

t <- .15
n <- 200
d <- 50
latent <- FALSE
nlam <- 30
countvar <- TRUE
# function for euclidean norm
euclid_norm <- function(x) sqrt(sum(x^2))
edgenumber <- function(mat = mat, cut = 0) {
    return(
        sum((abs(mat) > (0 + cut))[lower.tri((abs(mat) > (0 + cut)))])
    )
}


benchmark_run <- function(runs = 100, n = n, d = d, param = .1, nlam = 30, g, namevector = c("binary" = T, "ordinal" = F, "poisson" = F)) {
    data <- generate_data(n = n, d = d)
    # data_0 <- data[[1]]
    data_1 <- data[[2]]
    omega <- data[[3]]

    data <- NULL

    data_mixed <- make_ordinal_general(data_1, namevector = namevector)

    ### learn sample correlation matrix
    sigma_oracle <- hatR_poly(data_1)
    sigma_fan <- hatR_fan(data_mixed)
    sigma_poly <- hatR_poly(data_mixed)
    sigma_mle <- hatR_mle(data_mixed)

    data_mixed <- NULL
    corr_list <- list(sigma_oracle, sigma_fan, sigma_poly, sigma_mle)

    for (corr_mat in seq_len(length(corr_list))) {
        omega_path <- huge::huge(corr_mat, nlambda = nlam, method = "glasso", verbose = FALSE)
        omega_hat <- ebic_selection(x = omega_path, n = n, s = corr_mat, param = param)

        auc_hat <- auc(truth = omega, huge_obj = omega_path)
        frobenius <- base::norm(omega_hat - omega, type = "F")

        adj_estimate <- abs(omega_hat) > 0
        tpr <- tpr(truth = omega, estimate = adj_estimate)
        fpr <- fpr(truth = omega, estimate = adj_estimate)
    }
}

identity <- function(x) {
    return(x)
}

nthroot <- function(x, n) {
    return(x^(1 / n))
}


generate_data <- function(n = 200, d = 50, g = identity) {
    t <- .15
    if (d == 50) {
        c <- -.2798173
    } else if (d == 250) {
        c <- -0.03072694
    } else if (d == 750) {
        c <- -0.01697741
    } else if (d == 1500) {
        c <- -0.01227969
    }

    prob_vec <- function(c) {
        return(
            rbinom(1,
                1,
                prob = (1 / sqrt(2 * pi)) * exp(euclid_norm(
                    (runif(2, min = 0, max = 1) - runif(2, min = 0, max = 1))
                ) / (2 * c))
            )
        )
    }
    ##### Initiate the precision matrix \Omega = \Sigma^{-1}
    Omega <- diag(nrow = d, ncol = d)

    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            Omega[i, j] <- Omega[j, i] <- t * prob_vec(c)
        }
    }

    if (!corpcor::is.positive.definite(Omega)) {
        Omega <- as.matrix(Matrix::nearPD(Omega, corr = T, keepDiag = T)$mat)
    }
    # diag(Omega) <- 1

    edge_number <- edgenumber(Precision = Omega)

    #### retrieve Sigma ####
    Sigma <- Matrix::chol2inv(chol(Omega)) # is simply a lot faster!175
    # Sigma <- solve(Omega)
    Sigma_corr <- stats::cov2cor(Sigma)
    #### Done: Now we can simulate from multivariate normal ####
    data_0 <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = Sigma_corr)

    data_1 <- apply(data_0, 2, g)

    return(list(
        "data_original" = data_0,
        "data_transformed" = data_1,
        "Omega" = Omega,
        "n_E" = edge_number
    ))
}

make_ordinal_general <- function(
    data = data,
    proportion = .5,
    namevector = c("binary" = T, "ordinal" = T, "poisson" = T),
    unbalanced = .2,
    low = .05,
    high = .1,
    lambda = 6,
    num_breaks = round(runif(1, 3, 10))) {
    # Define local vars:
    d <- ncol(data)
    d_1 <- floor(d * proportion)
    ordinal <- data[, 1:d_1]

    if (namevector[1] == T && namevector[2] == T && namevector[3] == T) {
        p_devide <- diff(c(0, floor(d_1 / sum(namevector)), d_1))
        ordinal_binary <- ordinal[, 1:p_devide[1]]
        ordinal_ordinal <- ordinal[, (p_devide[1] + 1):p_devide[2]]
        ordinal_poisson <- ordinal[, (p_devide[2] + 1):d_1]
    }
    if (namevector[1] == T && namevector[2] == T && namevector[3] == F) {
        p_devide <- diff(c(0, floor(d_1 / sum(namevector))))
        ordinal_binary <- ordinal[, 1:p_devide[1]]
        ordinal_ordinal <- ordinal[, (p_devide[1] + 1):d_1]
        ordinal_poisson <- NULL
    }
    if (namevector[1] == T && namevector[2] == F && namevector[3] == T) {
        p_devide <- diff(c(0, floor(d_1 / sum(namevector))))
        ordinal_binary <- ordinal[, 1:p_devide[1]]
        ordinal_ordinal <- NULL
        ordinal_poisson <- ordinal[, (p_devide[1] + 1):d_1]
    }
    if (namevector[1] == F && namevector[2] == T && namevector[3] == T) {
        p_devide <- diff(c(0, floor(d_1 / sum(namevector))))
        ordinal_binary <- NULL
        ordinal_ordinal <- ordinal[, 1:p_devide[1]]
        ordinal_poisson <- ordinal[, (p_devide[1] + 1):d_1]
    }
    if (namevector[1] == F && namevector[2] == F && namevector[3] == T) {
        ordinal_binary <- NULL
        ordinal_ordinal <- NULL
        ordinal_poisson <- ordinal[, 1:d_1]
    }
    if (namevector[1] == F && namevector[2] == T && namevector[3] == F) {
        ordinal_binary <- NULL
        ordinal_ordinal <- ordinal[, 1:d_1]
        ordinal_poisson <- NULL
    }
    if (namevector[1] == T && namevector[2] == F && namevector[3] == F) {
        ordinal_binary <- ordinal[, 1:d_1]
        ordinal_ordinal <- NULL
        ordinal_poisson <- NULL
    }

    if (namevector["binary"] == T) {
        #### Xbinary
        #### split binary so as to control fraction of unbalanced binary data
        pbin1 <- runif(floor(ncol(ordinal_binary) * (1 - unbalanced)), .4, .6)
        pbin2 <- runif((ncol(ordinal_binary) - floor(ncol(ordinal_binary) * (1 - unbalanced))), low, high)
        pbin <- c(pbin1, pbin2)
        for (i in seq_len(ncol(ordinal_binary))) {
            ordinal_binary[, i] <- qbinom(pnorm(scale(ordinal_binary[, i])), size = 1, prob = pbin[i])
        }
    }
    if (namevector["ordinal"] == T) {
        cum_mat <- list()
        for (k in seq_len(ncol(ordinal_binary))) {
            breaks <- runif(num_breaks)
            sum_breaks <- sum(breaks)
            breaks_norm <- sort(breaks) / sum_breaks
            cum_mat[[k]] <- c(0, cumsum(breaks_norm))
        }
        for (i in seq_len(ncol(ordinal_binary))) {
            u <- pnorm(scale(ordinal_ordinal[, i]))
            ordinal_ordinal[, i] <- cut(
                u,
                breaks = cum_mat[[i]],
                include.lowest = T,
                ordered_result = T,
                labels = 1:(length(cum_mat[[i]]) - 1)
            )
        }
    }
    if (namevector["poisson"] == T) {
        #### Poisson Variables from Threshold generation
        lambda_sample <- abs(rnorm(ncol(ordinal_poisson), mean = lambda, sd = 2))
        for (i in seq_len(ncol(ordinal_poisson))) {
            ordinal_poisson[, i] <- qpois(
                pnorm(
                    scale(ordinal_poisson[, i])
                ),
                lambda = lambda_sample[i]
            )
        }
    }

    ordinal <- cbind(ordinal_binary, ordinal_ordinal, ordinal_poisson)
    continuous <- (data[, (d_1 + 1):d])
    data_mixed <- as.data.frame(cbind(ordinal, continuous))
    ### declare ordinal variables as factors!
    for (f in seq_len(d_1)) {
        data_mixed[, f] <- factor(data_mixed[, f], ordered = T)
    }
    return(data_mixed)
}

hatR_fan <- function(data = data, verbose = T) {
    if (sum(sapply(data, is.factor)) == 0 && verbose == T) {
        cat("Warning, there are no factors in the input data.
            Did you declare ordinal variables as factors?")
    }
    d <- ncol(data)
    hatR <- diag(d)

    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            if (is.numeric(data[, i]) && is.numeric(data[, j])) {
                ### Fan et.al.
                hatR[i, j] <- hatR[j, i] <- sin(pcaPP::cor.fk(data[, i], data[, j]) * pi / 2)
                ###
            }
            if ((is.factor(data[, i]) && is.numeric(data[, j])) || (is.numeric(data[, i]) && is.factor(data[, j]))) {
                hatR[i, j] <- hatR[j, i] <- fan_case_2(data[, i], data[, j])
            }
            if (is.factor(data[, i]) && is.factor(data[, j])) {
                hatR[i, j] <- hatR[j, i] <- fan.case.3(data[, i], data[, j])
                ###
            }
        }
    }
    if (!corpcor::is.positive.definite(hatR)) {
        hatR_pd <- as.matrix(Matrix::nearPD(hatR, corr = T, keepDiag = T)$mat)
    } else {
        hatR_pd <- hatR
    }
    return(hatR_pd)
}


adhoc_polyserial <- function(cont, disc, maxcor = 0.9999) {
    disc <- as.integer(disc)
    n <- length(disc)
    cummarg_propotions <- c(0, cumsum(table(disc) / n))
    threshold_estimate <- qnorm(cummarg_propotions)

    values <- sort(as.integer(unique(disc)))
    lambda <- sum(dnorm(head(threshold_estimate, -1)[-1] * diff(values)))
    s_disc <- sd(disc)
    r <- npn_pearson(cont, disc)
    corr_hat <- r * s_disc / lambda

    if (abs(corr_hat) >= 1) {
        corr_hat <- sign(corr_hat) * maxcor
    }
    if (is.null(corr_hat)) {
        corr_hat <- polycor::polyserial(cont, disc)
    }
    return(corr_hat)
}

npn_pearson <- function(cont, disc) {
    f_x <- f_hat(cont)
    y <- as.integer(disc)
    r <- cor(f_x, y, method = "pearson")
    return(r)
}

f_hat <- function(x) {
    n <- length(x)
    npn_thresh <- 1 / (4 * (n^(1 / 4)) * sqrt(pi * log(n)))
    x <- qnorm(pmin(
        pmax(rank(x) / n, npn_thresh),
        1 - npn_thresh
    ))
    x <- x / sd(x)
    return(x)
}


hatR_poly <- function(data = data, verbose = T) {
    if (sum(sapply(data, is.factor)) == 0 && verbose == T) {
        cat("Warning, there are no factors in the input data.
            Did you declare ordinal variables as factors?")
    }
    d <- ncol(data)
    hatR <- matrix(1, d, d)

    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            if (is.numeric(data[, i]) && is.numeric(data[, j])) {
                hatR[i, j] <- hatR[j, i] <- 2 * sin(pi / 6 * spearman(data[, i], data[, j]))
            }
            if ((is.factor(data[, i]) && is.numeric(data[, j])) || (is.numeric(data[, i]) && is.factor(data[, j]))) {
                if (is.factor(data[, j])) {
                    hatR[i, j] <- hatR[j, i] <- adhoc_polyserial(data[, i], data[, j])
                } else {
                    hatR[i, j] <- hatR[j, i] <- adhoc_polyserial(data[, j], data[, i])
                }
            }
            if (is.factor(data[, i]) && is.factor(data[, j])) {
                hatR[i, j] <- hatR[j, i] <- polycor::polychor(data[, i], data[, j])
            }
        }
    }

    if (!requireNamespace("stringr", quietly = TRUE)) {
        stop("Please install package \"stringr\".")
    }
    pair <- hatR[lower.tri(hatR)]
    if (any(abs(pair) > .9)) {
        sapply(seq_along(hatR[lower.tri(hatR)][which(abs(pair) > .95)]), function(k) {
            warning(
                paste0(
                    "Correlation of the pair ",
                    stringr::str_c(as.character(which(hatR[lower.tri(hatR)][which(abs(pair) > .9)][k] == hatR,
                        arr.ind = T
                    )[, 1]), collapse = ",")
                ),
                " in the nonparanormal estimate is close to the boundary. Inverse might be misleading. "
            )
        })
    }
    if (!corpcor::is.positive.definite(hatR)) {
        hatR_pd <- as.matrix(Matrix::nearPD(hatR, corr = T, keepDiag = T)$mat)
    } else {
        hatR_pd <- hatR
    }
    return(hatR_pd)
}



hatR_mle <- function(data = data, verbose = T) {
    if (sum(sapply(data, is.factor)) == 0 && verbose == T) {
        cat("Warning, there are no factors in the input data.
            Did you declare ordinal variables as factors?")
    }
    d <- ncol(data)
    hatR <- matrix(1, d, d)

    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            if (is.numeric(data[, i]) && is.numeric(data[, j])) {
                hatR[i, j] <- hatR[j, i] <- cor(data[, i], data[, j])
            }
            if ((is.factor(data[, i]) && is.numeric(data[, j])) || (is.numeric(data[, i]) && is.factor(data[, j]))) {
                if (is.factor(data[, j])) {
                    hatR[i, j] <- hatR[j, i] <- polycor::polyserial(data[, i], data[, j])
                } else {
                    hatR[i, j] <- hatR[j, i] <- polycor::polyserial(data[, j], data[, i])
                }
            }
            if (is.factor(data[, i]) && is.factor(data[, j])) {
                hatR[i, j] <- hatR[j, i] <- polycor::polychor(data[, i], data[, j])
            }
        }
    }

    if (!requireNamespace("stringr", quietly = TRUE)) {
        stop("Please install package \"stringr\".")
    }
    pair <- hatR[lower.tri(hatR)]
    if (any(abs(pair) > .9)) {
        sapply(seq_along(hatR[lower.tri(hatR)][which(abs(pair) > .95)]), function(k) {
            warning(
                paste0(
                    "Correlation of the pair ",
                    stringr::str_c(as.character(which(hatR[lower.tri(hatR)][which(abs(pair) > .9)][k] == hatR,
                        arr.ind = T
                    )[, 1]), collapse = ",")
                ),
                " in the nonparanormal estimate is close to the boundary. Inverse might be misleading. "
            )
        })
    }
    if (!corpcor::is.positive.definite(hatR)) {
        hatR_pd <- as.matrix(Matrix::nearPD(hatR, corr = T, keepDiag = T)$mat)
    } else {
        hatR_pd <- hatR
    }
    return(hatR_pd)
}


spearman <- function(x, y) {
    rankx <- rank(x)
    ranky <- rank(y)
    rankmean <- (length(rankx) + 1) / 2

    rho <- (sum((rankx - rankmean) * (ranky - rankmean))) /
        (sqrt(sum((rankx - rankmean)^2) * sum((ranky - rankmean)^2)))
    return(rho)
}


#### input truth: True Precision Matrix
#### huge_obj: huge object e.g. glasso output
auc <- function(truth = truth, huge_obj = huge_obj) {
    stopifnot((class(huge_obj) == "huge"))

    adj_estimate <- lapply(
        seq_len(nlam),
        function(q) abs(huge_obj$icov[[q]]) > 0
    )

    fpr_lambda <- sapply(
        seq_len(nlam),
        function(l) {
            fpr(
                truth = truth,
                estimate = adj_estimate[[l]]
            )
        }
    )
    tpr_lambda <- sapply(
        seq_len(nlam),
        function(l) {
            tpr(
                truth = truth,
                estimate = adj_estimate[[l]]
            )
        }
    )

    fp <- c(fpr_lambda, 1)
    tp <- c(tpr_lambda, 1)
    dfp <- c(diff(fp), 0)
    dtp <- c(diff(tp), 0)

    AUC <- sum(tp * dfp) + sum(dtp * dfp) / 2
    return(AUC)
}

tpr <- function(truth = truth, estimate = estimate) {
    n_E <- edgenumber(truth)
    return(
        sum((truth != 0 & estimate != 0)[lower.tri(
            (truth != 0 & estimate != 0)
        )]) / n_E
    )
}

fpr <- function(truth = truth, estimate = estimate) {
    n_E <- edgenumber(truth)
    d <- dim(truth)[1]

    return(
        sum((truth == 0 & estimate != 0)[lower.tri(
            (truth == 0 & estimate != 0)
        )]) / (d * (d - 1) / 2 - n_E)
    )
}


ebic_selection <- function(x = x, param = param, n = n, s = s, partial = F) {
    stopifnot((class(x) == "huge"))
    d <- dim(x$data)[1]
    nlambda <- length(x$lambda)
    eBIC <- rep(0, nlambda)
    for (ind in seq_len(nlambda)) {
        huge_path <- x$path[[ind]]
        edge <- edgenumber(huge_path)
        huge_path[upper.tri(huge_path, diag = T)] <- 1
        zero_mat <- which(huge_path == 0, arr.ind = T)
        loglik <- suppressWarnings(glasso::glasso(s = s, rho = 0, nobs = n, zero = zero_mat)$loglik)
        eBIC[ind] <- -2 * loglik + edge * log(n) + 4 * edge * param * log(d)
    }

    Omega_hat <- x$icov[[which.min(eBIC)]]
    if (partial == T) {
        Omega_hat.standardized <- -1 * stats::cov2cor(Omega_hat)
        diag(Omega_hat.standardized) <- -1 * diag(Omega_hat.standardized)
    } else {
        Omega_hat.standardized <- stats::cov2cor(Omega_hat)
    }
    return(Omega_hat.standardized)
}
