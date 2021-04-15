### just checking sth. out

margal_prop_r <- .2
margal_prop_r_1 <- .1
gamma_r <- qnorm(margal_prop_r)
gamma_r_1 <- qnorm(margal_prop_r_1)

rho <- -.3


X <- rnorm(100, 4, 4)

X_scale <- (X - mean(X))/sd(X)

ell <- vector(length = 100)

for (i in 1:length(X_scale)) {

gamma_hat_r <- (gamma_r - rho*X_scale[i])/sqrt(1-rho^2)
gamma_hat_r_1 <- (gamma_r_1 - rho*X_scale[i])/sqrt(1-rho^2)

p_xi_given_x_k <- pnorm(gamma_hat_r) - pnorm(gamma_hat_r_1) 

numerator <- (1-rho^2)^(3/2)
dell_p_xi_given_x_k <- numerator*(dnorm(gamma_hat_r)*(gamma_r*rho-X_scale[i]) - dnorm(gamma_hat_r_1)*(gamma_r_1*rho-X_scale[i]))
ell[i] <- dell_p_xi_given_x_k/p_xi_given_x_k

}

ell
sum(ell)
