library(cmdstanr)
library(posterior)
library(bayesplot)

## simulate data
mod_normal_hier <- cmdstan_model("./stan/hierarchical_sim.stan")

lme4::lmer(stan_data$y ~ 1 + (1 | stan_data$indiv_idx))

K <- 9
N_per_indiv <- c(10, 5, 1000, 10, 1, 5, 100, 10, 5)
indiv_idx <- as.factor(rep(1:K, N_per_indiv))
N <- length(indiv_idx)
sigma <- 10


mod_normal_hier_data <- mod_normal_hier$sample(
  data = list(
    N = N,
    K = K,
    indiv_idx = indiv_idx,
    sigma = sigma
  ),
  iter_sampling = 1,
  chains = 1,
  seed = 38243293,
  fixed_param = TRUE
)

y <- mod_normal_hier_data$draws("y", format = "matrix") |>
  as.vector()

stan_data <- list(
  y = y,
  N = N,
  K = K,
  indiv_idx = indiv_idx,
  sigma = sigma
)

max_N <- max(N_per_indiv)
D <- matrix(0, nrow = max_N, K)
mu_bar <- sapply(1:K, function(i) mean(stan_data$y[indiv_idx == i]))

for (k in 1:K) {
  D[1:N_per_indiv[k], k] <- stan_data$y[indiv_idx == k] 
  # if (N_per_indiv[k] != max_N) {
  #   N_samples <- max_N - N_per_indiv[k]
  #   D[(N_per_indiv[k] + 1):max_N, k] <- mu_bar[k] +
  #     rnorm(N_samples, mean = 0, sd = 0.001)
  # }
}

y.bar <- matrix(rowMeans(D), max_N, 1) 
# column vector of ones 
ones.vec <- matrix(1, K, 1) 
# sample scatter matrix (Equation 11) 
#D <- D + matrix(rnorm(max_N * K, mean = 0, sd = 1), nr = max_N, K)
S <- (D - y.bar %*% t(ones.vec)) %*% t((D - y.bar %*% t(ones.vec))) 

eigen_S <- eigen(S)
M <- length(which(eigen_S$values > 0))
lambda <- eigen_S$values[which(eigen_S$values > 0)]
Q <- eigen_S$vectors[, 1:M]

Q %*% diag(lambda) %*% t(Q)
LaplacesDemon::as.inverse(S)

#S <- cov(t(D))

out <- fastmatrix::ldl(S)

out$d[which(is.na(out$d))] <- 0.01
out$d[which(out$d <= 0)] <- 0.01
out$lower[is.infinite(out$lower)] <- 0.01
L_est <- diag(sqrt(out$d)) %*% out$lower
tcrossprod(L_est)

mod_marginal <- cmdstan_model("./stan/3_hier_marginal.stan")

mod_marginal_out <- mod_marginal$sample(
  data = list(
    P = 1000,
    J = K,
    y_bar = c(y.bar),
    X = Q %*% diag(lambda) %*% t(Q),
    sigma_indv = sigma,
    N = 40
  ),
  parallel_chains = 4
)
mod_marginal_out$summary(c("mu", "sigma_group"))

sqrt(10^2 + 3.5^2)

mod_hier_ncp <- cmdstan_model("./stan/hierarchical_ncp.stan")
mod_normal_hier_ncp <- mod_hier_ncp$sample(
  data = list(
    y = y,
    N = N,
    K = K,
    indiv_idx = indiv_idx,
    sigma = sigma
  ),
  seed = 12312312,
  parallel_chains = 4,
  adapt_delta = 0.95
)
mod_normal_hier_ncp$summary(c("mu", "tau", "theta"))



sigma_e <- 0.5
sigma_t <- 1.23
P <- 10

G <- matrix(sigma_t^2, P, P)
diag(G) <- sigma_e^2 + sigma_t^2
1/sqrt(det(G))

sqrt(det(solve(G)))

eta <- (sigma_e * sigma_t) / sqrt(P * sigma_t^2 + sigma_e^2)
eta / (sigma_t * sigma_e^P)

sigma_e^(1-P) / sqrt(P * sigma_t^2 + sigma_e^2)

