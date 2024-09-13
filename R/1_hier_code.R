library(cmdstanr)
library(posterior)
library(bayesplot)

## simulate data
mod_normal_hier <- cmdstan_model("./stan/hierarchical_sim.stan")

K <- 9
N_per_indiv <- c(10, 5, 1000, 10, 1, 5, 100, 10, 5)
indiv_idx <- rep(1:K, N_per_indiv)
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

## fit model centered
mod_hier_cp <- cmdstan_model("./stan/hierarchical_cp.stan")
mod_normal_hier_cp <- mod_hier_cp$sample(
  data = list(
    y = y,
    N = N,
    K = K,
    indiv_idx = indiv_idx,
    sigma = sigma
  ),
  seed = 12312312,
  parallel_chains = 4
)

np_cp <- nuts_params(mod_normal_hier_cp)
cp_draws <- mod_normal_hier_cp$draws(c("theta", "tau"), format = "array") 

## fit model non-centered
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
  parallel_chains = 4
)

np_ncp <- nuts_params(mod_normal_hier_ncp)
ncp_draws <- mod_normal_hier_ncp$draws(c("eta", "tau"), format = "array") 

## fit model partially-centered
mod_hier_pcp <- cmdstan_model("./stan/hierarchical_pcp.stan")
I_y <- by(y, indiv_idx, FUN = sd)
var_var <- var(I_y, na.rm = T)
# when w = 0 => centered
# when w = 1 => non-centered
# w <- (1 + ((N_per_indiv)/I_y)^2)^-1
w <-  (1 + (I_y^2 * N_per_indiv) / (var_var * sum(N_per_indiv)))^-1
w[is.na(w)] <- 1

# w[ncp_idx] <- 0
# w[cp_idx] <- 1

mod_normal_hier_pcp <- mod_hier_pcp$sample(
  data = list(
    y = y,
    N = N,
    K = K,
    indiv_idx = indiv_idx,
    sigma = sigma,
    w = as.vector(w)
  ),
  seed = 12312312,
 #iter_sampling = 5000,
  parallel_chains = 4
)

np_pcp <- nuts_params(mod_normal_hier_pcp)
pcp_draws <- mod_normal_hier_pcp$draws(c("eta", "tau"), format = "array") 

## fit model mixed-centered
mod_hier_mixed <- cmdstan_model("./stan/hierarchical_mixed.stan")
data = list(
  y = y,
  N = N,
  K = K,
  indiv_idx = indiv_idx,
  sigma = sigma)
thresh <- 25
ncp_idx <- which(table(data$indiv_idx) <= thresh)
cp_idx <- which(table(data$indiv_idx) > thresh)

data$K_ncp <- length(ncp_idx)
data$ncp_idx <- array(ncp_idx)

data$K_cp <- length(cp_idx)
data$cp_idx <- array(cp_idx)

mod_normal_hier_mixed  <- mod_hier_mixed $sample(
  data = data,
  seed = 12312312,
 # iter_sampling = 5000,
  parallel_chains = 4
)


## inspect diagnostics
mod_normal_hier_cp$sampler_diagnostics()
table(mod_normal_hier_cp$sampler_diagnostics()[,,6])
table(mod_normal_hier_cp$sampler_diagnostics()[,,5])

table(mod_normal_hier_ncp$sampler_diagnostics()[,,6])
table(mod_normal_hier_ncp$sampler_diagnostics()[,,5])

table(mod_normal_hier_mixed$sampler_diagnostics()[,,6])
table(mod_normal_hier_mixed $sampler_diagnostics()[,,5])

table(mod_normal_hier_pcp$sampler_diagnostics()[,,6])
table(mod_normal_hier_pcp$sampler_diagnostics()[,,5])

plot_scatter <- function(k, par, draws, nuts_params) {
  mcmc_scatter(
    draws, 
    pars = c(paste0(par, "[", k, "]", sep=""), "tau"), 
    transform = list(tau = "log"), 
    np = nuts_params, 
    size = 0.2
  )
}

plot_list_cp <- lapply(1:K, function(x) plot_scatter(x, "theta", cp_draws, np_cp))
plot_list_ncp <- lapply(1:K, function(x) plot_scatter(x, "eta", ncp_draws, np_ncp))
plot_list_pcp <- lapply(1:K, function(x) plot_scatter(x, "eta", pcp_draws, np_pcp))

k <- 6
bayesplot_grid(
  plots = list(plot_list_cp[[k]], plot_list_ncp[[k]], plot_list_pcp[[k]]),
  legends = FALSE
  )

## Workshop Exercise

# Verify that the pcp gives the same prior distribution as the centered and non-centered
# Bonus: do this using simulation based calibration

# Update the pcp model to estimate the weights. Did the model fit well or not? If not, can you diagnose why?
