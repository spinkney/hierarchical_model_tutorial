library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)

mod_funnel <- cmdstan_model("./stan/basic_funnel.stan")
mod_funnel_out <- mod_funnel$sample(parallel_chains = 4)

np_funnel <- nuts_params(mod_funnel_out)
funnel_draws <- mod_funnel_out$draws(c("phi", "tau"), format = "array") |>
  merge_chains() |>
  as_draws_df()

divergence_plot <- mcmc_scatter(
  funnel_draws, 
  pars = c("phi", "tau"), 
  np = np_funnel, 
  size = 0.2
)

divergence_plot

# compare to stans hessian
mod_funnel_out$init_model_methods(hessian = TRUE)

calculate_hessian <- function(y, x) {
  hessian_matrix <- 
    matrix(c((x^2 * exp(-y)) / 2 - 1/9,  x * exp(-y), 
              x * exp(-y)             , -1 / exp(y)), 
            nrow = 2, byrow = TRUE)
  return(hessian_matrix)
}

calculate_hessian(0.5, 1.2)
mod_funnel_out$hessian(c(0.5, 1.2))

# calculate the condition number
# see what values have large condition numbers

mod_funnel_repar <- cmdstan_model("./stan/basic_funnel_repar.stan")
mod_funnel_repar_out <- mod_funnel_repar$sample(parallel_chains = 4)

# what do you think the hessian will look like?

mod_funnel_repar_out$init_model_methods(hessian = TRUE)
mod_funnel_repar_out$hessian(c(0.5, 1.2))

np_funnel_repar <- nuts_params(mod_funnel_repar_out)
funnel_repar_draws <- mod_funnel_repar_out$draws(c("phi_raw", "tau_raw"), format = "array") |>
  merge_chains() |>
  as_draws_df()

divergence_plot_repar <- mcmc_scatter(
  funnel_repar_draws, 
  pars = c("phi_raw", "tau_raw"), 
  np = np_funnel_repar, 
  size = 0.2
)

divergence_plot_repar
