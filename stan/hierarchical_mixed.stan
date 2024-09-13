// the Stan code is from https://betanalpha.github.io/assets/case_studies/hierarchical_modeling.html#3_The_Fundamental_Degeneracies_of_Normal_Hierarchical_Models
// The code is copyrighted by Michael Betancourt and licensed under the new BSD (3-clause) license
data {
  int<lower=0> N; // Number of observations
  vector[N] y; // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Number of individual contexts in hierarchy
  int<lower=0> K;
  // Individual context from which each observation is generated
  array[N] int<lower=1, upper=K> indiv_idx;
  
  int<lower=0, upper=K> K_ncp; // Number of noncentered individuals
  array[K_ncp] int<lower=1, upper=K> ncp_idx; // Index of noncentered individuals
  
  int<lower=0, upper=K> K_cp; // Number of centered individuals
  array[K_cp] int<lower=1, upper=K> cp_idx; // Index of noncentered individuals
}
parameters {
  real mu; // Population location
  real<lower=0> tau; // Population scale
  vector[K_ncp] theta_ncp; // Non-centered individual parameters
  vector[K_cp] theta_cp; // Ccentered individual parameters
}
transformed parameters {
  // Recentered individual parameters
  vector[K] theta;
  theta[ncp_idx] = mu + tau * theta_ncp;
  theta[cp_idx] = theta_cp;
}
model {
  mu ~ normal(0, 5); // Prior model
  tau ~ normal(0, 5); // Prior model
  
  theta_ncp ~ normal(0, 1); // Non-centered hierarchical model
  theta_cp ~ normal(mu, tau); // Centered hierarchical model
  
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}