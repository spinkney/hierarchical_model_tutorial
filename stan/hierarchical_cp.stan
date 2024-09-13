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
}
parameters {
  real mu; // Population location
  real<lower=0> tau; // Population scale
  vector[K] theta; // Centered individual parameters
}
model {
  mu ~ normal(0, 5); // Prior model
  tau ~ normal(0, 5); // Prior model
  theta ~ normal(mu, tau); // Centered hierarchical model
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}
generated quantities {
  array[N] real y_post_pred = normal_rng(theta[indiv_idx], sigma);
}