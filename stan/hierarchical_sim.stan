// the Stan code is from https://betanalpha.github.io/assets/case_studies/hierarchical_modeling.html#3_The_Fundamental_Degeneracies_of_Normal_Hierarchical_Models
// The code is copyrighted by Michael Betancourt and licensed under the new BSD (3-clause) license

data {
  int<lower=1> N; // Number of observations
  int<lower=0> K; // Number of individuals in hierarchy
  array[N] int<lower=1, upper=K> indiv_idx; 
  real<lower=0> sigma; 
}
transformed data {
  real mu = 4.5;           // Population location
  real<lower=0> tau = 3.5; // Population scale
}
generated quantities {
  array[K] real theta = normal_rng(rep_vector(mu, K), tau);
  array[N] real y = normal_rng(theta[indiv_idx], sigma);
}