data {
  int<lower=0> N;           // number of hospitals
  array[N] int<lower=0> y;  // remissions
  array[N] int<lower=0> K;  // number of trials
  
  int<lower=0> J; // number of covariates
  matrix[N, J] X; // covariate matrix
}
parameters {
  real alpha;
  vector[J] beta;
}
model {
  y ~ binomial_logit_glm(K, X, alpha, beta);
}
generated quantities {
  vector[N] y_pred = inv_logit(alpha + X * beta);
}