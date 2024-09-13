data {
  int<lower=0> N; // Number of observations
  vector[N] y; // Observations
  real<lower=0> sigma; // Measurement variability
  int<lower=0> K;
  array[N] int<lower=1, upper=K> indiv_idx;
  vector<lower=0, upper=1>[K] w;
}
transformed data {
 // real tau = 3;
}
parameters {
  real mu;            // Population location
  real<lower=0> tau;  // Population scale
  vector[K] eta;      // Non-centered individual parameters
}
transformed parameters {
 vector[K] tau_wt = tau - w * (tau - 1);
 vector[K] theta = (eta * tau + mu * w) ./ tau_wt;
}
model {
  mu ~ normal(0, 5);                   // Prior model
  tau ~ normal(0, 5);                  // Prior model
  eta ~ normal(mu * (1 - w), tau_wt);
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}
generated quantities {
  array[N] real y_post_pred = normal_rng(theta[indiv_idx], sigma);
  real theta_p = normal_rng(mu, tau);
}