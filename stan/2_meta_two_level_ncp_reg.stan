data {
  int<lower=0> N;        // number of obs
  vector[N] y;           // study effect
  vector[N] study_sigma; // study variablity
  int<lower=0> D;        // number of subjects
  int<lower=0> G;        // number of grade levels
  array[N] int<lower=1, upper=D> subject_idx;
  array[N] int<lower=1, upper=G> grade_idx;
  vector[N] X;
}
parameters {
  real intercept;              
  real intercept_sigma;
  real beta_mu;
  real beta_sigma;
  vector[N] mu_study_raw;       // varying intercept by study
}
transformed parameters {
  vector[N] sigma_study = exp(intercept_sigma + beta_sigma * X);
  vector[N] mu_study = intercept + beta_mu * X + mu_study_raw .* sigma_study;
}
model {
  intercept ~ std_normal(); 
  // see distribution at https://rok-cesnovar.github.io/stan-distributions/skew_normal
  intercept_sigma ~ skew_normal(-1, 1, -3); 
  beta_mu ~ std_normal();
  beta_sigma ~ skew_normal(0, 1, -2); 

  mu_study_raw ~ std_normal();

  y ~ normal(mu_study, study_sigma); // Observational model
}
generated quantities {
  real mu_est = mean(mu_study);
  real intercept_est = mu_est - mean(beta_mu * X);
  real intercept_sigma_est = log(sd(mu_study) / mean(exp(beta_sigma * X)));
}