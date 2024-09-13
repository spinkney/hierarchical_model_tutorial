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
transformed data {
  // real sd_y = sd(y);
  // real sd_x = sd(X);
  // vector[N] X_std = (X - mean(X)) / sd_x;
}
parameters {
  real intercept;
  real intercept_sigma_study;
  real sigma_subject_raw;
  
  real beta_mu;
  real beta_sigma;
  
  vector[D] mu_subject_raw;
  vector[N] mu_study_raw;       
}
transformed parameters {
 real sigma_subject = exp(sigma_subject_raw);
 vector[N] sigma_study = exp(intercept_sigma_study + beta_sigma * X);
 vector[D] mu_subject = mu_subject_raw * sigma_subject;
 vector[N] mu_study = intercept + mu_subject[subject_idx] + beta_mu * X + mu_study_raw .* sigma_study;
}
model {
  intercept ~ std_normal(); 
  intercept_sigma_study ~ skew_normal(-1, 1, -2); 
  sigma_subject_raw ~ skew_normal(-1, 1, -4); 
  
  beta_mu ~ std_normal();
  beta_sigma ~ skew_normal(0, 1, -2); 

  mu_study_raw ~ std_normal();
  mu_subject_raw ~ std_normal();
  
  y ~ normal(mu_study, study_sigma); // Observational model
}
generated quantities {
  real mu_est = mean(mu_study);
  real intercept_est = mu_est - mean(beta_mu * X) - mean(mu_subject);
  real intercept_sigma_est = log(sd(mu_study) / mean(exp(beta_sigma * X)));
}
