data {
  int<lower=0> N;        // number of obs
  vector[N] y;           // study effect
  vector[N] study_sigma; // study variablity
  int<lower=0> D;        // number of subjects
  int<lower=0> G;        // number of grade levels
  array[N] int<lower=1, upper=D> subject_idx;
  array[N] int<lower=1, upper=G> grade_idx;
}
parameters {
  real intercept;               // pop. effect 
  real<lower=0> sigma;
  vector[N] mu_study;    // varying intercept by study
}
model {
  intercept ~ std_normal(); // Prior model
  sigma ~ exponential(1);
  mu_study ~ normal(intercept, sigma); // Centered hierarchical model
  y ~ normal(mu_study, study_sigma); // Observational model
}
