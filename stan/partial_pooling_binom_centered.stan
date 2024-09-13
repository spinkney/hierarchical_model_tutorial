data {
  int<lower=0> N;           // number of hospitals
  array[N] int<lower=0> y;  // remissions
  array[N] int<lower=0> K;  // number of trials
}
parameters {
  real mu;                  // population mean of success log-odds
  real<lower=0> sigma;      // population sd of success log-odds
  vector[N] alpha;          // success log-odds
}
model {
  mu ~ normal(-1, 1);
  sigma ~ std_normal();
  alpha ~ normal(mu, sigma);
  y ~ binomial_logit(K, alpha);
}
generated quantities {
  vector[N] phi = inv_logit(alpha);
}