data {
  int<lower=0> N;           // number of hospitals
  array[N] int<lower=0> y;  // remissions
  array[N] int<lower=0> K;  // number of trials
  vector<lower=0, upper=1>[N] w;
}
parameters {
  real mu;                  // population mean of success log-odds
  real<lower=0> sigma;      // population sd of success log-odds
  vector[N] alpha_raw;          // success log-odds
}
transformed parameters {
  vector[N] alpha_pcp = (1 - w) * mu + alpha_raw * sigma;
}
model {
  mu ~ normal(-1, 1);
  sigma ~ std_normal();
  alpha_raw ~ normal(w * mu, 1);
  y ~ binomial_logit(K, alpha_pcp);
}
generated quantities {
  vector[N] phi = inv_logit(alpha_pcp);
}