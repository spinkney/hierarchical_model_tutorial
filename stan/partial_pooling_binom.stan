data {
  int<lower=0> N;           // number of groups
  array[N] int<lower=0> y;  // binomial counts
  array[N] int<lower=0> K;  // number of trials
}
parameters {
  real mu;                  // population mean of success log-odds
  real<lower=0> sigma;      // population sd of success log-odds
  vector[N] alpha_std;      // success log-odds
}
transformed parameters {
  vector[N] alpha = mu + sigma * alpha_std;
}
model {
  mu ~ normal(-1, 1);
  sigma ~ std_normal();
  alpha_std ~ std_normal();
  y ~ binomial_logit(K, alpha);
}
generated quantities {
  vector[N] phi = inv_logit(alpha);
}