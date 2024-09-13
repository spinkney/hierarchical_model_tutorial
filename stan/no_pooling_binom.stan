data {
  int<lower=0> N;           // number of hospitals
  array[N] int<lower=0> y;  // remissions
  array[N] int<lower=0> K;  // number of trials
}
parameters {
  vector[N] alpha;
}
model {
  y ~ binomial_logit(K, alpha);
}
generated quantities {
  vector[N] phi = inv_logit(alpha);
}