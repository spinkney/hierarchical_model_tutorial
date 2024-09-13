parameters {
  real tau;
  real phi;
}
model {
  tau ~ normal(0, 3);
  phi ~ normal(0, exp(0.5 * tau));
}