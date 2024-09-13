parameters {
  real tau_raw;
  real phi_raw;
}
transformed parameters {
  real tau = 3 * tau_raw;
  real phi = exp(tau / 2) * phi_raw;
}
model {
  tau_raw ~ std_normal(); 
  phi_raw ~ std_normal();
}