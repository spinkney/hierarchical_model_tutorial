// Birnbaum, Z. and Saunders, S. (1969). A new family of life distributions. Journal of Applied Probability 6, 319327.
// https://www.smu.edu/-/media/Site/Cox/Faculty/Research/FoxEdward/BISA_J_Appl_Prob.pdf
// formula 4

functions {
//   real bs_lpdf (vector x, real beta, real alpha) {
//    int N = num_elements(x);
//    real a = N * (0.5 * log(beta) - log(2 * alpha) - 0.5 * (log2() + log(pi())));
// 
//    a += -beta * sum((1 - x / beta)^2 ./ (2 * alpha^2 * x));
// 
//    a += sum(log1p(x/beta) - (3.0/2) * log(x));
// 
//   return a;
// }
  real bs_lpdf (vector x, real beta, real alpha) {
   int N = num_elements(x);
   real a = N * (0.5 * log(beta) - log(2 * alpha) - 0.5 * (log2() + log(pi())));

   a += -beta * dot_product((1 - x / beta)^2 , 1 / (2 * alpha^2 * x));

   a += sum(log1p(x/beta) - (3.0/2) * log(x));

  return a;
}

}
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> beta;
  real<lower=0> alpha;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  beta ~ exponential(1);
  alpha ~ exponential(1);
  y ~ bs(beta, alpha);
}

