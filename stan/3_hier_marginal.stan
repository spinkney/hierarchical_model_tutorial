functions {
  real wishart_plus_lpdf (matrix X, real nu, matrix V_inv, real sigma_e, real sigma_t) {
    int P = rows(X);
    real log_det_Vinv = (1 - P) * log(sigma_e) - 0.5 * log(P * sigma_t^2 + sigma_e^2);
    return -0.5 * trace(V_inv * X) + 0.5 * nu * log_det_Vinv;
  }
}

data {
  int<lower=0> P;
  int<lower=0> J;
  vector[P] y_bar;
  matrix[P, P] X;
  real sigma_indv;
  real N;
}
parameters {
  real mu;
  real<lower=0> sigma_group;
//  real<lower=0> sigma_indv;
}
transformed parameters {
  // matrix[J, J] L;
      matrix[P, P] C;
  
  //   {
  //   //  matrix[J, J] C;
  //   real sigma_group_sq = square(sigma_group);
  //   real sigma_indv_sq = square(sigma_indv);
  //   real diag_elements = sigma_group_sq + sigma_indv_sq;
  // 
  //   C[1, 1] = diag_elements;
  // 
  //   for (i in 2:J) {
  //    C[i, i] = diag_elements;
  //     for (j in 1:i - 1) {
  //       C[i, j] = sigma_group_sq;
  //       C[j, i] = C[i, j];
  //     }
  //   }
  // 
  //   L = cholesky_decompose(C);
  // 
  // }
  
  {
    real sigma_group_sq = square(sigma_group);
    real sigma_indv_sq = square(sigma_indv);
    real xi = P * sigma_group_sq * sigma_indv_sq + sigma_indv_sq^2;
    real inv_diag_elements = ((P - 1) * sigma_group_sq + sigma_indv_sq) / xi;

    C[1, 1] = inv_diag_elements;

    for (i in 2:P) {
     C[i, i] = inv_diag_elements;
      for (j in 1:i - 1) {
        C[i, j] = -sigma_group_sq / xi;
        C[j, i] = C[i, j];
      }
    }

  }
  
}
model {
  sigma_group ~ normal(0, 5);
 // sigma_group ~ exponential(1);
  //L_est ~ wishart(J, C);
  X ~ wishart_plus(J - 1, C, sigma_indv, sigma_group);
  {y_bar} ~ multi_normal_prec(rep_vector(mu, P), C * J);
}
generated quantities {
//  matrix[J, J] Sigma = multiply_lower_tri_self_transpose(L);
}
