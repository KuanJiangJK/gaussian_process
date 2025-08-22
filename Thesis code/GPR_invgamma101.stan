## GP File from the Stan userbook

functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1,
                     array[] real x1,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    
    matrix[N1, 1] X1_mat;  // transform to matrix. 
      for (i in 1:N1){
        X1_mat[i,1] = x1[i];
      }
  
    matrix[N2, 1] X2_mat;  // transform to matrix
      for (i in 1:N2){
        X2_mat[i,1] = x2[i];
      }

    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N2, N2] k_x2_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      K = gp_exp_quad_cov(x1, alpha, rho);
      for (n in 1:N1) {
        K[n, n] = K[n, n] + square(sigma);
      }
      K += X1_mat*X1_mat';
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha, rho);
      k_x1_x2 += X1_mat*X2_mat';
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      k_x2_x2 = gp_exp_quad_cov(x2, alpha, rho);
      k_x2_x2 += X2_mat*X2_mat';
      cov_f2 = k_x2_x2 - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}
data {
  int<lower=1> N1;
  array[N1] real x1;
  vector[N1] y1;
  int<lower=1> N2;
  array[N2] real x2;
}
transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-9;

    matrix[N1, 1] X1_mat;  // transform to matrix. 
      for (i in 1:N1){
        X1_mat[i,1] = x1[i];
      }
  
    matrix[N2, 1] X2_mat;  // transform to matrix
      for (i in 1:N2){
        X2_mat[i,1] = x2[i];
      }
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N1, N1] L_K;
  {
    matrix[N1, N1] K = gp_exp_quad_cov(x1, alpha, rho);
    real sq_sigma = square(sigma);
    K += X1_mat*X1_mat';
    // diagonal elements
    for (n1 in 1:N1) {
      K[n1, n1] = K[n1, n1] + sq_sigma;
    }

    L_K = cholesky_decompose(K);
  }

  rho ~ inv_gamma(5, 5);
  alpha^2 ~ inv_gamma(1, 0.1);
  sigma ~ std_normal();

  y1 ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;

  f2 = gp_pred_rng(x2, y1, x1, alpha, rho, sigma, delta);
  for (n2 in 1:N2) {
    y2[n2] = normal_rng(f2[n2], sigma);
  }
}
