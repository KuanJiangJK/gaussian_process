functions {           // This function samples from the posterior predictive distribution at new input locations,
                      // conditioned on observed data
  vector gp_pred_rng(real[] x2,  //Stan syntax for a 1D array (not a vector!) of type real
                     vector y1,
                     real[] x1, //** what does [] stand for?
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = rows(y1); //** number of training data point. difference between rows and size function?
    int N2 = size(x2); //** number of training data point. difference between rows and size function?
    vector[N2] f2; //** to store the sampled predictive function values (but I don't want, can I drop?)
    {
      matrix[N1, N1] L_K; 
      vector[N1] K_div_y1;  
      matrix[N1, N2] k_x1_x2; 
      matrix[N1, N2] v_pred; 
      vector[N2] f2_mu; 
      matrix[N2, N2] cov_f2; // cov for testing point
      matrix[N2, N2] diag_delta; // I 
      matrix[N1, N1] K; // cov + sigma^2*I for training point
      K = cov_exp_quad(x1, alpha, rho);
      for (n in 1:N1)
        K[n, n] = K[n,n] + square(sigma);
      L_K = cholesky_decompose(K); // the following 3 lines are calculating inverse of K* y1. this is L_k which enables L_K*L_K' = K
                                    // here L_K*L_K' = K, and what we want is K^-1*y1 
      K_div_y1 = mdivide_left_tri_low(L_K, y1); // left division of y1 by a lower-triangular view of L_K. inverse(tri(L_K)) * y1
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)'; //
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho); // covariance between training and testing
      f2_mu = (k_x1_x2' * K_div_y1); // conditional distribution: https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Computational_methods
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2); //
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta); 
    }
    return f2;
  }
}
data {
  int<lower=1> N1;
  real x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  real x2[N2];
}
transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N1, N1] L_K;
  {
    matrix[N1, N1] K = cov_exp_quad(x1, alpha, rho);
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:N1)
      K[n1, n1] = K[n1, n1] + sq_sigma;

    L_K = cholesky_decompose(K);
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  y1 ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;

  f2 = gp_pred_rng(x2, y1, x1, alpha, rho, sigma, delta);
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f2[n2], sigma);
}
