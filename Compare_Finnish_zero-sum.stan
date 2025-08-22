# comparsion with or without zero-sum kernle
# with zero-sum kernel
data {
  int<lower=1> J;             // number of groups/categories
  int<lower=1> K;             // number of dimensions for GP
  
  // training
  int<lower=1> N1;
  array[N1] vector[K] X1;        // training inputs
  int<lower=1, upper=J> group1[N1]; // training group indices
  vector[N1] Y1;
  
  // testing
  int<lower=1> N2;
  array[N2] vector[K] X2;        // testing inputs
  int<lower=1, upper=J> group2[N2]; // test group indices
}

transformed data {
  vector[N1] mu_a = rep_vector(0, N1); // to be used in the marginalized thingy
  real delta = 1e-5; // to be added as a jitter for stability
  matrix[N1, K] X1_mat;
  for (i in 1:N1){
      X1_mat[i] = X1[i]';
      }
  
  matrix[N2, K] X2_mat;
  for (i in 1:N2){
      X2_mat[i] = X2[i]';
      }
}

parameters {
  real<lower=0> lambda;          // GP std (non-group)
  real<lower=0> c_global2;          // GP std (non-group)
  real<lower=0> rho;            // GP length-scale (non-group)
  vector[J] mu_j;                 // group effect (group-specific)
  vector<lower=0>[J] g_lambda;   // group effect (group-specific, kernel)
  vector<lower=0>[J] c_group2;   // group effect (group-specific, kernel)
  vector<lower=0>[J] g_rho;     // group effect (group-specific, kernel)
  real<lower=0> sigma;          // observation noise (non-group specific)
  vector[K] beta;               // linear effect coefs (non-group specific). this should be integrated out!!
}

transformed parameters{
  real<lower=0> alpha;
  alpha = sqrt((c_global2*(lambda^2)) / (c_global2 + lambda^2));
  vector<lower=0>[J] g_alpha;
      for (j in 1:J) {
      g_alpha[j] = sqrt((c_group2[j]*(g_lambda[j]^2)) / (c_group2[j] + g_lambda[j]^2));
    }
}
model {
  matrix[N1, N1] L_K; // cholesky composition for the training dataset
  
  // calculate total kernel K_X1_X1 = K_base + K_group * Kzerosum
  {
    matrix[N1, N1] K_X1_X1 = cov_exp_quad(X1, alpha, rho); // base kernel
    
    matrix[N1, N1] K_group;
    for (i in 1:N1) {
      for (j in i:N1) {
        if (group1[i] == group1[j]) {
          int g = group1[i];
          real d = sqrt(dot_self(X1[i] - X1[j]));
          real k_val = square(g_alpha[g]) * exp(-0.5 * square(d / g_rho[g]));
          K_group[i, j] = k_val;
          K_group[j, i] = k_val; // ensure symmetry
        } else {
          K_group[i, j] = 0;
          K_group[j, i] = 0;
        }
      }
    }    
    
    matrix[N1, N1] K_zero_sum; //zerosum kernel
    for (i in 1:N1) {
      K_zero_sum[i, i] = 1;
      for (j in (i+1):N1){
        if (group1[i] == group1[j]) {
          K_zero_sum[i, j] = 1;
        } else {
          K_zero_sum[i, j] = -1/(J-1);
        }
          K_zero_sum[j, i] = K_zero_sum[i, j];
        }
      }
  
    real sq_sigma = square(sigma); // gaussian noise
    
    K_X1_X1  +=  K_group .* K_zero_sum + X1_mat*diag_matrix(rep_vector(1, K))*X1_mat';
    // K_X1_X1  +=  K_group + X1_mat*diag_matrix(rep_vector(1, K))*X1_mat';

    // diagonal elements
    for (n1 in 1:N1) {
      K_X1_X1[n1, n1] = K_X1_X1[n1, n1] + sq_sigma + delta;
    }
    L_K = cholesky_decompose(K_X1_X1);
  }
  
  // mean function part: only mu_j, beta integrated out
  vector[N1] mu = mu_a + mu_j[group1];
  
  // priors
  lambda ~ cauchy(0, 1);
  rho ~ inv_gamma(5, 5);
  g_lambda ~ cauchy(0, 1);
  c_group2 ~ inv_gamma(2,1);
  c_global2 ~ inv_gamma(2,1);
  g_rho ~ inv_gamma(5, 5);
  mu_j ~ std_normal();
  sigma ~ std_normal();
  beta ~ std_normal();
  Y1 ~ multi_normal_cholesky(mu, L_K);
}
