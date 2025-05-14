data {
  int<lower=1> N_obs;
  int<lower=1> K;             // number of covariates
  int<lower=1> J;             // number of groups
  array[N_obs] vector[K] X1;        // training inputs
  vector[N_obs] Y1;           // training outputs
  int<lower=1, upper=J> group1[N_obs]; // training group indices

  int<lower=1> N_test;
  array[N_test] vector[K] X2;        // 
  int<lower=1, upper=J> group2[N_test]; // test group indices
}

transformed data{
  int<lower = 1> N = N_obs + N_test;
  array[N] vector[K] X;
  for (i in 1:N_obs) {
    X[i] = X1[i];
  }
  for (j in 1:N_test) {
    X[N_obs + j] = X2[j];
  }
  
  int<lower=1> group[N];
  for (i in 1:N_obs) {
    group[i] = group1[i];
  }
  for (j in 1:N_test) {
    group[N_obs + j] = group2[j];
  }
}

parameters {
  real<lower=0> alpha;        // GP marginal std
  real<lower=0> rho;          // GP length-scale
  real<lower=0> sigma;        // observation noise
  vector[J] mu;               // group effect
  vector[N] eta;            // latent GP at all inputs
}

transformed parameters {
  matrix[N, N] K_f = gp_exp_quad_cov(X, alpha, rho);
  vector[N] f;
  
  for (i in 1:N) {
    K_f[i, i] += 1e-6; // jitter for stability
  }
  
  matrix[N, N] L_Kernel;  // decomposition
  L_Kernel = cholesky_decompose(K_f);
  f = L_Kernel * eta;
  f = f - mean(f);
  vector[N] f_group;
  f_group = mu[group] + f; // this is the function value with group effect.
}

model {
  // Priors
  alpha ~ normal(0, 5);
  rho ~ inv_gamma(5, 5);
  sigma ~ cauchy(0, 0.1);
  mu ~ normal(0, 2);
  eta ~ normal(0,1);
  // Likelihood
  Y1 ~ normal(f_group[1:N_obs], sigma);
}


// generated quantities {
//   vector[N_test] f_test;
//   vector[N_test] y_test_pred;
// 
//   {
//     matrix[N_obs, N_obs] K = K_f;
//     matrix[N_test, N_obs] K_s;
//     matrix[N_test, N_test] K_ss;
// 
//     for (i in 1:N_test) {
//       for (j in 1:N_obs) {
//         real sq_dist = squared_distance(X2[i], X1[j]);
//         K_s[i, j] = alpha^2 * exp(-0.5 * sq_dist / rho^2);
//       }
//     }
// 
//     for (i in 1:N_test) {
//       for (j in 1:N_test) {
//         real sq_dist = squared_distance(X2[i], X2[j]);
//         K_ss[i, j] = alpha^2 * exp(-0.5 * sq_dist / rho^2);
//       }
//       K_ss[i, i] += 1e-6; // jitter
//     }
// 
//     matrix[N_test, N_obs] K_s_Kinv = K_s * inverse(K);
//     f_test = K_s_Kinv * f;
//     for (i in 1:N_test) {
//       y_test_pred[i] = normal_rng(f_test[i] + mu[group2[i]], sigma);
//     }
//   }
// }


