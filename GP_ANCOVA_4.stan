## GP_ANCOVA_4
## THis Stan file seeks to use zero-sum kernel to include categorical effect into the GP.
## yij = beta * X + f(x) + muj + eij. 
## f(x) ~ N(0, K)
## WHERE K = Kzero-sum * Ksquared-exp

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
  vector[K] beta;               // linear effect coefs.
  vector[J] mu;               // group effect
  vector[N] eta;            // latent GP at all inputs
}

transformed parameters {
  matrix[N, N] K_base = gp_exp_quad_cov(X, alpha, rho);
  vector[N] f;
  
  matrix[N, N] K_zero_sum;
  for (i in 1:N) {
    K_zero_sum[i, i] = 1;
    for (j in (i+1):N){
      if (group[i] == group[j]) {
        K_zero_sum[i, j] = 1;
      } else {
        K_zero_sum[i, j] = -1/(J-1);
      }
      K_zero_sum[j, i] = K_zero_sum[i, j];
      }
  }
  
  matrix[N, N] K_f = K_zero_sum .* K_base;
  
  for (i in 1:N) {
    K_f[i, i] += 1e-6; // jitter for stability
  }
  
  matrix[N, N] L_Kernel;  // decomposition
  L_Kernel = cholesky_decompose(K_f);
  f = L_Kernel * eta; // this is f(x)
  // f = f - mean(f); // suspicious: constrain f(x) to make it identifiable, or it would meddle with muj
  vector[N] F;
  vector[N] xbeta;
  
  for(i in 1: N){
    xbeta[i] = dot_product(X[i], beta);
  }
  F = mu[group] + f + xbeta ; // this is the function value with group effect.
}

model {
  // Priors
  alpha ~ student_t(1, 0, 1);
  rho ~ inv_gamma(5, 5);
  sigma ~ cauchy(0, 0.1);
  mu ~ normal(0, 2);
  eta ~ normal(0,1);
  beta ~ normal(0,1);
  // Likelihood
  Y1 ~ normal(F[1:N_obs], sigma);
}
