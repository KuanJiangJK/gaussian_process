functions {
  // Deterministic ANCOVA mean: mu_n = dot(beta, X[n]) + mu_j[group[n]]
  vector ancova_mean(array[] vector X, int[] group, vector beta, vector mu_j) {
    int N = size(X);
    vector[N] mu;
    for (n in 1:N)
      mu[n] = dot_product(beta, X[n]) + mu_j[group[n]];
    return mu;
  }

  // Posterior predictive RNG for ANCOVA (only usable in generated quantities)
  vector ancova_pred_rng(array[] vector X, int[] group, vector beta, vector mu_j, real sigma) {
    int N = size(X);
    vector[N] y;
    
    for (n in 1:N)
      y[n] = dot_product(beta, X[n]) + mu_j[group[n]];
    return y;
  }
}

data {
  int<lower=1> J;                 // number of groups
  int<lower=1> K;                 // number of covariates

  // training
  int<lower=1> N1;
  array[N1] vector[K] X1;         // training covariates
  int<lower=1, upper=J> group1[N1];
  vector[N1] Y1;

  // testing (for prediction/plots)
  int<lower=1> N2;
  array[N2] vector[K] X2;         // test covariates
  int<lower=1, upper=J> group2[N2];
}

parameters {
  vector[K] beta;                 // common slope(s)
  vector[J] mu_j;                 // group intercepts
  real<lower=0> sigma;            // residual SD
}

model {
  vector[N1] mu1;

  // Priors (weakly-informative; adjust as needed)
  beta  ~ normal(0, 1);           // common slope(s)
  mu_j  ~ normal(0, 1);           // group intercepts
  sigma ~ normal(0, 1);           // half-normal via <lower=0>

  // Likelihood
  mu1 = ancova_mean(X1, group1, beta, mu_j);
  Y1  ~ normal(mu1, sigma);
}

generated quantities {
  // Training fitted values and PPC
  vector[N1] mu1_hat = ancova_mean(X1, group1, beta, mu_j);
  vector[N1] y1_rep;
  for (n in 1:N1)
    y1_rep[n] = normal_rng(mu1_hat[n], sigma);

  // Test predictive mean and posterior predictive draws
  vector[N2] mu2 = ancova_mean(X2, group2, beta, mu_j);
  vector[N2] y2  = ancova_pred_rng(X2, group2, beta, mu_j, sigma);
}
