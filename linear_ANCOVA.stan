data {
  int<lower=1> N_obs;
  int<lower=0> K; // number of predictors (including intercept)
  int<lower=0> J; // number of groups
  matrix[N_obs, K] X1; // covariates for training points (including 1's in the first column)
  vector[N_obs] Y1; // training point
  int<lower=1, upper=J> group1[N_obs]; // group index for each observation
  int<lower=1> N_test;
  matrix[N_test, K] X2; // testing point
  int<lower=1, upper=J> group2[N_test];
}

// transformed data {
//   int<lower=1> N = N_obs + N_test;
//   vector[N] group = group1 + group2
// }

parameters {
  vector[K] beta; // from beta0 to betak (K predictors including intercept)
  real<lower=0> sigma; // variance for error. ANCOVA assume same for all groups.
  vector[J] mu; // from mu1 to muj. This is the constant for each group
}

// transformed parameters{
//   
//   vector[N] F;
//   F = X*beta + mu[group]
// }

model {
  beta ~ normal(0, 1);
  mu ~ normal(0, 1);
  sigma ~ cauchy(0, 2);
  // for (n in 1:N_test){
  //   target += normal_lpdf(Y1[n] | X1[n] * beta + mu[group1[n]], sigma);
  // }
  target += normal_lpdf(Y1 | X1 * beta + mu[group1], sigma);
}

generated quantities {
  vector[N_test] y_test_pred;
  for (n in 1:N_test) {
    y_test_pred[n] = normal_rng(X2[n] * beta + mu[group2[n]], sigma);
    }
}

