data {
  int<lower=1> N_obs;      // number of training points
  vector[N_obs] X1;        // 1 predictor (not including intercept)
  vector[N_obs] Y1;        // training outcomes
  int<lower=1> N_test;     // number of testing points
  vector[N_test] X2;       // test predictors
}

transformed data {
  real delta = 1e-9;       // jitter to stabilize cholesky_decompose
  int<lower=1> N = N_obs + N_test; // sum up all data points
  vector[N] X;             // combined input for training + testing
  for (i in 1:N_obs) {
    X[i] = X1[i];
  }
  for (j in 1:N_test) {
    X[N_obs + j] = X2[j];
  }
}

parameters {
  real alpha;              // intercept
  real beta;               // slope for X
  real<lower=0> sigmaf;     // GP magnitude
  real<lower=0> sigman;     // noise std
  vector[N] eta;            // standard normal for GP
}

transformed parameters {
  matrix[N, N] Kernel;
  matrix[N, N] L_Kernel;  // this is the kernel after cholesky_decompose
  vector[N] f; // non-linear part
  vector[N] F; // whole function. F = beta*x + alpha + f
  real<lower=0> sigmaf2 = square(sigmaf);
  real<lower=0> l = 1.0;     // fixed length-scale. could also be learned from bayesian hierarchical modelling. Here for simplicity

  // Build squared exponential kernel
  for (i in 1:(N - 1)) {
    Kernel[i, i] = sigmaf2 + delta;
    for (j in (i + 1):N) {
      real distance = X[i] - X[j];
      Kernel[i, j] = sigmaf2 * exp(-0.5 * square(distance) / square(l));
      Kernel[j, i] = Kernel[i, j];
    }
  }
  Kernel[N, N] = sigmaf2 + delta;

  L_Kernel = cholesky_decompose(Kernel);

  f = L_Kernel * eta;
  F = alpha + beta * X + f;
}

model {
  // Priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigmaf ~ student_t(1, 0, 1);
  sigman ~ normal(0, 1);
  eta ~ std_normal();

  // Likelihood
  Y1 ~ normal(F[1:N_obs], sigman);
}
