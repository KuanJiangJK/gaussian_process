data {
  int<lower=1> N_obs;
  int<lower=0> K; // number of predictors (including intercept, corresponding to 1's in X)
  matrix[N_obs, K] X1;
  vector[N_obs] Y1; // training point
  int<lower=1> N_test;
  matrix[N_test, K] X2; // testing point
}

transformed data {
  real delta = 1e-9;  // jitter to stabilize the cholesky_decompose
  int<lower=1> N = N_obs + N_test;
  matrix[N, K] X; //combine observation and test points
  for (i in 1:N_obs) {
    X[i,] = X1[i,];
  }
  for (j in 1:N_test) {
    X[N_obs + j,] = X2[j,];
  }
  matrix[N, K-1] XX = block(X, 1, 2, N, K-1); // here drop the 1st column (all 1's), might be helpful for calculating distance between two points
}

parameters {
  vector[K] beta;
  real<lower=0> sigmaf; // magnitude of the GP, the sigma2 before the exp...
  vector<lower=0>[K-1] rec_L2; // reciprocal of the square length scale. 1/(L^2)
  real<lower=0> sigman; // standard deviation of measurement error
  vector[N] eta; // this, together with L_kernel (cholesky decomposition), is a substitute of the sampled function f, making sampling more efficient. f = Lk*eta. When L_kernel is determined given parameters (sigmaf & length scale), eta~N(0,1) is the randomness lefted for a random function.
}

transformed parameters {
  matrix[K-1, K-1] M = diag_matrix(rec_L2); // each dimension has its own "length scale".
  vector[N] f;  // value of GP(0, KERNEL)
  vector[N] F;  // value of total function, Beta*X + GP = Yi - e
  matrix[N, N] Kernel;
  matrix[N, N] L_Kernel;  // decomposition
  real<lower=0> sigmaf2 = square(sigmaf);
  
  // Covariance matrix with squared exponential kernel
  for (i in 1:(N - 1)) {
    Kernel[i, i] = sigmaf2 + delta;
    for (j in (i + 1):N) {
      vector[K-1] distance = to_vector(XX[i,] - XX[j,]);
      Kernel[i, j] = sigmaf2* exp(-0.5 * distance' * M * distance);
      Kernel[j, i] = Kernel[i, j];
    }
  }
  Kernel[N, N] = sigmaf2 + delta;
  L_Kernel = cholesky_decompose(Kernel);
  f = L_Kernel * eta;
  F = f + X*beta;
  }

model {
  rec_L2 ~ multi_normal(rep_vector(0, K-1), diag_matrix(rep_vector(1, K-1)));
  beta ~ multi_normal(rep_vector(0, K), diag_matrix(rep_vector(1, K)));
  sigmaf ~ student_t(1, 0, 1);
  eta ~ multi_normal(rep_vector(0, N), diag_matrix(rep_vector(1, N))); // previously i used normal distribution here, what is the difference here? how does stan deal with vector & just a single value? When I place a normal instead of multivariate normal, what happened?
  sigman ~ normal(0, 1);
  Y1 ~ normal(F[1:N_obs], sigman); // this is the likelihood of data.
}




