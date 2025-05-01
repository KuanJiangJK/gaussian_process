data {
  int<lower=1> N;             // total number of points
  int<lower=1> J;             // number of groups
  vector[N] x;                // single input predictor
  vector[N] y;                // output
  int group[N]; // group indicator
}

transformed data {
  int nj[J];
  for (j in 1:J) {
    nj[j] = 0;
    }
  for (i in 1:N) {
    nj[group[i]] += 1; // count number of samples per group
    }
}

parameters {
  matrix[J,2] beta;
  vector<lower=0>[J] alpha2;     // GP magnitude per group
  vector<lower=0>[J] rho;        // length scale per group
  vector<lower=0>[J] sigma;      // noise std per group
}

model {
  // Priors
  for (j in 1:J) {
    beta[j] ~ normal(0, 5);    
    alpha2[j] ~ normal(0, 1);  
    sigma[j] ~ normal(0, 1);   
    rho[j] ~ normal(0, 1);     
  }
  // Likelihood
  for (j in 1:J){
    int nj_current = nj[j]; // the number of observations of the current group
    int idx[nj_current];  // index array for group j. it tells the corresponding location of each observation from a certain group, so that data won't be messed up. For a neatly arranged dataset, this line might be unnecessary.
    {
      int k = 1;
        for (i in 1:N) {  //go through the data, see if the point belongs to the current group
          if (group[i] == j) {
            idx[k] = i;
            k += 1;
            }
        }
    }
      // Build mean vector for group j
      vector[nj_current] mu_j; // mean function, a vector with size of nj
      for (k in 1:nj_current) {
        mu_j[k] = beta[j,1] + beta[j,2] * x[idx[k]];
      }

      // Build covariance matrix for group j
      matrix[nj_current, nj_current] K_j;
      for (k1 in 1:nj_current) {
        for (k2 in 1:nj_current) {
          real sq_dist = square(x[idx[k1]] - x[idx[k2]]);
          K_j[k1, k2] = alpha2[j] * exp(-0.5 * sq_dist / square(rho[j]));
        }
        K_j[k1, k1] += square(sigma[j]) + 1e-6;  // add noise variance + jitter on diagonal
      }

      // Group-specific multivariate normal likelihood
      target += multi_normal_lpdf(y[idx] | mu_j, K_j);
    }
  }


