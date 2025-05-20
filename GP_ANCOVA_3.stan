// This file model Yij = fj(x) + muj + beta*X + f(x) + e
// fj(x) + muj is the grouping effect. the effect could be varying on differnet x's
// beta*X + f(x) is the base effect (the effect excluded from grouping effect?? or the base effect for all the groups to add upon...)

data {
  int<lower=1> N_obs;
  int<lower=1> K;             // number of covariates
  int<lower=1> J;             // number of groups
  array[N_obs] vector[K] X1;        // training inputs
  vector[N_obs] Y1;           // training outputs
  int<lower=1, upper=J> group1[N_obs]; // training group indices

  int<lower=1> N_test;
  array[N_test] vector[K] X2;        
  int<lower=1, upper=J> group2[N_test]; // test group indices
}

transformed data{
  int<lower = 1> N = N_obs + N_test;  // total number of observations
  array[N] vector[K] X;               // total input. Store as "array" because gp_exp_quad_cov needs X to be an array..
  for (i in 1:N_obs) {
    X[i] = X1[i];
  }
  for (j in 1:N_test) {
    X[N_obs + j] = X2[j];
  }
  
  int<lower=1> group[N];               // total group membership identifier. 
  for (i in 1:N_obs) {
    group[i] = group1[i];
  }
  for (j in 1:N_test) {
    group[N_obs + j] = group2[j];
  }

  array[J, N] int idx_j;             // storing indices per group. We used N as length because stan doesn't allow dynamic indexing, and this seems to be a common trick, we actually only used first Nj (number of observations per group) spaces per vector, others are just peddings. intentionally left to be blank
                                     // it tells the location of the 1st, 2nd, 3rd, ....X'th group observations in the total dataset, 
  array[J] int Nj;                   // Nj is the number of observations per group (here for each group I didn't differentiate if the ob is training or testing)
  for (j in 1:J) {                   // because we know the location of them in the total dataset. eventually they will be mapped back to the reuslt vector of length of N. Since we knew which ones in the N-length vecotr are training or testing, so it doesn't matter at all
    Nj[j] = 0;                       // initialize all Nj to be 0
    } 
    // Fill in indices
   for (i in 1:N) {
     int g = group[i];              // extract the "group membership"" of the i'th observation
     Nj[g] += 1;                    // plus one for the counter for that specific group "g"
     idx_j[g][Nj[g]] = i;           // Store the index i for group g at the Nj[g]-th position,
                                    // Nj[g] represents the number of times we've encountered group g so far.
  }
}

parameters {
  real<lower=0> alpha;          // GP marginal std (non-group specific)
  real<lower=0> rho;            // GP length-scale (non-group specific)
  real<lower=0> sigma;          // observation noise (non-group specific)
  vector[K] beta;               // linear effect coefs (non-group specific)
  vector[J] mu_a;                 // group effect (group-specific, mean)
  vector<lower=0>[J] g_alpha;   // group effect (group-specific, kernel)
  vector<lower=0>[J] g_rho;     // group effect (group-specific, kernel)
  
  vector[N] eta;                    // eta ~ N(0,1), for randomness
  array[J] vector[max(Nj)] g_eta;   // g_eta ~ N(0,1), for randomness
}

transformed parameters {
  // Yij = f(x) + beta*X + muj + fj(x) + error
  // F = f(x) + beta*X + muj + fj(x)
  // grand trend: f(x) + beta*X
  // group specific trend: muj + fj(x)
  
  vector[J] mu = 3 * mu_a;  //this is non-centered parameterization. leads to mu ~ N(0, 3^2)
  
  // this is the block for grand trend part i.e. f(x) + beta*X
  vector[N] f;         // define f(x), the non-linear part
  matrix[N, N] K_f = gp_exp_quad_cov(X, alpha, rho); // define kernel for f(x), non-linear part
  
  for (i in 1:N) {
    K_f[i, i] += 1e-6; // jitter for stability
  }

  matrix[N, N] L_Kernel;  // C decomposition. seems tricky when sample size gets larger.
  L_Kernel = cholesky_decompose(K_f);
  f = L_Kernel * eta; // this is f(x) before constraint
  f = f - mean(f); // constrain f(x) to make it identifiable, or it would meddle with muj
  
  vector[N] xbeta;     // define beta * X, linear part
  for(i in 1: N){      // dot product, returning a singular. NOTICE: different from ".*", which is elementvise multiplication, returning a vector. GOOGLE IT.
    xbeta[i] = dot_product(X[i], beta);
  }                    // the reason why I HAVE TO do this is, I used "array" as input for X (see data declaration blcok).
                       // I use "array" for X because I wanna use "gp_exp_quad_cov", assuming it might be the best practice for efficency to generate KERNEL MATRIX in STAN
  
  // this is the block for group-specific trend part i.e. muj + fj(x)
  vector[N] f_j;
  
  for (j in 1:J) {
    // Covariance matrix
    // Prepare group X 
    int Nj_j = Nj[j];               // temporarily save the present number of group observations
    array[Nj_j] vector[K] Xj;      // subset Xj for group j. only inputs relevent to j would be considered.
    for (n in 1: Nj_j){             // extract the observation for group j to a subset Xj.
      Xj[n] = X[idx_j[j,n]];        // We must do this one extra loop. WE cannot use the brilliant X[idx_j[j,]] loopping from 1 to J. Why, because we have 0 on the padding area. Stan would read it and take it fucking seriously.
    }
    matrix[Nj_j, Nj_j] K_j = gp_exp_quad_cov(Xj, g_alpha[j], g_rho[j])
                             + 1e-6 * diag_matrix(rep_vector(1.0, Nj_j));

    // Cholesky factorization for GP prior
    matrix[Nj_j, Nj_j] L_Kj = cholesky_decompose(K_j);

    // prepare temporary group specific eta vector. why this? same reason with Xj!  L_Kj * g_eta[j]; won't be done because of padding 0's!!!!!
    vector[Nj_j] gg_eta; // extract the g_eta for group g. So gg_eta
    for (n in 1: Nj_j){
      gg_eta[n] = g_eta[j,n];
    }
    vector[Nj_j] f_jj = L_Kj * gg_eta;

    // Store f_jj into f_j vector. MAPPING back to the position.
    for (i in 1:Nj_j) {
      int n = idx_j[j, i];
      f_j[n] = f_jj[i];
    }
  }
  
  // Total function
  vector[N] F;
  F =  mu[group] + f_j + f + xbeta ; // this is the function value with group effect.
}
// mu[group] performs indexing of mu by position according to the group membership of each observation.

model {
  // Priors

  alpha ~ cauchy(0, 1);
  rho ~ inv_gamma(5, 5);
  eta ~ normal(0,1);
  beta ~ normal(0,1);
  sigma ~ normal(0, 1);

  mu_a ~ normal(0, 1); //auxiliary std.normal for mu
  g_alpha ~ cauchy(0, 1);
  g_rho ~ inv_gamma(5, 5);
  for (j in 1:J) {
    g_eta[j] ~ normal(0, 1);  // standard normal for GP prior
  }
  // Likelihood
  Y1 ~ normal(F[1:N_obs], sigma);
}



