// function for predicted function value
functions {
  // Start ########################### FUNCTION 1 #######################################
  vector gp_pred_rng(array[] vector X2, // testing point
                     vector Y1, // training output
                     int J, // number of groups
                     int K, // number of input dimensions
                     array[] vector X1, // training
                     int[] group1, // training group
                     int[] group2, // testing group
                     real alpha, // GP magnitude
                     real rho, // Length-scale
                     real sigma, // gaussian noise
                     real delta, // jitter
                     // vector g_alpha, // group GP magnitude (vector whose length is the same with the number of group)
                     // vector g_rho, // group GP length scale
                     vector mu_j) { // group constant shift
    int N1 = rows(Y1); 
    int N2 = size(X2); 
    matrix[N1, K] X1_mat;  // transform to matrix. 
    // I set input of testing and training point as array to use the stan function cov_exp_quad but if i want to calculate X'%*%X (which would be needed after marginalizing out the lienar coef Beta), I need matrix form.
      for (i in 1:N1){
        X1_mat[i] = X1[i]';
  }
  
  matrix[N2, K] X2_mat;  // transform to matrix
  for (i in 1:N2){
    X2_mat[i] = X2[i]';
      }
      
    vector[N2] f2; // testing point function values
    {
      matrix[N1, N1] L_K; 
      vector[N1] Y1_centered;
      vector[N1] K_div_Y1; // component for posterior predictive distribution of P(f*|Y)
      matrix[N1, N2] v_pred; // component for posterior predictive distribution of P(f*|Y)
      matrix[N1, N1] K_X1_X1; // Cov within training
      matrix[N1, N2] K_X1_X2; // Cov between training and testing
      matrix[N2, N2] K_X2_X2; // Cov within testing
      vector[N2] f2_mu; // mean function of P(f*|Y)
      matrix[N2, N2] cov_f2; // Cov of P(f*|Y)
      matrix[N2, N2] diag_delta; // jitter for stability
      
      Y1_centered = Y1 - mu_j[group1];
  // calculate total kernel K_x1_x1 and K^-1*y1 ############# 
  {
    K_X1_X1 = cov_exp_quad(X1, alpha, rho); // base kernel (grand non-linear trend)
    
    // matrix[N1, N1] K_group; // group specifc kernel, each element calculated based on its corresponding group kernel parameters
    // for (i in 1:N1) {
    //   for (j in i:N1) {
    //     if (group1[i] == group1[j]) {
    //       int g = group1[i];
    //       real d = sqrt(dot_self(X1[i] - X1[j]));
    //       real k_val = square(g_alpha[g]) * exp(-0.5 * square(d / g_rho[g]));
    //       K_group[i, j] = k_val;
    //       K_group[j, i] = k_val; // ensure symmetry
    //     } else {
    //       K_group[i, j] = 0;
    //       K_group[j, i] = 0;
    //     }
    //   }
    // }
    
    // matrix[N1, N1] K_zero_sum; //zerosum kernel, to be multiplied by the group specific kernel
    // for (i in 1:N1) {
    //   K_zero_sum[i, i] = 1;
    //   for (j in (i+1):N1){
    //     if (group1[i] == group1[j]) {
    //       K_zero_sum[i, j] = 1;
    //     } else {
    //       K_zero_sum[i, j] = -1.0 / (J - 1.0);
    //     }
    //       K_zero_sum[j, i] = K_zero_sum[i, j];
    //     }
    //   }

    // K_X1_X1  +=  K_group .* K_zero_sum + X1_mat*X1_mat'; // add them up. X1_mat*X1_mat' comes from marginalization over beta. originally should be [X1_mat %*% Cov_rho %*% X1_mat']. But drop the cov since I use prior beta ~ N(0,1). 
    K_X1_X1  +=  X1_mat*X1_mat'; // add them up. X1_mat*X1_mat' comes from marginalization over beta. originally should be [X1_mat %*% Cov_rho %*% X1_mat']. But drop the cov since I use prior beta ~ N(0,1). 
    // diagonal elements, add the gaussian noise and add the jitter (delta, to be declared in transformed data, as it was done in the user guide)
    real sq_sigma = square(sigma); // gaussian noise
    for (n1 in 1:N1) {
      K_X1_X1[n1, n1] = K_X1_X1[n1, n1] + sq_sigma + delta;
    }
    L_K = cholesky_decompose(K_X1_X1);
  }

  // end calculating total kernel K_x1_x1 #########################
  
  // calcule total kernel K_x1_x2 and predictive distribution variance component #########################
  {
    K_X1_X2 = cov_exp_quad(X1, X2, alpha, rho); // base kernel, grand non-linear
    
    // matrix[N1, N2] K_group_X1_X2; // group specific kernel
    // for (i in 1:N1) {
    //   for (j in 1:N2) { // this one is not a square matrix N*N
    //     if (group1[i] == group2[j]) {
    //       int g = group1[i];
    //       real d = sqrt(dot_self(X1[i] - X2[j]));
    //       real k_val = square(g_alpha[g]) * exp(-0.5 * square(d / g_rho[g]));
    //       K_group_X1_X2[i, j] = k_val;
    //     } else {
    //       K_group_X1_X2[i, j] = 0;
    //     }
    //   }
    // }
    // 
    // matrix[N1, N2] K_zero_sum_X1_X2; //zerosum kernel
    // for (i in 1:N1) {
    //   for (j in 1:N2){
    //     if (group1[i] == group2[j]) {
    //       K_zero_sum_X1_X2[i, j] = 1;
    //     } else {
    //       K_zero_sum_X1_X2[i, j] = -1.0 /(J - 1.0);
    //     }
    //   }
    // }
    // K_X1_X2  += K_group_X1_X2 .* K_zero_sum_X1_X2 + X1_mat*X2_mat';

    K_X1_X2  +=  X1_mat*X2_mat';

  }

 // END calculating total kernel K_x1_x2 #########################
 
 // calculate total kernel K_x2_x2 to be used in cov_f2 ############### Initiall I was wrong this part, forgetting combine grand non-linear kernel with group-specific thingy, as well as zero-sum, thus resulting in a large number of Non-positive definite matrix warnings, effectiveless sampling.
   {
    K_X2_X2 = cov_exp_quad(X2, alpha, rho); // base kernel
    
    // matrix[N2, N2] K_group; // group kernel, each element generated based on group kernel parameters
    // for (i in 1:N2) {
    //   for (j in i:N2) {
    //     if (group2[i] == group2[j]) {
    //       int g = group2[i];
    //       real d = sqrt(dot_self(X2[i] - X2[j]));
    //       real k_val = square(g_alpha[g]) * exp(-0.5 * square(d / g_rho[g]));
    //       K_group[i, j] = k_val;
    //       K_group[j, i] = k_val; // ensure symmetry
    //     } else {
    //       K_group[i, j] = 0;
    //       K_group[j, i] = 0;
    //     }
    //   }
    // }
    
    // matrix[N2, N2] K_zero_sum; //zerosum kernel
    // for (i in 1:N2) {
    //   K_zero_sum[i, i] = 1;
    //   for (j in (i+1):N2){
    //     if (group2[i] == group2[j]) {
    //       K_zero_sum[i, j] = 1;
    //     } else {
    //       K_zero_sum[i, j] = -1.0 / (J - 1.0);
    //     }
    //       K_zero_sum[j, i] = K_zero_sum[i, j];
    //     }
    //   }

    // K_X2_X2  +=  K_group .* K_zero_sum + X2_mat*X2_mat'; // X2_mat*X2_mat' comes from marginalization over beta
    K_X2_X2  += X2_mat*X2_mat'; // X2_mat*X2_mat' comes from marginalization over beta

  }
  // END calculating total kernel K_x2_x2 #########################

      K_div_Y1 = mdivide_left_tri_low(L_K, Y1_centered); // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
      K_div_Y1 = mdivide_right_tri_low(K_div_Y1', L_K)'; // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
      f2_mu = (K_X1_X2' * K_div_Y1); // mean of P(f*|Y)
      v_pred = mdivide_left_tri_low(L_K, K_X1_X2); // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
      cov_f2 = K_X2_X2 - v_pred' * v_pred; // cov of P(f*|Y)
      
      f2_mu += mu_j[group2];
      
      diag_delta = diag_matrix(rep_vector(delta, N2)); // add the jitter to K_X2_X2
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta); 
    }
    return f2;
  }

// END ########################### FUNCTION 1 #######################################
// Start ########################### FUNCTION 2 #######################################

  vector gp_glob_pred_rng(array[] vector X2, // testing point
                     vector Y1, // training output
                     int J, // number of groups
                     int K, // number of input dimensions
                     array[] vector X1, // training
                     int[] group1, // training group
                     int[] group2, // testing group
                     real alpha, // GP magnitude
                     real rho, // Length-scale
                     real sigma,
                     real delta,
                     vector mu_j// jitter
                     ) { // group constant shift
    int N1 = rows(Y1); 
    int N2 = size(X2); 
    matrix[N1, K] X1_mat;  // transform to matrix. 
                           // I set input of testing and training point as array to use the stan function cov_exp_quad but if i want to calculate X'%*%X (which would be needed after marginalizing out the lienar coef Beta), I need matrix form.
for (i in 1:N1){
  X1_mat[i] = X1[i]';
      }
  
    matrix[N2, K] X2_mat;  // transform to matrix
      for (i in 1:N2){
        X2_mat[i] = X2[i]';
}

vector[N1] Y1_centered;
Y1_centered = Y1 - mu_j[group1];

vector[N2] f2glob; // testing point function values
{
  matrix[N1, N1] L_K; 
  vector[N1] K_div_Y1; // component for posterior predictive distribution of P(f*|Y)
  matrix[N1, N2] v_pred; // component for posterior predictive distribution of P(f*|Y)
  matrix[N1, N1] K_X1_X1; // Cov within training
  matrix[N1, N2] K_X1_X2; // Cov between training and testing
  matrix[N2, N2] K_X2_X2; // Cov within testing
  vector[N2] f2_mu; // mean function of P(f*|Y)
  matrix[N2, N2] cov_f2; // Cov of P(f*|Y)
  matrix[N2, N2] diag_delta; // jitter for stability
  
  // calculate total kernel K_x1_x1 and K^-1*y1 ############# 
  {
    K_X1_X1 = cov_exp_quad(X1, alpha, rho); // base kernel (grand non-linear trend)
    K_X1_X1  +=  X1_mat*X1_mat';  

    // diagonal elements, add the gaussian noise and add the jitter (delta, to be declared in transformed data, as it was done in the user guide)
    real sq_sigma = square(sigma); // gaussian noise
    for (n1 in 1:N1) {
      K_X1_X1[n1, n1] = K_X1_X1[n1, n1] + sq_sigma + delta;
    }
    L_K = cholesky_decompose(K_X1_X1);
  }

  // end calculating total kernel K_x1_x1 #########################
  
  // calcule total kernel K_x1_x2 and predictive distribution variance component #########################
  {
    K_X1_X2 = cov_exp_quad(X1, X2, alpha, rho); // base kernel, grand non-linear
    K_X1_X2  += X1_mat*X2_mat';
    
  }
  
  // END calculating total kernel K_x1_x2 #########################
  
  // calculate total kernel K_x2_x2 to be used in cov_f2 ############### Initiall I was wrong this part, forgetting combine grand non-linear kernel with group-specific thingy, as well as zero-sum, thus resulting in a large number of Non-positive definite matrix warnings, effectiveless sampling.
  {
    K_X2_X2 = cov_exp_quad(X2, alpha, rho); // base kernel
    K_X2_X2  +=  X2_mat*X2_mat'; // X2_mat*X2_mat' comes from marginalization over beta
  }
  // END calculating total kernel K_x2_x2 #########################
  
  K_div_Y1 = mdivide_left_tri_low(L_K, Y1_centered); // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
  K_div_Y1 = mdivide_right_tri_low(K_div_Y1', L_K)'; // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
                                   f2_mu = (K_X1_X2' * K_div_Y1); // mean of P(f*|Y)
  v_pred = mdivide_left_tri_low(L_K, K_X1_X2); // this step is to compute components of the posterior predictive distribution P(f*|Y), same thing was done as well in the stan user guide
  cov_f2 = K_X2_X2 - v_pred' * v_pred; // cov of P(f*|Y)
  diag_delta = diag_matrix(rep_vector(delta, N2)); // add the jitter to K_X2_X2
  f2glob = multi_normal_rng(f2_mu, cov_f2 + diag_delta); 
}
return f2glob;
  }
}

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
  real<lower=0> alpha;          // GP std (non-group)
  real<lower=0> rho;            // GP length-scale (non-group)
  vector[J] mu_j;                 // group effect (group-specific)
  // vector<lower=0>[J] g_alpha;   // group effect (group-specific, kernel)
  // vector<lower=0>[J] g_rho;     // group effect (group-specific, kernel)
  real<lower=0> sigma;          // observation noise (non-group specific)
  vector[K] beta;               // linear effect coefs (non-group specific). this should be integrated out!!
}

model {
  matrix[N1, N1] L_K; // cholesky composition for the training dataset
  
  // calculate total kernel K_X1_X1 = K_base + K_group * Kzerosum
  {
    matrix[N1, N1] K_X1_X1 = cov_exp_quad(X1, alpha, rho); // base kernel
    
    // matrix[N1, N1] K_group;
    // for (i in 1:N1) {
    //   for (j in i:N1) {
    //     if (group1[i] == group1[j]) {
    //       int g = group1[i];
    //       real d = sqrt(dot_self(X1[i] - X1[j]));
    //       real k_val = square(g_alpha[g]) * exp(-0.5 * square(d / g_rho[g]));
    //       K_group[i, j] = k_val;
    //       K_group[j, i] = k_val; // ensure symmetry
    //     } else {
    //       K_group[i, j] = 0;
    //       K_group[j, i] = 0;
    //     }
    //   }
    // }    
    
    // matrix[N1, N1] K_zero_sum; //zerosum kernel
    // for (i in 1:N1) {
    //   K_zero_sum[i, i] = 1;
    //   for (j in (i+1):N1){
    //     if (group1[i] == group1[j]) {
    //       K_zero_sum[i, j] = 1;
    //     } else {
    //       K_zero_sum[i, j] = -1/(J-1);
    //     }
    //     K_zero_sum[j, i] = K_zero_sum[i, j];
    //   }
    // }
    
    real sq_sigma = square(sigma); // gaussian noise
    
    // K_X1_X1  +=  K_group .* K_zero_sum + X1_mat*diag_matrix(rep_vector(1, K))*X1_mat';
    K_X1_X1  +=  X1_mat*diag_matrix(rep_vector(1, K))*X1_mat';

    // diagonal elements
    for (n1 in 1:N1) {
      K_X1_X1[n1, n1] = K_X1_X1[n1, n1] + sq_sigma + delta;
    }
    L_K = cholesky_decompose(K_X1_X1);
  }
  
  // mean function part: only mu_j, beta integrated out
  vector[N1] mu = mu_a + mu_j[group1];
  
  // priors
  alpha ~ cauchy(0, 1);
  rho ~ inv_gamma(5, 5);
  // g_alpha ~ cauchy(0, 1);
  // g_rho ~ inv_gamma(5, 5);
  mu_j ~ std_normal();
  sigma ~ std_normal();
  beta ~ std_normal();
  Y1 ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;
  vector[N2] f2glob;
  // f2glob = gp_glob_pred_rng(X2, Y1, J, K, X1, group1, group2, alpha, rho, sigma, delta);
  // f2 = gp_pred_rng(X2, Y1, J, K, X1, group1, group2, alpha, rho, sigma, delta, g_alpha, g_rho, mu_j);
  f2glob = gp_glob_pred_rng(X2, Y1, J, K, X1, group1, group2, alpha, rho, sigma, delta, mu_j);
  f2 = gp_pred_rng(X2, Y1, J, K, X1, group1, group2, alpha, rho, sigma, delta, mu_j);

  for (n2 in 1:N2){
    y2[n2] = normal_rng(f2[n2], sigma);
  }
}

// First run, takes about 1000 seconds for 2000 iterations over 450 training and 100 testing (2 dimensions of input, 2 groups), no warnings by the end of the sampling,
// excpet one warming in the middle of sampling (approximaltly around the 600th iteration), saying there occurs a non-positive definite matrix (and the warning message as well said it doesn't matter if it is not frequently occured)

