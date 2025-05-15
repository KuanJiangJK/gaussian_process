## number of test points affect sampling speed. 
# it also affects the plotting function, change the number accodingly.  
# when preparing data for plot, I have ####sample_f = rep(seq(1,10),each = 200)#### this works when N_test = 200. 
## change it accordingly when change N_test
rm(list = ls())
## GP_ANCOVA

library(MASS)     # for mvrnorm
library(rstan)
library(ggplot2)

set.seed(123)

# Parameters
N_obs <- 70
N_test <- 100
K <- 1            # number of covariates
J <- 3             # number of groups

# Covariates: X1 (training), X2 (test)
X1 <- matrix(runif(N_obs * K, -2, 2), ncol = K)
X2 <- matrix(runif(N_test * K, -2, 2), ncol = K)

# Group assignments
group1 <- sample(1:J, N_obs, replace = TRUE)
group2 <- sample(1:J, N_test, replace = TRUE)

# True parameters
true_mu <- c(1.5, -1, 0.5)
alpha <- 1.0
rho <- 1.2
sigma <- 0.3

# Squared exponential kernel function
kernel <- function(x, y, alpha, rho) {
  sq_dist <- as.matrix(dist(rbind(x, y)))^2
  K <- alpha^2 * exp(-0.5 * sq_dist[1:nrow(x), (nrow(x)+1):nrow(sq_dist)] / rho^2)
  return(K)
}

# Compute covariance matrix for X1
sq_exp_cov <- function(X, alpha, rho) {
  D <- as.matrix(dist(X))^2
  K <- alpha^2 * exp(-0.5 * D / rho^2)
  K + diag(1e-6, nrow(K))  # jitter for numerical stability
}

K_f <- sq_exp_cov(X1, alpha, rho)
f <- as.vector(mvrnorm(1, mu = rep(0, N_obs), Sigma = K_f))

# Generate training data
Y1 <- f + true_mu[group1] + rnorm(N_obs, 0, sigma)

# Prepare Stan data list
stan_data <- list(
  N_obs = N_obs,
  N_test = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2
)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

######### Session for GP_ANOCVA.stan (without horseshoe prior) ####################
fit <- stan(
  file = "GP_ANCOVA.stan",
  data = stan_data,
  iter = 2000,
  chains = 1,
  seed = 123
)


library(rstan)
library(ggplot2)

post <- rstan::extract(fit)
posterior_df <- as.data.frame(fit)

# Parameters to trace
par(mfrow=c(2,3))
plot(fit, pars = "alpha", plotfun = "stan_trace")
plot(fit, pars = "rho", plotfun = "stan_trace")
plot(fit, pars = "sigma", plotfun = "stan_trace")
plot(fit, pars = "mu", plotfun = "stan_trace")

# observed data
obs_df <- data.frame(
  x = X1,
  y = Y1,
  group = as.factor(group1)
)
# test inputs
N_test <- nrow(X2)
f_samples <- post$f_group  # dimension: iterations x N_obs + N_test
# Extract only test part
f_test_samples <- f_samples[, (nrow(X1) + 1):(nrow(X1) + N_test)]
# sampling 10 per group
draws <- sample(1e3, 10)

mean_df <- data.frame(x = X2,
                      y = colMeans(f_test_samples),
                      group = factor(group2))


plot.df <- data.frame(x = rep(X2, times= 10),
                      y = as.vector(t(f_test_samples[draws,])),
                      group = factor(group2),
                      sample_f = rep(seq(1,10),each = 100))



# Plot everything
ggplot() +
  geom_point(data = obs_df, aes(x = x, y = y, color = group), size = 1.5) +
  geom_line(data = mean_df, aes(x = x, y = y, color = group), linewidth = 1.2) + 
  geom_line(data = plot.df, aes(x = x, y = y, group = interaction(sample_f, group), color = group), linewidth = 0.5, alpha = 0.5) +
  labs(title = "Observed Data and Posterior Samples by Group",
       x = "X", y = "f(x)") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

######### Session for GP_ANOCVA_2.stan (with horseshoe prior) ####################
fit2 <- stan(
  file = "GP_ANCOVA_2.stan",
  data = stan_data,
  iter = 2000,
  chains = 1,
  seed = 123
)

post2 <- rstan::extract(fit2)
posterior_df2 <- as.data.frame(fit2)

# Parameters to trace
par(mfrow=c(2,3))
plot(fit2, pars = "alpha", plotfun = "stan_trace")
plot(fit2, pars = "rho", plotfun = "stan_trace")
plot(fit2, pars = "sigma", plotfun = "stan_trace")
plot(fit2, pars = "beta", plotfun = "stan_trace")
plot(fit2, pars = "mu", plotfun = "stan_trace")
par(mfrow=c(1,1))
plot(density(post2$alpha))
abline(v = mean(post2$alpha), col = "red")
abline(v = median(post2$alpha), col = "blue")

# observed data
obs_df <- data.frame(
  x = X1,
  y = Y1,
  group = as.factor(group1)
)
# test inputs
N_test <- nrow(X2)
f_samples <- post2$f_group  # dimension: iterations x N_obs + N_test
# Extract only test part
f_test_samples <- f_samples[, (nrow(X1) + 1):(nrow(X1) + N_test)]
# sampling 10 per group
draws <- sample(1e3, 10)

mean_df <- data.frame(x = X2,
                      y = colMeans(f_test_samples),
                      group = factor(group2))


plot.df <- data.frame(x = rep(X2, times= 10),
                      y = as.vector(t(f_test_samples[draws,])),
                      group = factor(group2),
                      sample_f = rep(seq(1,10),each = 100))



# Plot everything
ggplot() +
  geom_point(data = obs_df, aes(x = x, y = y, color = group), size = 1.5) +
  geom_line(data = mean_df, aes(x = x, y = y, color = group), linewidth = 1.2) + 
  geom_line(data = plot.df, aes(x = x, y = y, group = interaction(sample_f, group), color = group), linewidth = 0.5, alpha = 0.5) +
  labs(title = "Observed Data and Posterior Samples by Group",
       x = "X", y = "f(x)") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")


## Compare the two files ----
post <- rstan::extract(fit)
posterior_df <- as.data.frame(fit)
post2 <- rstan::extract(fit2)
posterior_df2 <- as.data.frame(fit2)

par(mfrow=c(2,2))
plot(density(post$alpha))
abline(v = mean(post$alpha), col = "red")
abline(v = median(post$alpha), col = "blue")
plot(density(post2$alpha))
abline(v = mean(post2$alpha), col = "red")
abline(v = median(post2$alpha), col = "blue")

## Rdata saved as > save.image("~/Documents/R/gaussian_process/GP_ANCOVA_output.RData")

