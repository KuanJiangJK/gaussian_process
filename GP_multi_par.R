library(rstan)
rstan_options(auto_write = TRUE)
N_obs <- 30
K <- 4 # number of parameters (including the 1's colum for the intercept)
x11 <- rnorm(N_obs) # the first index 1 means observation data, thus 2 means testing data point
x12 <- rnorm(N_obs)
x13 <- runif(N_obs, 0, 10)
X1 <- matrix(data = c(rep(1,N_obs), x11, x12, x13), nrow = N_obs, ncol = K)
Y1 <- 2+  0.5 * x11 - 0.4 * x12 + 2* sin(x13)  + rnorm(N_obs, 0,1)


N_test <- 100
x21 <- rnorm(N_test)
x22 <- rnorm(N_test)
x23 <- runif(N_test, 0, 10)
X2 <- matrix(data = c(rep(1,N_test), x21, x22, x23), nrow = N_test, ncol = K)

stan.data <- list(N_obs= N_obs,
                  K = K,
                  X1 = X1,
                  Y1 = Y1,
                  N_test = N_test ,
                  X2 = X2)

sample_STAN <- stan(file = 'GP_multi_par.stan',
                     data = stan.data,
                     chains = 1, 
                     warmup = 1e3,
                     iter = 2e3)
# Chain 1: Gradient evaluation took 0.005385 seconds
# Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 53.85 seconds.
# Chain 1: Adjust your expectations accordingly!
# Chain 1:  Elapsed Time: 303.453 seconds (Warm-up)
# Chain 1:                585.11 seconds (Sampling)
# Chain 1:                888.563 seconds (Total)
# Chain 2: Gradient evaluation took 0.004339 seconds
# Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 43.39 seconds.
# Chain 2: Adjust your expectations accordingly!
# Chain 2:  Elapsed Time: 345.672 seconds (Warm-up)
# Chain 2:                603.312 seconds (Sampling)
# Chain 2:                948.984 seconds (Total)
# Warning messages:
# 1: There were 30 divergent transitions after warmup. See
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: Examine the pairs() plot to diagnose sampling problems


plot(sample_STAN, pars = "beta", plotfun = "stan_trace", nrow = 4)
plot(sample_STAN, pars = "rec_L2", plotfun = "stan_trace", nrow = 3)
plot(sample_STAN, pars = "sigman", plotfun = "stan_trace", nrow = 1)
plot(sample_STAN, pars = c("sigmaf", "sigmaf2"), plotfun = "stan_trace", nrow = 2)

posterior_samples <- rstan::extract(sample_STAN)
eta_draws <- posterior_samples$eta
plot(density(rnorm(1e5,0,1)))
lines(density(eta_draws[sample(1e3,1),]), col = "red")
beta_draws <- posterior_samples$beta
rec_L2_draws <- posterior_samples$rec_L2
sigmaf2_draws <- posterior_samples$sigmaf2
sigman_draws <- posterior_samples$sigman

par(mfrow=c(4,1))
plot(density(beta_draws[,1]))
abline(v = mean(beta_draws[,1]),col = "red")
plot(density(beta_draws[,2]))
abline(v = mean(beta_draws[,2]),col = "red")
plot(density(beta_draws[,3]))
abline(v = mean(beta_draws[,3]),col = "red")
plot(density(beta_draws[,4]))
abline(v = mean(beta_draws[,4]),col = "red")

par(mfrow=c(3,1))
plot(density(rec_L2_draws[,1]))
abline(v = mean(rec_L2_draws[,1]),col = "red")
plot(density(rec_L2_draws[,2]))
abline(v = mean(rec_L2_draws[,2]),col = "red")
plot(density(rec_L2_draws[,3]))
abline(v = mean(rec_L2_draws[,3]),col = "red")

par(mfrow=c(2,1))

plot(density(sigmaf2_draws))
abline(v = mean(sigmaf2_draws),col = "red")
plot(density(sigman_draws))
abline(v = mean(sigman_draws),col = "red")

alpha2_draws_2 <- posterior_samples_62$alpha2
plot(density(alpha2_draws_1), xlim = c(0,2), ylim = c(0, 10) )
lines(density(alpha2_draws_2), xlim = c(0,2), col=2)
plot(alpha2_draws_1, ylim = c(0,10))
plot(alpha2_draws_2, ylim = c(0,10))

