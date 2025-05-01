# F = alpha + beta * X + f
# Yi = Fi + error
# f ~ GP(0, Kernel)
# error ~ N(0, sigma_n^2)

set.seed(123)

# install.packages("rstan")
library(rstan)

# Generate observed data
N_obs <- 30
X_obs <- sort(runif(N_obs, -3, 3))  # sorted for better plots
f_true <- function(x) { sin(x) + 0.5 * x }  # true underlying function
sigma_noise <- 0.5
Y_obs <- f_true(X_obs) + rnorm(N_obs, 0, sigma_noise)

# Generate test data
N_test <- 100
X_test <- seq(-3.5, 3.5, length.out = N_test)

# Prepare data list
data_list <- list(
  N_obs = N_obs,
  X1 = X_obs,
  Y1 = Y_obs,
  N_test = N_test,
  X2 = X_test
)

# Compile and sample
rstan_options(auto_write = TRUE) # write the .rds into the directory. no need to compile next time
options(mc.cores = parallel::detectCores())

fit <- stan(
  file = "GP_single_par.stan",
  data = data_list,
  iter = 2000,
  warmup = 1000,
  chains = 1,
  seed = 123
)

posterior <- rstan::extract(fit) # extract function from stan would easily clushed by the one from tidyverse
posterior_beta <- posterior$beta
posterior_alpha <- posterior$alpha
posterior_F <- posterior$F # this is the posterior function value. It is different from posterior samples, wihch would be N(F, measurement error)

# Visual check
plot(X_obs, Y_obs, pch = 19, main = "Generated Data", xlab = "X", ylab = "Y")
curve(f_true, add = TRUE, col = "blue") # plot the true line

GP_meanf <- colMeans(posterior_F)[31: 130] # gp mean function. (only for test point) 
lines(X_test, GP_meanf, col = "red", lwd = 2)

for( i in 1:30){
  lines(X_test, posterior_F[sample(1000, 1),][31: 130], col=adjustcolor("red", 0.3))
} # 30 samples from posterior distribution of F


lm <- lm(Y_obs ~ X_obs)
abline(a = coef(lm)[[1]], b = coef(lm)[[2]], col = "black") # plot the linear regression line


