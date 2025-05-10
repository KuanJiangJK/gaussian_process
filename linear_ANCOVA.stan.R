set.seed(123)

# Simulate data
N_per_group <- 30
J <- 3
N <- J * N_per_group

group <- rep(1:J, each = N_per_group)
x <- rnorm(N, mean = 0, sd = 1)

# True parameters
beta <- 2.0       # slope
mu <- c(1.0, 2.0, 3.0)  # group intercepts
sigma <- 1.0

# Generate outcome
y <- mu[group] + beta * x + rnorm(N, 0, sigma)

# Create design matrix with intercept
X1 <- matrix(x[train_idx], ncol = 1)
X2 <- matrix(x[test_idx], ncol = 1)

K <- ncol(X)

# Split into training and test
train_idx <- sample(1:N, size = 70)
test_idx <- setdiff(1:N, train_idx)

# Prepare data for Stan
data_list <- list(
  N_obs = length(train_idx),
  K = K,
  J = J,
  X1 = X1,
  Y1 = y[train_idx],
  group1 = group[train_idx],
  N_test = length(test_idx),
  X2 = X2,
  group2 = group[test_idx]
)

library(rstan)

library(ggplot2)

# Compile model

fit <- stan(file = "linear_ANCOVA.stan",
  data = data_list,
  seed = 123,
  chains = 4,
  iter  = 1000,
  warmup = 500
)


# Extract predictions
draws <- extract(fit)
plot(fit, pars = "beta", plotfun = "stan_trace")
plot(fit, pars = "mu", plotfun = "stan_trace")
