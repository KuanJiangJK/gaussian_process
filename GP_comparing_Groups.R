library(rstan)
rstan_options(auto_write = TRUE)
set.seed(123)
n1 <- 30
n2 <- 30
n3 <- 50
x1 <- sort(runif(n1, 0, 5))
x2 <- sort(runif(n2, 0, 5))
x3 <- sort(runif(n3, 0, 5))

beta1 <- c(1.0, 0.5)  # intercept and slope for group 1
beta2 <- c(0.5, 0.8)  # intercept and slope for group 2
beta3 <- c(0.4, -0.1)  # intercept and slope for group 3
sigma1 <- 0.1
sigma2 <- 0.1
sigma3 <- 0.1
alpha2 <- 0.5
M1 <- 1.0
M2 <- 2.0
M3 <- 3.0

# Construct covariance matrices
make_cov <- function(x, M) {
  distmat <- as.matrix(dist(x))
  exp(-0.5 * M * distmat^2)
}
K1 <- make_cov(x1, M1)
K2 <- make_cov(x2, M2)
K3 <- make_cov(x3, M3)

# when added the measurement error, there will be divergent transition problems
y1 <- as.vector(mvrnorm(1, mu = beta1[1] + beta1[2] * x1, Sigma = alpha2 * K1 + diag(sigma1^2, n1)))
y2 <- as.vector(mvrnorm(1, mu = beta2[1] + beta2[2] * x2, Sigma = alpha2 * K2 + diag(sigma2^2, n2)))
y3 <- as.vector(mvrnorm(1, mu = beta3[1] + beta3[2] * x3, Sigma = alpha2 * K3 + diag(sigma3^2, n3)))

# Combine
group <- c(rep(1, n1), rep(2, n2),rep(3, n3))
x_all <- c(x1, x2, x3)
y_all <- c(y1, y2, y3)

data_list <- list(
  N = n1 + n2 + n3,
  x = x_all,
  y = y_all,
  group = group,
  J = length(unique(group))
)


fit <- stan(
  file  =  "GP_comparing_Groups.stan",
  data = data_list,
  iter = 2000,
  chains = 2,
  seed = 123
)

draw <- rstan::extract(fit)
draw$beta[,1,1]
head(draw$beta)
print(sqrt(1/c(M1,M2,M3)))
plot(fit, pars = "beta", plotfun = "stan_trace", nrow = length(unique(group)))
plot(fit, pars = "rho", plotfun = "stan_trace", nrow = length(unique(group)))
plot(fit, pars = "alpha2", plotfun = "stan_trace", nrow = length(unique(group)))
plot(fit, pars = "sigma", plotfun = "stan_trace", nrow = length(unique(group)))
