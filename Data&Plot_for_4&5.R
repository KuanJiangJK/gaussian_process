rm(list = ls())
library(MASS)  # for mvrnorm
library(rstan)
library(ggplot2)

################# Data generation ##########################
set.seed(123)
N_obs <- 70 # observations
N_test <- 100 # test points
K < - 2
# K <- 1  # number of covariates
J <- 3  # number of groups

# Generate covariates and groups
X1 <- matrix(runif(N_obs * K, -2, 2), ncol = K) # Note: X1 refers to training, and X2 refers to testing
X2 <- matrix(runif(N_test * K, -2, 2), ncol = K)
group1 <- sample(1:J, N_obs, replace = TRUE)  # Group1 is the grouping vector for training. Membership of an input and corresponding output.
group2 <- sample(1:J, N_test, replace = TRUE) # Group2 is the grouping vector for testing. Membership of an input and corresponding output.

# Kernel function for non-linear scenarios
sq_exp_cov <- function(X, alpha = 1, rho = 1) {
  D <- as.matrix(dist(X))^2 # Distance
  K <- alpha^2 * exp(-0.5 * D / rho^2) # Kernel calculation
  K + diag(1e-6, nrow(K))  # numerical stability
}

### Scenario 1: Linear ANCOVA (y = x + mu) ###
true_mu <- c(1.5, -1, 0.5) # mu is intercepts for each groups. 
beta <- 1.2 
Y1 <- beta * X1[,1] + true_mu[group1] + rnorm(N_obs, 0, 0.3)

# Prepare the dataset for Stan
data_linear <- list(
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

### Scenario 2: Non-linear ANCOVA with constant intercept (y = f(x) + mu) ###
K_f <- sq_exp_cov(X1, alpha = 1, rho = 1.2)
f <- as.vector(mvrnorm(1, rep(0, N_obs), K_f)) # generate f
Y1 <- f + true_mu[group1] + rnorm(N_obs, 0, 0.3) # add the grouping intercept. Stan has similar expression like true_mu[group1]  

data_nl_constant <- list(
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

### Scenario 3: Non-linear ANCOVA with varying intercept (mu varies with x) ###
mu_func <- list(                ## define a list of functions, it is a compact way to define group-specific mean functions
  function(x) 1.5 * sin(0.8*x),  # group 1 
  function(x) -1.0 * cos(1.2*x), # group 2
  function(x) 0.5 * x            # group 3
)                               ## mu_func[[1]](1.57/0.8)
mu_values <- sapply(1:N_obs, function(i) mu_func[[group1[i]]](X1[i,1])) # For each observation i, use the function corresponding to its group to compute a mean value based on its x value
# example of sapply
# sapply(1:5, function(x) x^2)
# Output: 1 4 9 16 25

K_f <- sq_exp_cov(X1, alpha = 1, rho = 1.2)
f <- as.vector(mvrnorm(1, rep(0, N_obs), K_f))
Y1 <- f + mu_values + rnorm(N_obs, 0, 0.3)

data_nl_varying <- list(
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

################# MODEL FITTING... ######################

# GP_ANCOVA 4
# Lienar 
fit_linear <- stan(
  file = "GP_ANCOVA_4.stan",
  data = data_linear,
  iter = 2000,
  chains = 1,
  seed = 123
)

# non-linear fixed intercept
fit_nl_constant <- stan(
  file = "GP_ANCOVA_4.stan",
  data = data_nl_constant,
  iter = 2000,
  chains = 1,
  seed = 123
)

# non-linear varying intercept
fit_nl_varying <- stan(
  file = "GP_ANCOVA_4.stan",
  data = data_nl_varying,
  iter = 2000,
  chains = 1,
  seed = 123
)


# GP_ANCOVA 5
# Lienar 
fit_linear <- stan(
  file = "GP_ANCOVA_5.stan",
  data = data_linear,
  iter = 2000,
  chains = 1,
  seed = 123
)

# non-linear fixed intercept
fit_nl_constant <- stan(
  file = "GP_ANCOVA_5.stan",
  data = data_nl_constant,
  iter = 2000,
  chains = 1,
  seed = 123
)

# non-linear varying intercept
fit_nl_varying <- stan(
  file = "GP_ANCOVA_5.stan",
  data = data_nl_varying,
  iter = 2000,
  chains = 1,
  seed = 123
)


#### PLOT PLOT PLOT ####
### 1. Extract Posterior Samples 
post_linear <- rstan::extract(fit_linear)
post_nl_constant <- rstan::extract(fit_nl_constant)
post_nl_varying <- rstan::extract(fit_nl_varying)

### 2. Plots for alpha 
plot_alpha <- function(object, col = "black", title, olim = max(object, na.rm =TRUE)){
  plot(object, main=title, 
       ylab="alpha", xlab="Iteration", pch=20, col = col,
       ylim = c(0, olim))
  plot(density(object), main= title, 
       xlab="alpha", col= col , lwd=2,
       xlim = c(0, olim))
  abline(v=median(object), col="red", lty=2)
}


plot_g_alpha <- function(object, main = NULL, col = "black", olim = max(object, na.rm = TRUE)){
  nc <- ncol(object)
  for (i in 1:nc){
    plot(object[,i], main = paste0(main, " g_alpha ", i),
         xlab="Iteration", ylab= "g_alpha^2",
         ylim = c(0, olim), pch=20, col = col)
  }
}

plot_g_alpha_density <- function(object, main = NULL, col = "black", olim = max(object, na.rm = TRUE)){
  nc <- ncol(object)
  for (i in 1:nc){
    plot(density(object[,i]), main = paste0(main, " g_alpha ", i),
         xlab="g_alpha",
         xlim = c(0, olim), lwd=2, col = col)
    abline(v=median(object[,i]), col="red", lty=2)
  }
}

par(mfrow=c(3,2))
plot_alpha(post_linear$alpha^2, col = "blue", title = "linear", 4)
plot_alpha(post_nl_constant$alpha^2, col = "red", title = "nl_constant", 4)
plot_alpha(post_nl_varying$alpha^2, col = "green", title = "nl_varying", 4)

par(mfrow=c(3,3))
plot_g_alpha(post_linear$g_alpha^2, main = "nl", col = "blue", olim = 4)
plot_g_alpha(post_nl_constant$g_alpha^2, main = "nl", col = "red", olim = 4)
plot_g_alpha(post_nl_varying$g_alpha^2, main = "nl", col = "green", olim = 4)
plot_g_alpha_density(post_linear$g_alpha^2, main = "nl", col = "blue", olim = 4)
plot_g_alpha_density(post_nl_constant$g_alpha^2, main = "nl", col = "red", olim = 4)
plot_g_alpha_density(post_nl_varying$g_alpha^2, main = "nl", col = "green", olim = 4)

### TRACE PLOT

plot(fit_linear, pars = "alpha", plotfun = "stan_trace") + ggtitle ("Linear ANCOVA")
plot(fit_nl_constant, pars = "alpha", plotfun = "stan_trace") + ggtitle ("Non-Linear Constant")
plot(fit_nl_varying, pars = "alpha", plotfun = "stan_trace") + ggtitle ("Non-Linear Varying")
plot(fit_linear, pars = "g_alpha", plotfun = "stan_trace", ncol = 3) + ggtitle ("Linear ANCOVA")
plot(fit_nl_constant, pars = "g_alpha", plotfun = "stan_trace", ncol = 3) + ggtitle ("Non-Linear Constant")
plot(fit_nl_varying, pars = "g_alpha", plotfun = "stan_trace", ncol = 3) + ggtitle ("Non-Linear Varying")

### 4. Posterior Predictive
# Linear ANCOVA
obs_df_linear <- data.frame(
  x = data_linear$X1,
  y = data_linear$Y1,
  group = as.factor(data_linear$group1)
)

f_test_linear <- post_linear$F[, (nrow(data_linear$X1)+1):(nrow(data_linear$X1)+data_linear$N_test)]
draws_linear <- sample(1:nrow(f_test_linear), 10)

mean_df_linear <- data.frame(
  x = data_linear$X2,
  y = colMeans(f_test_linear),
  group = factor(data_linear$group2)
)

plot_df_linear <- data.frame(
  x = rep(data_linear$X2, 10),
  y = as.vector(t(f_test_linear[draws_linear,])),
  group = factor(data_linear$group2),
  sample_f = rep(1:10, each = data_linear$N_test)
)

plot_linear <- ggplot() +
  geom_point(data = obs_df_linear, aes(x = x, y = y, color = group), size = 1.5) +
  geom_line(data = mean_df_linear, aes(x = x, y = y, color = group), linewidth = 1.2) +
  geom_line(data = plot_df_linear, aes(x = x, y = y, group = interaction(sample_f, group)), 
            color = "gray40", linewidth = 0.5, alpha = 0.5) +
  labs(title = "Linear ANCOVA",
       x = "X", y = "y") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

### Non-linear Constant Intercept ANCOVA 
obs_df_nl_constant <- data.frame(
  x = data_nl_constant$X1,
  y = data_nl_constant$Y1,
  group = as.factor(data_nl_constant$group1)
)

f_test_nl_constant <- post_nl_constant$F[, (nrow(data_nl_constant$X1)+1):(nrow(data_nl_constant$X1)+data_nl_constant$N_test)]
draws_nl_constant <- sample(1:nrow(f_test_nl_constant), 10)

mean_df_nl_constant <- data.frame(
  x = data_nl_constant$X2,
  y = colMeans(f_test_nl_constant),
  group = factor(data_nl_constant$group2)
)

plot_df_nl_constant <- data.frame(
  x = rep(data_nl_constant$X2, 10),
  y = as.vector(t(f_test_nl_constant[draws_nl_constant,])),
  group = factor(data_nl_constant$group2),
  sample_f = rep(1:10, each = data_nl_constant$N_test)
)

plot_nl_constant <- ggplot() +
  geom_point(data = obs_df_nl_constant, aes(x = x, y = y, color = group), size = 1.5) +
  geom_line(data = mean_df_nl_constant, aes(x = x, y = y, color = group), linewidth = 1.2) +
  geom_line(data = plot_df_nl_constant, aes(x = x, y = y, group = interaction(sample_f, group)), 
            color = "gray40", linewidth = 0.5, alpha = 0.5) +
  labs(title = "Non-linear Constant Intercept",
       x = "X", y = "y") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

### Non-linear Varying Intercept ANCOVA 
obs_df_nl_varying <- data.frame(
  x = data_nl_varying$X1,
  y = data_nl_varying$Y1,
  group = as.factor(data_nl_varying$group1)
)

f_test_nl_varying <- post_nl_varying$F[, (nrow(data_nl_varying$X1)+1):(nrow(data_nl_varying$X1)+data_nl_varying$N_test)]
draws_nl_varying <- sample(1:nrow(f_test_nl_varying), 10)

mean_df_nl_varying <- data.frame(
  x = data_nl_varying$X2,
  y = colMeans(f_test_nl_varying),
  group = factor(data_nl_varying$group2)
)

plot_df_nl_varying <- data.frame(
  x = rep(data_nl_varying$X2, 10),
  y = as.vector(t(f_test_nl_varying[draws_nl_varying,])),
  group = factor(data_nl_varying$group2),
  sample_f = rep(1:10, each = data_nl_varying$N_test)
)

plot_nl_varying <- ggplot() +
  geom_point(data = obs_df_nl_varying, aes(x = x, y = y, color = group), size = 1.5) +
  geom_line(data = mean_df_nl_varying, aes(x = x, y = y, color = group), linewidth = 1.2) +
  geom_line(data = plot_df_nl_varying, aes(x = x, y = y, group = interaction(sample_f, group)), 
            color = "gray40", linewidth = 0.5, alpha = 0.5) +
  labs(title = "Non-linear Varying Intercept",
       x = "X", y = "y") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

### Plot together!!
library(gridExtra)
grid.arrange(
  plot_linear, plot_nl_constant, plot_nl_varying,
  ncol = 1  
)
