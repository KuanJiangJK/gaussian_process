##################### PLOT PLOT PLOT Functions ###############
# Plots for alpha/sigma/length_scale
plot_para <- function(object, col = "black", breaks = 1e2, title, minlim = min(object, na.rm =TRUE), maxlim = max(object, na.rm =TRUE)){
  par(mar=c(2, 2, 1.5, 0.5))
  if(is.na(ncol(object))){
    object <- as.matrix(object)
    nc <- 1
    }else{
    nc <- ncol(object)
  }
  par(mfrow=c(2,nc))
  for (i in 1:nc){
    plot(object[,i], main = paste0(title,"  ", i),
         xlab="Iteration", ylab = title,
         ylim = c(minlim, maxlim), pch=20, col = col)
  }
  for (i in 1:nc){
    hist(object[,i], breaks, main= paste0(title,"  ", i),
         xlab = title, col = col,
         xlim = c(minlim, maxlim),
         freq = FALSE)
    abline(v=median(object[,i]), col="red", lty=2)
  }
}

########################### kidiq Dataset ###################
# rm(list = ls())
# standardize function
std_col_stan <- function(object, cols = seq(1:ncol(object))){
  for (i in 1:length(cols)){
    n <- cols[i]
    mean <- mean(object[,n])
    s <- sd(object[,n])
    object[,n] <- (object[,n]-mean)/s
  }
  return(object)
}

# data
kidiq <- read.csv("kidiq.csv")
kidiq <- std_col_stan(kidiq, cols = c(1,3,4))
N_obs <- nrow(kidiq)
sub_kidiq <- seq(1, N_obs,1)
K <- 2 # mom IQ, mom age
J <- 2 # highshcool or not
X1 <- as.matrix(kidiq[,c("mom_iq","mom_age")][sub_kidiq,])
colnames(X1) <- NULL
Y1 <- kidiq[,"kid_score"][sub_kidiq]
group1 <- kidiq[, "mom_hs"][sub_kidiq]+1 ## this means, group1 is no highschool, group2 is highschool
N_test <- 200
group2 <- sample(1:2, size = N_test, replace = TRUE)
X2_momiq <- seq(from = min(kidiq$mom_iq), to = max(kidiq$mom_iq), length.out = N_test)
X2_momage <- seq(from = min(kidiq$mom_age), to = max(kidiq$mom_age), length.out = N_test)
X2 <- as.matrix(cbind(X2_momiq, X2_momage))
colnames(X2) <- NULL

# Prepare data for rstan. Without c^2 regularization
kidiq_stan <- list(
  N1 = N_obs,
  N2 = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2
)

kidiq_stan_fit_7 <- stan(
  file = "GP_ANCOVA_7.stan",
  data = kidiq_stan,
  iter = 2000,
  chains = 1,
  seed = 123
)

kidiq_stan_fit_8 <- stan(
  file = "GP_ANCOVA_8.stan",
  data = kidiq_stan,
  iter = 2000,
  chains = 1,
  seed = 123
)

# Prepare user-defined regularization, by setting different c^2, you get different cap for the alpha2
# c^2 should be determined relative to the magnitude of tau^2*lambda^2 (according to the paper)
# where by defination, lambda ~ half-Cauchy(0, 1), truncated to larger than 0.
# Here WE don't have tau^2, which is global shrinkage parameter, so assumed to be just 1.
# the quantile values of Cauchy^2 might be helpful for determining a valid value for c^2
sample_cauchy2 <- abs(rcauchy(1e5, 0,1))^2
quantile(sample_cauchy2, probs = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
#         1%          10%          25%          50%          75%          90%          95%          99% 
#   2.488663e-04 2.656248e-02 1.740023e-01 1.005081e+00 5.806097e+00 3.942509e+01 1.559432e+02 3.761669e+03 

### With tunned c2 for regularized horseshoe prior
kidiq_stan_5 <- list(  # Cap top 5% value of the alpha2 (I am not sure if it really did so. Need verify)
  N1 = N_obs,
  N2 = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2,
  c_global2 =  160,
  c_group2 = 160
)

kidiq_stan_25 <- list( # Cap top 25% value of the alpha2
  N1 = N_obs,
  N2 = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2,
  c_global2 =  5,
  c_group2 = 5
)

kidiq_stan_50 <- list( # Cap top 50% value of the alpha2
  N1 = N_obs,
  N2 = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2,
  c_global2 =  1,
  c_group2 = 1
)

kidiq_stan_75 <- list( # Cap top 75% value of the alpha2
  N1 = N_obs,
  N2 = N_test,
  K = K,
  J = J,
  X1 = X1,
  X2 = X2,
  Y1 = Y1,
  group1 = group1,
  group2 = group2,
  c_global2 =  0.17,
  c_group2 = 0.17
)

kidiq_stan_fit_5 <- stan(
  file = "GP_ANCOVA_9.stan",
  data = kidiq_stan_5,
  iter = 2000,
  chains = 1,
  seed = 123
)
kidiq_stan_fit_25 <- stan(
  file = "GP_ANCOVA_9.stan",
  data = kidiq_stan_25,
  iter = 2000,
  chains = 1,
  seed = 123
)
kidiq_stan_fit_50 <- stan(
  file = "GP_ANCOVA_9.stan",
  data = kidiq_stan_50,
  iter = 2000,
  chains = 1,
  seed = 123
)
kidiq_stan_fit_75 <- stan(
  file = "GP_ANCOVA_9.stan",
  data = kidiq_stan_75,
  iter = 2000,
  chains = 1,
  seed = 123
)


## PLOT kidiq #####

# Extract stan fit
post_kid <- rstan::extract(kidiq_stan_fit_7)
post_kid <- rstan::extract(kidiq_stan_fit_8)

post_kid <- rstan::extract(kidiq_stan_fit_5) # cap top 5%
post_kid <- rstan::extract(kidiq_stan_fit_25) # cap top 25%
post_kid <- rstan::extract(kidiq_stan_fit_50) # cap top 50%
post_kid <- rstan::extract(kidiq_stan_fit_75) # cap top 75%

# plot_para(post_kid$lambda_global^2, breaks = 1e4, title = "lambda_global^2", maxlim = 1) # for GP_ANCOVA_8 only, where c is free to estimated
# plot_para(post_kid$lambda_group^2, breaks = 1e4, title = "lambda_group^2", maxlim = 1) # for GP_ANCOVA_8 only, where c is free to estimated
plot_para(post_kid$sigma^2, breaks = 1e2, title = "sigma^2")
plot_para(post_kid$alpha^2, breaks = 3e2, title = "alpha^2")
plot_para(post_kid$rho, breaks = 1e2, title = "rho")
plot_para(post_kid$mu_j, title = "mu_j group")
plot_para(post_kid$g_alpha^2, breaks = 1e2, title = "alpha^2 group", maxlim = 1.51)
plot_para(post_kid$g_rho^2, breaks = 1e2, title = "rho^2 group", maxlim = 20)

plot(kidiq_stan_fit_7, pars = "alpha", plotfun = "stan_trace")
plot(kidiq_stan_fit_7, pars = "g_alpha", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_7, pars = "rho", plotfun = "stan_trace") 
plot(kidiq_stan_fit_7, pars = "g_rho", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_7, pars = "mu_j", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_7, pars = "sigma", plotfun = "stan_trace", nrow = 3) 

plot(kidiq_stan_fit_8, pars = "alpha", plotfun = "stan_trace")
plot(kidiq_stan_fit_8, pars = "g_alpha", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_8, pars = "rho", plotfun = "stan_trace") 
plot(kidiq_stan_fit_8, pars = "g_rho", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_8, pars = "mu_j", plotfun = "stan_trace", nrow = 3)
plot(kidiq_stan_fit_8, pars = "sigma", plotfun = "stan_trace", nrow = 3) 

## training point for plot
obs_df <- data.frame(
  x = kidiq_stan$X1,
  y = kidiq_stan$Y1,
  group = factor(kidiq_stan$group1,
                 levels = c(1,2),
                 labels = c("no highschool", "highschool"))
)

## testing prediction (function values & predictive samples) for plot
f_test <- post_kid$f2 # function value
y_test <- post_kid$y2 # sampleing value

n_draws <- 6 # number of functions drawn (irrespective of groups. They will be splitted ultimately in ggplot by  "interaction", and thus each group 20 drwas )
draws <- sample(1:nrow(f_test), n_draws)

mean_df <- data.frame(       # mean function/samples to be plot
  x = kidiq_stan$X2,
  f = colMeans(f_test),
  y = colMeans(y_test),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("no highschool", "highschool"))
)

plot_df <- data.frame(       ## functions to be plot (only plot function values, because sample values are wiggly and ugly)
  x = rep(kidiq_stan$X2, n_draws),
  f = as.vector(t(f_test[draws,])),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("no highschool", "highschool")),
  sample_f = rep(1:n_draws, each = kidiq_stan$N2)
)

plot_kidiq <- ggplot() +
  geom_point(data = obs_df, aes(x = x.1, y = y, color = group), size = 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_line(data = mean_df, aes(x = x.1, y = f, color = group), linetype = 1, linewidth = 1.5, alpha = 0.9) +
  geom_line(data = mean_df, aes(x = x.1, y = y, linetype = group), color = "blue", linewidth = 0.8, alpha = 0.7) +
  geom_line(data = plot_df, aes(x = x, y = f, group = interaction(sample_f, group), 
                                color = group), linewidth = 0.5, alpha = 0.5) +
  labs(title = "Non-linear ANCOVA",
       x = "X", y = "y",
       linetype = "predictive sampling mean",
       color = "predictive function mean") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")



