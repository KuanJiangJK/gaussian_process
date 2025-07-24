## GP Prior comparision
library(tidyverse)
library(rstan)
library(ggh4x) #Useful for boxplot: https://stackoverflow.com/questions/4106614/how-to-control-ylim-for-a-faceted-plot-with-different-scales-in-ggplot2

# Simulate data
generate_gp_data <- function(t, N = 50, sigma = 1) {
  x <- seq(-1, 1, length.out = N)
  linear_part <- function(x) 3 * x
  nonlinear_part <- function(x) 3 * sin(3 * x) + x^2
  w <- (t - 1) / 5
  f <- (1 - w) * linear_part(x) + w * nonlinear_part(x)
  y <- f + rnorm(N, 0, sigma)
  tibble(x = x, y = y, t = paste0("t=", t))
}
df_all <- map_dfr(1:6, generate_gp_data)

prepare_stan_data <- function(df) {
  N <- nrow(df)
  idx_train <- sample(N, size = floor(0.5 * N))
  idx_test <- setdiff(seq_len(N), idx_train)
  df_train <- df[idx_train, ]
  df_test <- df[idx_test, ]
  
  x_mean <- mean(df_train$x)
  x_sd <- sd(df_train$x)
  y_mean <- mean(df_train$y)
  y_sd <- sd(df_train$y)
  
  standardize <- function(x, mu, sd) (x - mu) / sd
  
  df_train <- df_train %>%
    mutate(x_std = standardize(x, x_mean, x_sd),
           y_std = standardize(y, y_mean, y_sd))
  df_test <- df_test %>%
    mutate(x_std = standardize(x, x_mean, x_sd),
           y_std = standardize(y, y_mean, y_sd))
  
  list(
    N1 = nrow(df_train),
    x1 = df_train$x_std,
    y1 = df_train$y_std,
    N2 = nrow(df_test),
    x2 = df_test$x_std,
    x_mean = x_mean, x_sd = x_sd,
    y_mean = y_mean, y_sd = y_sd,
    test = df_test
  )
}

stan_data_list <- map(1:6, ~{
  df_t <- df_all %>% filter(t == paste0("t=", .x))
  prepare_stan_data(df_t)
})

ggplot(df_all, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 1, color = "#1d336c") +
  facet_wrap(~t, ncol = 3) +
  labs(
    title = "Simulated Datasets with Increasing Nonlinearity",
    x = "x", y = "y"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))


model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_uniform.stan")
model_labels <- c("horseshoe", "normal", "uniform")

fit_results <- list()

for (m in seq_along(model_files)) {
  model_name <- model_labels[m]
  fit_results[[model_name]] <- list()
  for (t in 1:6) {
    cat("Fitting", model_name, "model, dataset t =", t, "\n")
    stan_input <- stan_data_list[[t]]
    fit <- stan(
      file = model_files[m],
      data = stan_input[1:5], # only pass N1, x1, y1, N2, x2
      iter = 2000, chains = 1, seed = 1000 + m * 10 + t
    )
    fit_results[[model_name]][[t]] <- fit
  }
}

######################################################
posterior_df <- data.frame()
for (model in names(fit_results)) {
  for (t in 1:6) {
    post <- rstan::extract(fit_results[[model]][[t]])
    posterior_df <- bind_rows(
      posterior_df,
      tibble(value = post$alpha^2, param = "alpha^2", t = paste0("t=", t), model = model),
      tibble(value = post$sigma^2, param = "sigma^2", t = paste0("t=", t), model = model),
      tibble(value = post$rho, param = "rho", t = paste0("t=", t), model = model)
    )
  }
}

posterior_df$t <- factor(posterior_df$t, levels = paste0("t=", 1:6))
posterior_df$model <- factor(posterior_df$model, levels = model_labels)

theme_classic_box <- theme_classic(base_size = 8) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    # axis.line = element_line(colour = "black", linewidth = 0),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "top"
  )

ggplot(posterior_df, aes(x = model, y = value)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5) +
  # facet_grid(param ~ t, scales = "free_y") +
  ggh4x::facet_grid2(param ~ t, scales = "free_y") +
  ggh4x::facetted_pos_scales(
    y = list(
      param == "rho" ~ scale_y_continuous(limits = c(0, 5)),
      param == "alpha^2" ~ scale_y_continuous(limits = c(0, 5))
    )) +
  theme_classic_box +
  labs(title = "Posterior Comparison", y = "Posterior Value", x = "Model") 

gp_draws_df <- data.frame()
for (model in names(fit_results)) {
  for (t in 1:6) {
    post <- rstan::extract(fit_results[[model]][[t]])
    x2 <- stan_data_list[[t]]$test$x_std
    y_mean <- stan_data_list[[t]]$y_mean
    y_sd <- stan_data_list[[t]]$y_sd
    x_orig <- stan_data_list[[t]]$test$x
    f2_mean <- apply(post$f2, 2, mean) * y_sd + y_mean
    y2_low <- apply(post$y2, 2, quantile, probs = 0.025) * y_sd + y_mean
    y2_high <- apply(post$y2, 2, quantile, probs = 0.975) * y_sd + y_mean
    
    gp_draws_df <- bind_rows(gp_draws_df,
                             tibble(x = x_orig, y = f2_mean,
                                    lower = y2_low, upper = y2_high,
                                    model = model, t = paste0("t=", t))
    )
  }
}

ggplot(gp_draws_df, aes(x = x, y = y, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 0.7) +
  facet_grid(model~ t) +
  labs(title = "Posterior Predictive Functions", y = "Predicted y", x = "x") +
  theme_classic_box

