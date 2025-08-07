rm(list = ls())
library(tidyverse)
library(rstan)
library(parallel)

Nplot <- 8 # we want 8 data scenarios
generate_ancova_scenario_set <- function(t, N = 300, J = 3, sigma = 0.5) {
  set.seed(1000+t)
  scenarios <- list()
  X <- runif(N, -3, 3)
  group <- sample(1:J, N, replace = TRUE)
  
  # Global effect: transition from linear to nonlinear (t = 1 to 3)
  grand_linear <- function(x) 0.5*x
  grand_nonlinear <- function(x)    0.5 * (sin(pi * x) +  0.3* x^2) 
  w_global <- if (t <= 4) (t-1) / 3 else 1  # 0, 1/3, 2/3, 1
  global_effect <- function(x) (1 - w_global) * grand_linear(x) + w_global * grand_nonlinear(x)
  
  # constant group shifts
  constant_shift <- c(-3, 0, 3)  
  
  # Group-specific nonlinear deviations
  w_group <- if (t <= 4) 0 else (t - 4)  
  nonlinear_strength <- c(-1, 0, 1) * w_group  
  
  # Group-specific nonlinear function
  group_shift <- function(x, g) {
    nonlinear_strength[g] * 0.5* (cos(pi * x / 2) + 0.3*x^2) + constant_shift[g]
  }
  
  # Generate response
  y <- map_dbl(seq_along(X), ~ group_shift(X[.x], group[.x]) + global_effect(X[.x]) + rnorm(1, 0, sigma))
  scenarios[[t]] <- tibble(
    x = X,
    y = y,
    group = group,
    scenario = paste0("t=", t)
  )
  bind_rows(scenarios)
}

ancova_scenarios <- map_dfr(1:Nplot, generate_ancova_scenario_set)
theme_classic_box <- theme_classic(base_size = 8) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    # axis.line = element_line(colour = "black", linewidth = 0),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "right"
  )
# === Standardize x and y for each scenario ===
ancova_scenarios_group_std <- ancova_scenarios %>%
  group_by(scenario) %>%
  mutate(
    x_std = (x - mean(x)) / sd(x),
    y_std = (y - mean(y)) / sd(y)
  ) %>%
  ungroup()

## true values to overlay
# Generate true curves over a grid
true_curve_df <- map_dfr(1:Nplot, function(t) {
  set.seed(1000 + t)
  x_seq <- seq(-3, 3, length.out = 300)
  group_ids <- 1:3
  
  grand_linear <- function(x) 0.5 * x
  grand_nonlinear <- function(x) 0.5 * (sin(pi * x) + 0.3 * x^2)
  w_global <- if (t <= 4) (t - 1) / 3 else 1
  global_effect <- function(x) (1 - w_global) * grand_linear(x) + w_global * grand_nonlinear(x)
  
  constant_shift <- c(-3, 0, 3)
  w_group <- if (t <= 4) 0 else (t - 4)
  nonlinear_strength <- c(-1, 0, 1) * w_group
  group_shift <- function(x, g) {
    nonlinear_strength[g] * 0.5 * (cos(pi * x / 2) + 0.3 * x^2) + constant_shift[g]
  }
  
  # Expand grid: each group across all x
  expand_grid(
    x = x_seq,
    group = group_ids
  ) %>%
    mutate(
      y = map2_dbl(x, group, ~ global_effect(.x) + group_shift(.x, .y)),
      scenario = paste0("t=", t)
    )
})

# Extract means and sds per scenario
means_sds <- ancova_scenarios_group_std %>%
  group_by(scenario) %>%
  summarise(
    x_mean = mean(x),
    x_sd = sd(x),
    y_mean = mean(y),
    y_sd = sd(y),
    .groups = "drop"
  )

# Join and standardize
true_curve_std <- true_curve_df %>%
  left_join(means_sds, by = "scenario") %>%
  mutate(
    x_std = (x - x_mean) / x_sd,
    y_std = (y - y_mean) / y_sd
  )


ggplot(ancova_scenarios, aes(x = x, y = y, color = factor(group))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~scenario, ncol = 3) +
  labs(title = "ANCOVA Scenarios: Linear to Nonlinear Grand & Group Effects",
       x = "x", y = "y") +
  theme_classic_box

ggplot(ancova_scenarios_group_std, aes(x = x_std, y = y_std, color = factor(group))) +
  geom_point(alpha = 0.6) +
  geom_line(
    data = true_curve_std,
    aes(x = x_std, y = y_std, group = interaction(scenario, group), color = factor(group)),
    linewidth = 0.8, linetype = "solid", alpha = 0.8, inherit.aes = FALSE
  ) +
  facet_wrap(~scenario, ncol = 4) +
  scale_color_manual(
    name = "Group",
    values = c("1" = "#c23726", "2" = "#1d336c", "3" = "#e8bf4d"),
    labels = c("Group 1", "Group 2", "Group 3")
  ) +
  labs(title = NULL, x = "std x", y = "std y") +
  theme_classic_box


### Fit to Stan using mclapply ------
# Stan file path
# stan_file <- "GP_ANCOVA_7.stan"
stan_file <- "GP_ANCOVA_finnish.stan"
# stan_file <- "GP_ANCOVA_constant.stan"

# Split data by scenario (returns a named list)
scenario_list <- split(ancova_scenarios_group_std, ancova_scenarios_group_std$scenario)

# Prepare Stan data
prepare_stan_data <- function(df) {
  N <- nrow(df)
  train_idx <- sample(N, size = floor(0.5 * N))
  test_idx <- setdiff(seq_len(N), train_idx)
  df_train <- df[train_idx, ]
  df_test <- df[test_idx, ]
  
  list(
    J = 3,
    K = 1,
    N1 = nrow(df_train),
    X1 = array(df_train$x_std, dim = c(nrow(df_train), 1)),
    group1 = df_train$group,
    Y1 = df_train$y_std,
    N2 = nrow(df_test),
    X2 = array(df_test$x_std, dim = c(nrow(df_test), 1)),
    group2 = df_test$group
  )
}
# Prepare and store all Stan data objects
stan_data_list <- lapply(scenario_list, prepare_stan_data)
names(stan_data_list) <- names(scenario_list)

# Fit model function for one scenario
fit_one_scenario <- function(name) {
  stan_data <- stan_data_list[[name]]
  stan(
    file = stan_file,
    data = stan_data,
    iter = 2000,
    chains = 1,
    seed = 123,
    refresh = 0
  )
}
# Use multicore apply
options(mc.cores = parallel::detectCores())
fit_results <- mclapply(names(scenario_list), fit_one_scenario, mc.cores = 8)
names(fit_results) <- names(scenario_list)
# save(fit_results, file = "fit_results_ANCOVA_finnish313.RData")
# save(fit_results, file = "fit_results_ANCOVA_finnish21.RData")
# save(fit_results, file = "fit_results_ANCOVA_finnish11.RData")
# save(fit_results, file = "fit_results_ANCOVA_hs.RData")
# save(fit_results, file = "fit_results_ANCOVA_constant.RData")
# load("fit_results_ANCOVA_constant.RData")
# load("fit_results_ANCOVA_hs.RData")
# load("fit_results_ANCOVA_finnish.RData")
###################
# Store the alpha ----
alpha <- lapply(fit_results, function(x) rstan::extract(x, "alpha")$alpha)
alpha2_df <- do.call(rbind, lapply(seq_along(alpha), function(i) {
  data.frame(
    scenario = paste0("t=", i),
    # scenario = i,
    alpha2 = alpha[[i]]^2
  )
}))
# Store the group_alpha ----
g_alpha <- lapply(fit_results, function(x) rstan::extract(x, "g_alpha")$g_alpha)
g_alpha2_df <- do.call(rbind, lapply(seq_along(g_alpha), function(i) {
  data.frame(
    scenario = paste0("t=", i),
    # scenario = i,
    g_alpha2 = g_alpha[[i]]^2
  )
}))
colnames(g_alpha2_df)[2:4] <- c("group1", "group2", "group3")
g_alpha2_long <- tidyr::pivot_longer(
  g_alpha2_df,
  cols = starts_with("group"),
  names_to = "group",
  values_to = "g_alpha2"
)
# Store the group_mu ----
mu_j <- lapply(fit_results, function(x) rstan::extract(x, "mu_j")$mu_j)
mu_j_df <- do.call(rbind, lapply(seq_along(mu_j), function(i) {
  data.frame(
    scenario = paste0("t=", i),
    # scenario = i,
    mu_j = mu_j[[i]]
  )
}))
colnames(mu_j_df)[2:4] <- c("group1", "group2", "group3")
# Reshape to long format so each row is one g_alpha2 value
mu_j_long <- tidyr::pivot_longer(
  mu_j_df,
  cols = starts_with("group"),
  names_to = "group",
  values_to = "mu_j"
)

## plot 
ggplot(alpha2_df, aes(x = scenario, y = alpha2)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(title = expression("Posterior of " * alpha^2),
       x = "Scenario", y = expression(alpha^2)) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_classic_box

ggplot(g_alpha2_long, aes(x = group, y = g_alpha2)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  facet_wrap(~scenario, nrow = 2, scales = "free_y") +
  coord_cartesian(ylim = c(-1, 9)) +
  scale_x_discrete(
    name = "Group",
    labels = c("group1" = "Group 1", "group2" = "Group 2", "group3" = "Group 3")
  ) +
  labs(
    title = NULL,
    y = expression(alpha[group]^2)
  ) +
  theme_classic_box

ggplot(mu_j_long, aes(x = group, y = mu_j)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  facet_wrap(~scenario, nrow = 2, scales = "free_y") +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_x_discrete(
    name = "Group",
    labels = c("group1" = "Group 1", "group2" = "Group 2", "group3" = "Group 3")
  ) +
  labs(
    title = NULL,
    y = expression(mu[group])
  ) +
  theme_classic_box

#################################plot function
all_plots_df <- data.frame()
band_df_all <- data.frame()
band_glob_df_all <- data.frame()

for (t in names(fit_results)) {
  stan_data <- stan_data_list[[t]]
  post <- rstan::extract(fit_results[[t]])
  
  f_test <- post$f2
  y_test <- post$y2
  f2glob_test <- post$f2glob
  draws <- sample(1:nrow(f_test), 10)
  
  obs_df <- data.frame(
    x = stan_data$X1[, 1],
    y = stan_data$Y1,
    group = factor(stan_data$group1),
    model = "model",
    t = t,
    type = "obs"
  )
  
  mean_df <- data.frame(
    x = stan_data$X2[, 1],
    y = colMeans(f_test),
    group = factor(stan_data$group2),
    model = "model",
    t = t,
    type = "mean"
  )
  
  sample_df <- data.frame(
    x = rep(stan_data$X2[, 1], 10),
    y = as.vector(t(f_test[draws, ])),
    group = factor(rep(stan_data$group2, times = 10)),
    draw = rep(1:10, each = stan_data$N2),
    model = "model",
    t = t,
    type = "draw"
  )
  
  glob_df <- data.frame(
    x = stan_data$X2[, 1],
    y = colMeans(f2glob_test),
    group = "global",
    model = "model",
    t = t,
    type = "global"
  )
  
  band_df <- data.frame(
    x = stan_data$X2[,1],
    lower = apply(f_test, 2, quantile, probs = 0.025),
    upper = apply(f_test, 2, quantile, probs = 0.975),
    group = factor(stan_data$group2),
    model = "model",
    t = t
  )
  
  band_glob <- data.frame(
    x = stan_data$X2[,1],
    lower = apply(f2glob_test, 2, quantile, probs = 0.025),
    upper = apply(f2glob_test, 2, quantile, probs = 0.975),
    model = "model",
    t = t
  )
  
  band_df_all <- bind_rows(band_df_all, band_df)
  band_glob_df_all <- bind_rows(band_glob_df_all, band_glob)
  
  all_plots_df <- bind_rows(
    all_plots_df,
    obs_df,
    mean_df,
    sample_df,
    glob_df
  )
}
main_plot <- ggplot() +
  geom_ribbon(data = band_df_all,
              aes(x = x, ymin = lower, ymax = upper, fill = group),
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
  # geom_ribbon(data = band_glob_df_all,
  #             aes(x = x, ymin = lower, ymax = upper),
  #             fill = "black", alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = filter(all_plots_df, type == "obs"),
             aes(x = x, y = y, color = group),
             alpha = 0.5, size = 1.3, show.legend = FALSE) +
  geom_line(data = filter(all_plots_df, type == "draw"),
            aes(x = x, y = y, group = interaction(draw, group), color = group),
            alpha = 0.4, linewidth = 0.4) +
  geom_line(data = filter(all_plots_df, type == "mean"),
            aes(x = x, y = y, color = group),
            linewidth = 1.1, alpha = 0.8) +
  geom_line(data = filter(all_plots_df, type == "global"),
            aes(x = x, y = y),
            color = "black", linewidth = 1, linetype = 1) +
  facet_wrap(~t, ncol = 4) +
  labs(x = "Standardized X", y = "Standardized Y") +
  scale_color_manual(values = c("1" = "#c23726", "2" = "#1d336c", "3" = "#e8bf4d")) +
  scale_fill_manual(values = c("1" = "#c23726", "2" = "#1d336c", "3" = "#e8bf4d")) +
  theme_classic_box +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = c(-5, 5))

print(main_plot)

