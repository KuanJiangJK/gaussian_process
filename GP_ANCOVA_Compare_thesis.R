
rm(list = ls())
load(N150itr3e3Compare.RData)
library(tidyverse)
library(rstan)

Nplot <- 6 
generate_ancova_scenario_set <- function(t, N = 150, J = 3, sigma = 1) {
  scenarios <- list()
  X <- runif(N, -1, 1)
  group <- sample(1:J, N, replace = TRUE)
  
  # Global effect: transition from linear to nonlinear (t = 1 to 3)
  grand_linear <- function(x) 5*x
  grand_nonlinear <- function(x)    5*sin(2 * x) + 2*x^2
  w_global <- if (t <= 3) (t - 1) / 2 else 1  # 0, 0.5, 1
  global_effect <- function(x) (1 - w_global) * grand_linear(x) + w_global * grand_nonlinear(x)
  
  # constant group shifts
  constant_shift <- c(-3, 0, 3)  
  
  # Group-specific nonlinear deviations
  w_group <- if (t <= 3) 0 else (t - 3)  
  nonlinear_strength <- c(-1, 0, 1) * w_group  
  
  # Group-specific nonlinear function
  group_shift <- function(x, g) {
    nonlinear_strength[g] * (5*cos(pi * x / 2) + 2*x^2) + constant_shift[g]
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
generate_scenario_type1 <- function(t, N = 150, J = 3, sigma = 1) {
  X <- runif(N, -1, 1)
  group <- sample(1:J, N, replace = TRUE)
  
  w_global <- (t - 1) / 5
  w_group <- 0
  
  grand_linear <- function(x) 5 * x
  grand_nonlinear <- function(x) 5 * sin(2 * x) + 2 * x^2
  constant_shift <- c(-3, 0, 3)
  nonlinear_strength <- c(-1, 0, 1)
  
  global_effect <- function(x) (1 - w_global) * grand_linear(x) + w_global * grand_nonlinear(x)
  group_shift <- function(x, g) constant_shift[g]
  
  y <- sapply(seq_along(X), function(i) group_shift(X[i], group[i]) + global_effect(X[i]) + rnorm(1, 0, sigma))
  
  tibble(x = X, y = y, group = group, scenario = paste0("t=", t), type = "type1")
}
generate_scenario_type2 <- function(t, N = 150, J = 3, sigma = 1) {
  X <- runif(N, -1, 1)
  group <- sample(1:J, N, replace = TRUE)
  
  w_global <- 0
  w_group <- (t - 1) / 5
  
  grand_linear <- function(x) 5 * x
  grand_nonlinear <- function(x) 5 * sin(2 * x) + 2 * x^2
  constant_shift <- c(-3, 0, 3)
  nonlinear_strength <- c(-3, 0, 3)
  
  global_effect <- function(x) grand_linear(x)
  group_shift <- function(x, g) nonlinear_strength[g] * (5 * cos(pi * x / 2) + 2 * x^2) * w_group + constant_shift[g]
  
  y <- sapply(seq_along(X), function(i) group_shift(X[i], group[i]) + global_effect(X[i]) + rnorm(1, 0, sigma))
  
  tibble(x = X, y = y, group = group, scenario = paste0("t=", t), type = "type2")
}
generate_scenario_type3 <- function(t, N = 150, J = 3, sigma = 1) {
  X <- runif(N, -1, 1)
  group <- sample(1:J, N, replace = TRUE)
  
  w_global <- 1
  w_group <- (t - 1) / 5
  
  grand_linear <- function(x) 5 * x
  grand_nonlinear <- function(x) 5 * sin(2 * x) + 2 * x^2
  constant_shift <- c(-3, 0, 3)
  nonlinear_strength <- c(-3, 0, 3)
  
  global_effect <- function(x) grand_nonlinear(x)
  group_shift <- function(x, g) nonlinear_strength[g] * (5 * cos(pi * x / 2) +  x^2) * w_group + constant_shift[g]
  
  y <- sapply(seq_along(X), function(i) group_shift(X[i], group[i]) + global_effect(X[i]) + rnorm(1, 0, sigma))
  
  tibble(x = X, y = y, group = group, scenario = paste0("t=", t), type = "type3")
}

ancova_scenarios <- map_dfr(1:Nplot, generate_ancova_scenario_set)
ggplot(ancova_scenarios, aes(x = x, y = y, color = factor(group))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~scenario, ncol = 3) +
  labs(title = "ANCOVA Scenarios: Linear to Nonlinear Grand & Group Effects",
       x = "x", y = "y") +
  theme_minimal()

df_all <- bind_rows(
  map_dfr(1:Nplot, generate_scenario_type1),
  map_dfr(1:Nplot, generate_scenario_type2),
  map_dfr(1:Nplot, generate_scenario_type3)
)
ggplot(df_all, aes(x = x, y = y, color = factor(group))) +
  geom_point(alpha = 0.6) +
  facet_grid(type ~ scenario) +
  labs(title = "Three ANCOVA Simulation Types over 6 Transitions",
       x = "x", y = "y") +
  theme_minimal()

# ----------------------- 2. Create stan_data list for t = 1:4 ----------------------- #
prepare_stan_data <- function(df) {
  N <- nrow(df)
  train_idx <- sample(N, size = floor(0.5 * N))
  test_idx <- setdiff(seq_len(N), train_idx)
  df_train <- df[train_idx, ]
  df_test <- df[test_idx, ]
  
  # Standardize based on training
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

# Generate and prepare data
stan_data_list <- list()
for (t in 1:Nplot) {
  set.seed(100 + t)
  df <- generate_scenario_type2(t)
  stan_data_list[[t]] <- prepare_stan_data(df)
  assign(paste0("stan_data_", t), stan_data_list[[t]])
}

# ----------------------- 3. Run Stan models for each t ----------------------- #
model_files <- c("GP_ANCOVA_uniform.stan", "GP_ANCOVA_7.stan", "GP_ANCOVA_normal.stan")
model_labels <- c("uniform", "horseshoe", "normal")

fit_results <- list()

for (m in seq_along(model_files)) {
  model_name <- model_labels[m]
  fit_results[[model_name]] <- list()
  
  for (t in 1:Nplot) {
    cat("Fitting model:", model_name, ", scenario t =", t, "\n")
    
    fit <- stan(
      file = model_files[m],
      data = stan_data_list[[t]],
      iter = 2e3,
      chains = 1,
      seed = 1000 + 10 * m + t
    )
    
    fit_results[[model_name]][[t]] <- fit
  }
}

# ----------------------- 4. Plotting diagnostic and functional results ----------------------- #
library(ggplot2)
# Parameters to extract
param_names <- c("alpha", "sigma", "rho", "g_alpha", "mu_j")
# posterior data frame
posterior_data <- data.frame()
for (model in names(fit_results)) {
  for (t in 1:Nplot) {
    post <- rstan::extract(fit_results[[model]][[t]])
    # alpha^2
    posterior_data <- rbind(posterior_data, data.frame(
      value = post$alpha^2,
      param = "alpha^2",
      t = paste0("t=", t),
      model = model
    ))
    # sigma^2
    posterior_data <- rbind(posterior_data, data.frame(
      value = post$sigma^2,
      param = "sigma^2",
      t = paste0("t=", t),
      model = model
    ))
    # rho
    posterior_data <- rbind(posterior_data, data.frame(
      value = post$rho,
      param = "rho",
      t = paste0("t=", t),
      model = model
    ))
    # g_alpha^2[j]
    for (j in 1:ncol(post$g_alpha)) {
      posterior_data <- rbind(posterior_data, data.frame(
        value = post$g_alpha[, j]^2,
        param = paste0("g_alpha^2[", j, "]"),
        t = paste0("t=", t),
        model = model
      ))
    }
    # mu_j[j]
    for (j in 1:ncol(post$mu_j)) {
      posterior_data <- rbind(posterior_data, data.frame(
        value = post$mu_j[, j],
        param = paste0("mu_j[", j, "]"),
        t = paste0("t=", t),
        model = model
      ))
    }
  }
}

# Plot
library(ggplot2)
# Convert t and model to factor
posterior_data$t <- factor(posterior_data$t, levels = paste0("t=", 1:Nplot))
posterior_data$model <- factor(posterior_data$model, levels = c("uniform", "horseshoe", "normal"))
params <- unique(posterior_data$param) # unique para names
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

# for (p in params) {
#   df_param <- posterior_data[posterior_data$param == p, ]
#   
#   range_x <- quantile(df_param$value, 0.97) - quantile(df_param$value, 0)
#   desired_bins <- 100
#   binwidth_val <- range_x / desired_bins
#   bw_val <- range_x / 15
#   
#   p_plot <- ggplot(df_param, aes(x = value, color = model)) +
#     geom_histogram(aes(y = ..density.., fill = model, col = model), binwidth = binwidth_val, alpha = 0.4, position = "identity", linewidth = 0.3) +
#     geom_density(linewidth = 0.7, alpha = 1, bw = bw_val) +
#     geom_vline(data = aggregate(value ~ model + t, data = df_param, FUN = median),
#                aes(xintercept = value, color = model),
#                linetype = "dashed", linewidth = 0.5) +
#     facet_wrap(~t, ncol = 2, scales = "free_y") +
#     coord_cartesian(xlim = c(0, quantile(df_param$value, 0.97))) +
#     scale_color_manual(values = c("uniform" = "#c23726", "horseshoe" = "#1d336c", "normal" = "#0e0e08")) + # #e8bf4d ##0e0e08 ## #6E6E6E
#     scale_fill_manual(values  = c("uniform" = "#c23726", "horseshoe" = "#1d336c", "normal" = "#0e0e08")) +
#     theme_classic_box +
#     theme(legend.position = "top",
#           strip.text = element_text(face = "bold")) +
#     labs(title = paste("Posterior Comparison for", p),
#          x = "Posterior sample", y = "Density")
#   
#   print(p_plot)
# }

## box plot
quantileshown <- c(0.97, 0.98, 0.98, 0.97, 0.97, 0.97,1,1,1)
names(quantileshown) <- params
for (p in params) {
  df_param <- posterior_data[posterior_data$param == p, ]
  # p_plot <- ggplot(df_param, aes(x = model, y = value, color = model)) +
  p_plot <- ggplot(df_param, aes(x = model, y = value)) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.6, alpha = 0.5) +
    facet_wrap(~t, ncol = 2, scales = "free_y") +
    coord_cartesian(ylim = c(quantile(df_param$value, 0), quantile(df_param$value, quantileshown[[p]]))) +
    # scale_color_manual(values = c("uniform" = "#c23726", "horseshoe" = "#1d336c", "normal" = "#0e0e08")) +
    theme_classic_box +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold")) +
    labs(title = paste("Posterior Comparison (Boxplot) for", p),
         x = "Model", 
         y = "Posterior Sample")
  print(p_plot)
}

#### plot functions #### 
library(tidyverse)
all_plots_df <- data.frame()  # will hold everything for ggplot
for (t in 1:Nplot) {
  stan_data <- stan_data_list[[t]]
  obs_df <- data.frame(
    x = stan_data$X1[, 1],
    y = stan_data$Y1,
    group = factor(stan_data$group1)
  )
  
  for (model in c("uniform", "horseshoe", "normal")) {
    post <- rstan::extract(fit_results[[model]][[t]])
    
    f_test <- post$f2 # predicted function value
    y_test <- post$y2 # predicted sampling value
    f2glob_test <- post$f2glob # predicted global component
    draws <- sample(1:nrow(f_test), 10)
    
    mean_df <- data.frame(
      x = stan_data$X2[, 1],
      f = colMeans(f_test),
      y = colMeans(y_test),
      group = factor(stan_data$group2)
    )
    
    # Global-only function mean (black line)
    glob_df <- data.frame(
      x = stan_data$X2[, 1],
      y = colMeans(f2glob_test),
      group = "global",  # use special group label, instead of 1, 2, 3
      model = model,
      t = paste0("t=", t),
      type = "global"
    )
    
    sample_df <- data.frame(
      x = rep(stan_data$X2[, 1], 10),
      f = as.vector(t(f_test[draws, ])),
      group = factor(rep(stan_data$group2, times = 10)),
      draw = rep(1:10, each = stan_data$N2)
    )
    # uncertainty band
    band_df <- data.frame(
      x = stan_data$X2[,1],
      lower = apply(f_test, 2, quantile, probs = 0.025),
      upper = apply(f_test, 2, quantile, probs = 0.975),
      group = factor(stan_data$group2),
      model = model,
      t = paste0("t=", t)
    )
    # Append to cumulative uncertainty band data
    if (!exists("band_df_all")) {
      band_df_all <- band_df ## for the first iteration, create the band_df_all
    } else { ## for the following, append
      band_df_all <- bind_rows(band_df_all, band_df)
    }
    # Global-only ribbon
    band_glob <- data.frame(
      x = stan_data$X2[,1],
      lower = apply(f2glob_test, 2, quantile, probs = 0.025),
      upper = apply(f2glob_test, 2, quantile, probs = 0.975),
      model = model,
      t = paste0("t=", t)
    )
    # Append to cumulative global ribbon data
    if (!exists("band_glob_df_all")) {
      band_glob_df_all <- band_glob
    } else {
      band_glob_df_all <- bind_rows(band_glob_df_all, band_glob)
    }
    
    obs_df$model <- model # add another column documenting the present model name
    obs_df$t <- paste0("t=", t) # add another column documenting the present dataset
    
    mean_df$model <- model # same 
    mean_df$t <- paste0("t=", t) # same
    
    sample_df$model <- model # same
    sample_df$t <- paste0("t=", t) # same
    
    all_plots_df <- bind_rows( # append all the data, with label by mutate
      all_plots_df,
      obs_df %>% mutate(type = "obs"),
      mean_df %>% select(-y) %>% rename(y = f) %>% mutate(type = "mean"),
      sample_df %>% rename(y = f) %>% mutate(type = "draw"),
      glob_df
    )
  }
}

# Plot functions 
main_plot <- ggplot() +
  geom_ribbon(data = band_df_all,
              aes(x = x, ymin = lower, ymax = upper, fill = group),
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
  geom_ribbon(data = band_glob_df_all,
              aes(x = x, ymin = lower, ymax = upper),
              fill = "black", alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE) +
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
  facet_grid(t ~ model) + ## a plot matrix with two discrete variable t and model. these two should be stored in the dataset
  # labs(title = "Posterior Function Samples and Mean by Scenario and Prior",
  #     x = "Standardized X", y = "Standardized Y") +
  labs(x = "Standardized X", y = "Standardized Y") +
  scale_color_manual(values = c("1" = "#c23726", "2" = "#1d336c", "3" = "#e8bf4d")) +
  scale_fill_manual(values = c("1" = "#c23726", "2" = "#1d336c", "3" = "#e8bf4d")) + #"#0e0e08"
  # scale_color_brewer(palette = "Dark2") +
  # scale_fill_brewer(palette = "Dark2") +  # add this line
  theme_classic_box +
  theme(legend.position = "right", legend.key.size = unit(0.8, "lines"), legend.title = element_text(size = 7), legend.text = element_text(size = 5)) +
  ylim(-5,5)
print(main_plot)

# traceplot
library(ggplot2)
trace_params <- c("alpha", "sigma", "rho", "g_alpha", "mu_j")
for (model in names(fit_results)) {
  for (t in 1:Nplot) {
    cat("\nTraceplot for model:", model, "t =", t, "\n")
    fit_trace <- fit_results[[model]][[t]]
    tp <- traceplot(fit_trace, pars = trace_params, inc_warmup = TRUE, col = "#0e0e08") +
      ggtitle(paste("Traceplot for model:", model, ", scenario t =", t)) +
      theme(legend.position = "none") +
      labs(x = "iteration", y = "value") +
      theme_classic_box
    print(tp)
  }
}

