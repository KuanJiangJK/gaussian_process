rm(list = ls())
parallel::detectCores()
library(ggplot2)
library(tidyverse)
library(rstan)
library(parallel)
library(R.utils)  # for withTimeout
# options(mc.cores = 5, buildtools.check = function(action) TRUE)
rstan::rstan_options(auto_write = TRUE)

set.seed(42)

# Parameters
sample_sizes <- c(10, 20, 40, 80)
sigmas <- c(0.1, 0.3, 0.5, 1)
replicates <- 100
nonlinearity_levels <- 0:3 # 4 nonlinear levels, from total linear to total nonlinear
beta <- 0.5 # this is linear effect
x_test <- seq(-3, 3, length.out = 1) # only one test point, this is to ensure efficiency. Better practice is to cancel this part, and also change the Stan code. But I am lazy to do so

model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_invgamma.stan") # Stan files
model_labels <- c("horseshoe", "normal", "invgamma") # lables 

generate_y <- function(x, level, beta = 0.5, sigma = 1) {
  w <- level / max(nonlinearity_levels)  # normalized nonlinearity level: 0, 1/3, 2/3, 1
  f_linear <- beta * x
  f_nonlinear <- 0.5 * (sin(pi * x) +  0.3* x^2) # consistent scale nonlinear part
  f <- (1 - w) * f_linear + w * f_nonlinear
  y <- f + rnorm(length(x), sd = sigma)
  return(y)
}

standardize <- function(x) list(z = (x - mean(x)) / sd(x), mu = mean(x), sd = sd(x)) # before passing to Stan we need to standardize
# show the sample of the code ----

# Settings
set.seed(123)
x_vals <- seq(-3, 3, length.out = 300)
n_draws <- 1

# Generate replicated, standardized data
draws_df <- expand.grid(
  level = nonlinearity_levels,
  sigma = sigmas,
  draw = seq_len(n_draws)
) |>
  rowwise() |>
  mutate(data = list({
    y <- generate_y(x_vals, level = level, beta = beta, sigma = sigma)
    y_std <- standardize(y)$z
    x_std <- standardize(x_vals)$z
    tibble(x = x_std, y_std = y_std, level = level, sigma = sigma, draw = draw)
  })) |>
  ungroup() |>
  select(-level, -sigma, -draw) |>
  unnest(cols = data)

# True function (standardized) for each level
true_df <- expand.grid(
  x = x_vals,
  level = nonlinearity_levels,
  sigma = sigmas
) |>
  group_by(level, sigma) |>
  mutate(
    f = (1 - level / max(nonlinearity_levels)) * beta * x +
      (level / max(nonlinearity_levels)) * 0.5 * (sin(pi * x) + 0.3 * x^2),
    x_std = standardize(x)$z,
    y_std = standardize(f)$z
  ) |>
  ungroup()

# Plot
ggplot(draws_df, aes(x = x, y = y_std, group = interaction(draw, level, sigma))) +
  # geom_line(alpha = 0.2, color = "#6c757d") +
  geom_point(alpha = 0.4, color = "#6c757d") +
  geom_line(data = true_df, aes(x = x_std, y = y_std), color = "#c23726", linewidth = 1, inherit.aes = FALSE) +
  facet_grid(level ~ sigma, labeller = label_both) +
  labs(
    title = element_blank(),
    x = "Standardized x", y = "Standardized y"
  ) +
  theme_classic_box +
  theme(strip.text = element_text(face = "bold"))

## -----

# Compile models once
compiled_models <- lapply(model_files, stan_model) # compile them all using stan function stan_model
names(compiled_models) <- model_labels

# All combinations
task_df <- expand.grid( # we want to generate a task matrix detailing each tast
  model_index = seq_along(model_labels), # a column indicting the model being used
  nonlin = nonlinearity_levels, # a column indicting the level of linearity
  n = sample_sizes, # size
  sigma = sigmas, # noise
  rep = seq_len(replicates) # the sequential number of the replicate for this scenario
)
task_df$id <- seq_len(nrow(task_df))  # add sequential index 

run_sim <- function(task_row) { # this function is to retrieve the task info to pass to the data generating and Stan modelling
  model_index <- task_row$model_index
  nonlin <- task_row$nonlin
  n <- task_row$n
  sigma <- task_row$sigma
  rep <- task_row$rep
  task_id <- task_row$id
  total_tasks <- nrow(task_df)
  
  model_name <- model_labels[model_index]
  model_fit <- compiled_models[[model_name]]
  set.seed(10000 + model_index * 1e5 + nonlin * 1e4 + n * 1e2 + as.integer(sigma * 100) * 10 + rep)
  
  # x <- rnorm(n)
  x <- seq(-3, 3, length.out = n)
  y <- generate_y(x, level = nonlin, beta = beta, sigma = sigma)
  
  x_std <- standardize(x)
  y_std <- standardize(y)
  
  stan_data <- list(
    N1 = n,
    x1 = as.array(x_std$z),
    y1 = as.vector(y_std$z),
    N2 = length(x_test),
    x2 = as.array((x_test - x_std$mu) / x_std$sd) # standardize the testing
  )
  
  log_line <- sprintf("[%s] %d/%d | PID=%d | Time=%s | N=%d σ=%.2f rep=%d nonlin=%d\n",
                      model_name, task_id, total_tasks, Sys.getpid(), format(Sys.time(), "%H:%M:%S"),
                      n, sigma, rep, nonlin)
  write(log_line, file = "progress.log", append = TRUE)
  
  fit <- tryCatch({
    withTimeout({
      sampling(model_fit, data = stan_data, iter = 3e3, chains = 1, refresh = 0)
    }, timeout = 1200, onTimeout = "error")  # 20 mins timeout
  }, error = function(e) {
    message(sprintf("  [ERROR task %d]: %s", task_id, e$message))
    return(NULL)
  })
  
  if (!is.null(fit)) {
    post <- rstan::extract(fit, pars = "alpha")
    alpha2_samples <- post$alpha^2
    return(tibble(
      model = model_name,
      N = n,
      sigma = sigma,
      replicate = rep,
      nonlinearity = paste0("level_", nonlin),
      alpha2_mean = mean(alpha2_samples),
      alpha2_median = median(alpha2_samples)
    ))
  } else {
    return(NULL)
  }
}

# Parallel execution, check progress: tail -f progress.log
results <- mclapply(split(task_df, seq_len(nrow(task_df))), run_sim, mc.cores = 6) # use all fck cores

# Combine results
posterior_summary_df <- bind_rows(results)
# Save
posterior_summary_df$model <- factor(posterior_summary_df$model, levels = model_labels)
posterior_summary_df$nonlinearity <- factor(posterior_summary_df$nonlinearity, levels = paste0("level_", nonlinearity_levels))
session_info <- sessionInfo()
save(posterior_summary_df, session_info, file = "posterior_summary.RData")

# rm(list = ls())
load("posterior_summary.RData")
# Set factor level order explicitly
posterior_summary_df$model <- factor(posterior_summary_df$model, levels = c("horseshoe", "normal", "invgamma"))

# Custom color palette
model_colors <- c(
  "horseshoe" = "#c23726",
  "normal" = "#1d336c",
  "invgamma" = "#e8bf4d"
)
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
# Generate the plots, 1 per nonlinearity level
plots <- lapply(levels(posterior_summary_df$nonlinearity), function(nl) {
  df_sub <- posterior_summary_df %>% filter(nonlinearity == nl)
  ggplot(df_sub, aes(x = model, y = alpha2_median, fill = model)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
    facet_grid(sigma ~ N, labeller = label_both) +
    scale_fill_manual(values = model_colors) +
    coord_cartesian(ylim = c(0, 5)) +
    # coord_cartesian(ylim = c(0, 3)) +
    # coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = element_blank(),
      # x = paste0("prior used on ", expression(alpha)," or ",expression(alpha^2)),
      x = expression("Prior placed on " * alpha ~ " or " ~ alpha^2),
      y = expression(median(alpha^2))
    ) +
    theme_classic_box +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
})
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
# plot mean
plots <- lapply(levels(posterior_summary_df$nonlinearity), function(nl) {
  df_sub <- posterior_summary_df %>% filter(nonlinearity == nl)
  ggplot(df_sub, aes(x = model, y = alpha2_mean, fill = model)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
    facet_grid(sigma ~ N, labeller = label_both) +
    scale_fill_manual(values = model_colors) +
    coord_cartesian(ylim = c(0, 8)) +
    # coord_cartesian(ylim = c(0, 4)) +
    # coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = element_blank(),
      # x = paste0("prior used on ", expression(alpha)," or ",expression(alpha^2)),
      x = expression("Prior placed on " * alpha ~ " or " ~ alpha^2),
      y = expression(mean(alpha^2))
    ) +
    theme_classic_box +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
})
# Combine all four with patchwork

