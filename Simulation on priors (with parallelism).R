library(tidyverse)
library(rstan)
options(buildtools.check = function(action) TRUE)

# Settings
sample_sizes <- c(10, 20, 40, 80)
sigmas <- c(0.01, 0.05, 0.1, 0.5, 1)
beta <- 0.5
x_test <- seq(-3, 3, length.out = 50)

# Helper function
standardize <- function(x) list(z = (x - mean(x)) / sd(x), mu = mean(x), sd = sd(x))

# Generate all (N, sigma) scenarios
sim_data_list <- list()
for (n in sample_sizes) {
  for (sigma in sigmas) {
    set.seed(1000 + n * 100 + round(sigma * 100))
    x <- rnorm(n)
    y <- beta * x + rnorm(n, sd = sigma)
    
    x_std <- standardize(x)
    y_std <- standardize(y)
    
    x1 <- as.array(x_std$z)
    y1 <- as.vector(y_std$z)
    x2 <- as.array((x_test - x_std$mu) / x_std$sd)
    
    label <- sprintf("N%d_sigma%.2f", n, sigma)
    sim_data_list[[label]] <- list(
      N1 = n,
      x1 = x1,
      y1 = y1,
      N2 = length(x2),
      x2 = x2,
      x_mean = x_std$mu,
      x_sd = x_std$sd,
      y_mean = y_std$mu,
      y_sd = y_std$sd,
      test = tibble(x = x_test)
    )
  }
}

# Stan setup
model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_uniform.stan")
model_labels <- c("horseshoe", "normal", "uniform")

# Fit
fit_results <- list()
for (m in seq_along(model_files)) {
  model_name <- model_labels[m]
  fit_results[[model_name]] <- list()
  
  for (label in names(sim_data_list)) {
    cat("Fitting", model_name, "model on", label, "\n")
    
    stan_input <- sim_data_list[[label]][1:5]  # only pass required parts
    fit <- tryCatch({
      stan(
        file = model_files[m],
        data = stan_input,
        iter = 2000, chains = 1,
        seed = 1000 + m * 10 + as.integer(as.factor(label)),
        refresh = 100
      )
    }, error = function(e) {
      message("  [ERROR] ", e$message)
      NULL
    })
    
    fit_results[[model_name]][[label]] <- fit
  }
}

# ==== 4. Extract posterior alpha values ====
posterior_df <- data.frame()

for (model in names(fit_results)) {
  for (label in names(fit_results[[model]])) {
    post <- rstan::extract(fit_results[[model]][[label]])
    
    # Ensure alpha exists
    if (!"alpha" %in% names(post)) next
    
    df <- tibble(
      alpha = post$alpha,
      model = model,
      data_label = label
    ) %>%
      separate(data_label, into = c("N", "sigma"), sep = "_") %>%
      mutate(
        N = as.integer(str_remove(N, "N")),
        sigma = as.numeric(str_remove(sigma, "sigma"))
      )
    
    posterior_df <- bind_rows(posterior_df, df)
  }
}

# ==== 5. Plotting ====
posterior_df$model <- factor(posterior_df$model, levels = c("horseshoe", "normal", "uniform"))

ggplot(posterior_df, aes(x = model, y = alpha^2, fill = model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5) +
  facet_grid(sigma ~ N, labeller = label_both) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = expression("Posterior Distribution of " * alpha^2),
    x = "Model",
    y = expression(alpha^2)
  )


###################
rm(list = ls())
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores(), buildtools.check = function(action) TRUE)

# Parameters
sample_sizes <- c(10, 20, 40, 80)
sigmas <- c(0.01, 0.05, 0.1, 0.5, 1)
replicates <- 100
beta <- 0.5
x_test <- seq(-3, 3, length.out = 100)
true_alpha2 <- 0  # ground truth

# Stan model files
model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_uniform.stan")
model_labels <- c("horseshoe", "normal", "uniform")

# Helper to standardize
standardize <- function(x) list(z = (x - mean(x)) / sd(x), mu = mean(x), sd = sd(x))

# Output collector
mse_df <- tibble()

for (m in seq_along(model_files)) {
  model_name <- model_labels[m]
  cat("Compiling model:", model_name, "\n")
  model_fit <- stan_model(file = model_files[m])  # compile once
  
  for (n in sample_sizes) {
    for (sigma in sigmas) {
      for (rep in 1:replicates) {
        set.seed(1000 + n * 100 + round(sigma * 100) + rep)
        x <- rnorm(n)
        y <- beta * x + rnorm(n, sd = sigma)
        
        x_std <- standardize(x)
        y_std <- standardize(y)
        
        x1 <- as.array(x_std$z)
        y1 <- as.vector(y_std$z)
        x2 <- as.array((x_test - x_std$mu) / x_std$sd)
        
        stan_data <- list(
          N1 = n,
          x1 = x1,
          y1 = y1,
          N2 = length(x2),
          x2 = x2
        )
        
        cat(sprintf("Fitting %s | N=%d sigma=%.2f rep=%d\n", model_name, n, sigma, rep))
        fit <- tryCatch({
          sampling(model_fit, data = stan_data, iter = 1000, chains = 1,
                   seed = 3000 + m * 10 + rep, refresh = 0)
        }, error = function(e) {
          message("  [ERROR]: ", e$message)
          return(NULL)
        })
        
        if (!is.null(fit)) {
          post <- rstan::extract(fit)
          mse <- mean((post$alpha^2 - true_alpha2)^2)
          
          mse_df <- bind_rows(mse_df, tibble(
            model = model_name,
            N = n,
            sigma = sigma,
            replicate = rep,
            mse = mse
          ))
        }
      }
    }
  }
}


# Plot2
ggplot(mse_df, aes(x = model, y = mse, fill = model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.5) +
  facet_grid(sigma ~ N, labeller = label_both) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = expression("Posterior MSE of " * alpha^2 * " Across Replicates"),
    x = "Model",
    y = "Posterior MSE"
  )

########################################
rm(list = ls())
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores(), buildtools.check = function(action) TRUE)
rstan::rstan_options(auto_write = TRUE)

set.seed(42)
# Parameters
sample_sizes <- c(10, 20, 40, 80)
sigmas <- c(0.01, 0.05, 0.1, 0.5, 1)
replicates <- 200
nonlinearity_levels <- 0:3  # 0 = linear, 3 = most nonlinear
beta <- 0.5
x_test <- seq(-3, 3, length.out = 1)

# Stan model files
model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_uniform.stan")
model_labels <- c("horseshoe", "normal", "uniform")

# Function to simulate y with increasing nonlinearity
generate_y <- function(x, level, beta = 0.5, sigma = 1) {
  w <- level / max(nonlinearity_levels)  # [0, 1]
  base <- beta * x
  wiggle <- w * 2 * sin(pi * x) + x^2 # fewer waves, larger amplitude
  f <- base + wiggle
  f + rnorm(length(x), sd = sigma)
}

# Standardization helper
standardize <- function(x) list(z = (x - mean(x)) / sd(x), mu = mean(x), sd = sd(x))

# Results collector
posterior_summary_list <- list()
row_id <- 1

# Main simulation loop
for (m in seq_along(model_files)) {
  model_name <- model_labels[m]
  cat("Compiling model:", model_name, "\n")
  model_fit <- stan_model(file = model_files[m])
  
  for (nonlin in nonlinearity_levels) {
    for (n in sample_sizes) {
      for (sigma in sigmas) {
        for (rep in 1:replicates) {
          set.seed(10000 + m * 1e5 + nonlin * 1e4 + n * 1e2 + as.integer(sigma * 100) * 10 + rep)
          x <- rnorm(n)
          y <- generate_y(x, level = nonlin, beta = beta, sigma = sigma)
          
          x_std <- standardize(x)
          y_std <- standardize(y)
          
          stan_data <- list(
            N1 = n,
            x1 = as.array(x_std$z),
            y1 = as.vector(y_std$z),
            N2 = length(x_test),
            x2 = as.array((x_test - x_std$mu) / x_std$sd)
          )
          
          cat(sprintf("Fitting %s | N=%d σ=%.2f rep=%d nonlin=%d\n", model_name, n, sigma, rep, nonlin))
          fit <- tryCatch({
            sampling(model_fit, data = stan_data, iter = 5e3, chains = 1, refresh = 0)
          }, error = function(e) {
            message("  [ERROR]: ", e$message)
            NULL
          })
          
          if (!is.null(fit)) {
            post <- rstan::extract(fit, pars = "alpha")
            alpha2_samples <- post$alpha^2
            posterior_summary_list[[row_id]] <- tibble(
              model = model_name,
              N = n,
              sigma = sigma,
              replicate = rep,
              nonlinearity = paste0("level_", nonlin),
              alpha2_mean = mean(alpha2_samples),
              alpha2_median = median(alpha2_samples)
            )
            row_id <- row_id + 1
          }
        }
      }
    }
  }
}
Sys.time()


# Final summary dataframe
posterior_summary_df <- bind_rows(posterior_summary_list)
posterior_summary_df$model <- factor(posterior_summary_df$model, levels = model_labels)
posterior_summary_df$nonlinearity <- factor(posterior_summary_df$nonlinearity,
                                            levels = paste0("level_", nonlinearity_levels))
session_info <- sessionInfo()
save(posterior_summary_df, session_info, file = "posterior_summary.RData")

for (level in levels(posterior_summary_df$nonlinearity)) {
  cat("Plotting for nonlinearity:", level, "\n")
  
  ggplot(
    posterior_summary_df %>% filter(nonlinearity == level),
    aes(x = model, y = alpha2_mean, fill = model)
  ) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.5) +
    facet_grid(sigma ~ N, labeller = label_both) +
    coord_cartesian(ylim = c(0, 5)) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste("Posterior Mean of alpha^2 –", level),
      x = "Model",
      y = expression("Posterior Mean of " * alpha^2)
    ) -> p
  
  print(p)
}

####### parallel running (for multi core machine)
getwd()
rm(list = ls())
parallel::detectCores()

library(tidyverse)
library(rstan)
library(parallel)
library(R.utils)  # for withTimeout
options(mc.cores = 0.9*parallel::detectCores(), buildtools.check = function(action) TRUE)
rstan::rstan_options(auto_write = TRUE)

set.seed(42)

# Parameters
sample_sizes <- c(10, 20, 40, 80)
sigmas <- c(0.01, 0.05, 0.1, 0.5, 1)
replicates <- 200
nonlinearity_levels <- 0:3
beta <- 0.5
x_test <- seq(-3, 3, length.out = 1)

model_files <- c("GPR_horseshoe.stan", "GPR_normal.stan", "GPR_uniform.stan")
model_labels <- c("horseshoe", "normal", "uniform")

generate_y <- function(x, level, beta = 0.5, sigma = 1) {
  w <- level / max(nonlinearity_levels)
  base <- beta * x
  wiggle <- w * 2 * sin(pi * x) + x^2
  f <- base + wiggle
  f + rnorm(length(x), sd = sigma)
}

standardize <- function(x) list(z = (x - mean(x)) / sd(x), mu = mean(x), sd = sd(x))

# Compile models once
compiled_models <- lapply(model_files, stan_model)
names(compiled_models) <- model_labels

# All combinations
task_df <- expand.grid(
  model_index = seq_along(model_labels),
  nonlin = nonlinearity_levels,
  n = sample_sizes,
  sigma = sigmas,
  rep = seq_len(replicates)
)
task_df$id <- seq_len(nrow(task_df))  # add sequential index

run_sim <- function(task_row) {
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
  
  x <- rnorm(n)
  y <- generate_y(x, level = nonlin, beta = beta, sigma = sigma)
  
  x_std <- standardize(x)
  y_std <- standardize(y)
  
  stan_data <- list(
    N1 = n,
    x1 = as.array(x_std$z),
    y1 = as.vector(y_std$z),
    N2 = length(x_test),
    x2 = as.array((x_test - x_std$mu) / x_std$sd)
  )
  
  log_line <- sprintf("[%s] %d/%d | PID=%d | Time=%s | N=%d σ=%.2f rep=%d nonlin=%d\n",
                      model_name, task_id, total_tasks, Sys.getpid(), format(Sys.time(), "%H:%M:%S"),
                      n, sigma, rep, nonlin)
  write(log_line, file = "progress.log", append = TRUE)
  
  fit <- tryCatch({
    withTimeout({
      sampling(model_fit, data = stan_data, iter = 5000, chains = 1, refresh = 0)
    }, timeout = 1200, onTimeout = "error")  # 1 hour timeout
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

# Parallel execution
results <- mclapply(split(task_df, seq_len(nrow(task_df))), run_sim, mc.cores = parallel::detectCores())

# Combine results
posterior_summary_df <- bind_rows(results)

# Save
posterior_summary_df$model <- factor(posterior_summary_df$model, levels = model_labels)
posterior_summary_df$nonlinearity <- factor(posterior_summary_df$nonlinearity, levels = paste0("level_", nonlinearity_levels))
session_info <- sessionInfo()
save(posterior_summary_df, session_info, file = "posterior_summary.RData")
