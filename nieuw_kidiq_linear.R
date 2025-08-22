########################### kidiq Dataset ###################
library(patchwork)
library(rstan)
library(tidyverse)
library(ggplot2)

rm(list = ls())

kidiq <- read.csv("kidiq.csv") # load the data
vars_to_std <- c("mom_iq", "kid_score", "mom_age") # variables to be standardized.
std_stats <- sapply(kidiq[vars_to_std], function(x) c(mean = mean(x), sd = sd(x))) # lapply: column-wise for DataFrames
std_stats <- std_stats[, c("mom_iq", "mom_age", "kid_score")]
kidiq[vars_to_std] <- scale(kidiq[vars_to_std])

N_obs <- nrow(kidiq)
K <- 2 # mom IQ, mom age
J <- 2 # highshcool or not
X1 <- as.matrix(kidiq[,c("mom_iq","mom_age")])
Y1 <- kidiq[,"kid_score"]
group1 <- kidiq[, "mom_hs"]+1 ## this means, group1 is no highschool, group2 is highschool
N_test <- 100
group2 <- sample(1:J, size = N_test, replace = TRUE)
X2_momiq <- seq(from = min(kidiq$mom_iq), to = max(kidiq$mom_iq), length.out = N_test)
X2_momage <- seq(from = min(kidiq$mom_age), to = max(kidiq$mom_age), length.out = N_test)
X2 <- cbind(X2_momiq, X2_momage)
X1 <- unname(X1)
X2 <- unname(as.matrix(X2))

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

options(mc.cores = parallel::detectCores())
kidiq_stan_fit_7 <- stan(
  file = "GP_ANCOVA_linear.stan",
  data = kidiq_stan,
  iter = 4000,
  chains = 4,
  seed = 123
)

# ---- Extract posterior + plot parameters & predictive functions ----

library(scales)  # for alpha()

# Extract stan fit
post_kid <- rstan::extract(kidiq_stan_fit_7)

# Which parameters to pull if present
expected_names <- c("alpha","rho","sigma","g_alpha","g_rho","mu_j","beta","mu2","y2")
present_names  <- intersect(expected_names, names(post_kid))

# Tidy extraction (no NA warnings)
para_post <- purrr::map_dfr(present_names, function(para) {
  arr <- post_kid[[para]]
  # scalar case: draws x 1 vector
  if (is.null(dim(arr))) {
    tibble::tibble(value = as.numeric(arr), para = para, idx = NA_integer_)
  } else {
    mat <- as.matrix(arr)  # draws x P
    purrr::map_dfr(seq_len(ncol(mat)), function(j)
      tibble::tibble(value = mat[, j], para = para, idx = j)
    )
  }
})

# Make a clean `group`/`name` column for plotting facets/labels
para_post <- para_post |>
  dplyr::mutate(
    group = dplyr::case_when(
      para == "mu_j" ~ idx,
      para == "g_alpha" ~ idx,
      para == "g_rho" ~ idx,
      TRUE ~ NA_integer_
    ),
    name = dplyr::case_when(
      para == "beta" & idx == 1 ~ "beta[mom_iq]",
      para == "beta" & idx == 2 ~ "beta[mom_age]",
      para == "mu_j" & idx == 1 ~ "mu[no_hs]",
      para == "mu_j" & idx == 2 ~ "mu[hs]",
      para %in% c("alpha","rho","sigma") ~ para,
      TRUE ~ paste0(para, "[", idx, "]")
    )
  )


# A compact, consistent theme
theme_classic_box <- theme_classic(base_size = 14) +
  theme(
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    axis.line      = element_blank(),
    strip.background = element_blank(),
    strip.text     = element_text(face = "bold", size = 14),
    plot.title     = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position= "top"
  )

# ---- (A) Key hyperparameters: global + group-specific ----

# Global: alpha, rho, sigma
p_globals <- ggplot(
  dplyr::filter(para_post, para %in% c("alpha","rho","sigma")),
  aes(x = para, y = value)
) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.4) +
  facet_wrap(~ para, ncol = 3, scales = "free",
             labeller = labeller(
               para = c("alpha"="alpha", "rho"="rho", "sigma"="sigma"),
               .default = label_parsed)) +
  scale_x_discrete(labels = c("alpha"=expression(alpha),
                              "rho"  =expression(rho),
                              "sigma"=expression(sigma))) +
  theme_classic_box +
  labs(title = element_blank(), y = "Posterior values", x = NULL)

# Group-specific: g_alpha, g_rho, mu_j
p_groups <- para_post |>
  dplyr::filter(para %in% c("g_alpha","g_rho","mu_j")) |>
  dplyr::mutate(group = factor(group, levels = c(1,2),
                               labels = c("j = 1 (Highschool False)", "j = 2 (Highschool True)"))) |>
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.4) +
  facet_wrap(~ para, ncol = 3, scales = "free",
             labeller = labeller(
               para = c("g_alpha"="alpha[group]",
                        "g_rho"  ="rho[group]",
                        "mu_j"   ="mu[group]"),
               .default = label_parsed)) +
  theme_classic_box +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(title = element_blank(), y = "Posterior values", x = NULL)

param_plot <- p_globals / p_groups
print(param_plot)
ggsave("kidiq_parameters.png", param_plot, width = 8, height = 6, dpi = 300)

# ---- (B) Predictive function values (f, y) with uncertainty bands ----

## training points for plot (back-transform to original scale)
obs_df <- data.frame(
  x = kidiq_stan$X1,
  y = kidiq_stan$Y1,
  group = factor(kidiq_stan$group1,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True"))
)

obs_df[, 1:3] <- mapply(function(x, mean, sd) x*sd + mean,
                        obs_df[, 1:3],
                        std_stats[1,],
                        std_stats[2,],
                        SIMPLIFY = TRUE)

# Copy X2 and back-transform predictors for plotting
X2_revert <- X2
X2_revert[,1:2] <- mapply(function(x, mean, sd) x*sd + mean,
                          as.data.frame(X2_revert),
                          std_stats[1,1:2],
                          std_stats[2,1:2],
                          SIMPLIFY = TRUE)

# Pull predictive objects from Stan and back-transform to original y-scale
# f2: latent function values; y2: posterior predictive; f2glob: global effect
stopifnot(all(c("y2") %in% names(post_kid)))

f_test  <- post_kid$y2

f_test  <- f_test * std_stats["sd","kid_score"] + std_stats["mean","kid_score"]

# Choose a handful of posterior function draws for display
set.seed(1)
n_draws <- 6
draws <- sample(seq_len(nrow(f_test)), n_draws)

# Build plotting frames
mean_df <- data.frame(
  x = X2_revert,
  f = colMeans(f_test),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True"))
)
f_test_ci <- apply(f_test, 2, quantile, probs = c(0.025, 0.975))
mean_df$lower <- f_test_ci[1, ]
mean_df$upper <- f_test_ci[2, ]


plot_df <- data.frame(
  x = rep(X2_revert[,1], n_draws),
  f = as.vector(t(f_test[draws, ])),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True")),
  sample_f = rep(seq_len(n_draws), each = kidiq_stan$N2)
)

# Master plot (x-axis = Mother's IQ; change to x.2 for Mom's Age)
plot_kidiq <- ggplot() +
  # observed points
  geom_point(data = obs_df, aes(x = x.1, y = y, color = group),
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  # pointwise 95% interval for f(x)
  geom_ribbon(data = mean_df, aes(x = x.1, ymin = lower, ymax = upper, fill = group),
              alpha = 0.20) +
  # posterior mean function
  geom_line(data = mean_df, aes(x = x.1, y = f, color = group),
            linewidth = 1.25, alpha = 0.9) +
  # a few posterior function draws for wiggle
  geom_line(data = plot_df, aes(x = x, y = f, group = interaction(sample_f, group), color = group),
            linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
  # global effect (averaged over groups)
  labs(title = element_blank(),
       x = "Mother's IQ", y = "Child Test Scores",
       color = "Predictive Function",
       fill  = "95% Interval (f)") +
  scale_color_manual(
    values = c("Highschool False" = "#c23726",
               "Highschool True"  = "#1d336c",
               "Global Effect"    = "#fac901")
  ) +
  scale_fill_manual(
    values = c("Highschool False" = "#c23726",
               "Highschool True"  = "#1d336c")
  ) +
  theme_classic_box +
  theme(
    legend.position = "inside",
    legend.key.size = unit(1, 'cm'),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.4, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.text = element_text(size = 10),
    legend.position.inside = c(0.8, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.box.background = element_rect(color = "transparent")
  ) + 
  guides(
    color = guide_legend(override.aes = list(
      linewidth = 1,        # consistent line width
      linetype = c("solid",  "solid")
    )),
    fill = guide_legend(override.aes = list(alpha = 0.4))
  ) +
  theme(
    legend.spacing.y = unit(0.5, "cm"),              # more vertical space between items
    legend.key.height = unit(0.4, "cm"),             # taller keys
    legend.key.width = unit(1.2, "cm")               # longer legend lines
  )

print(plot_kidiq)
ggsave("kidiq_predictive_functions.png", plot_kidiq, width = 7.5, height = 4.8, dpi = 300)

