########################### kidiq Dataset ###################
library(patchwork)
library(rstan)
library(tidyverse)
library(ggplot2)


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

kidiq_stan_fit_7 <- stan(
  file = "GP_ANCOVA_7.stan",
  data = kidiq_stan,
  iter = 4000,
  chains = 1,
  seed = 123
)

# Extract stan fit
post_kid <- rstan::extract(kidiq_stan_fit_7)
para_names <- names(post_kid)[1:6]
para_post <- data.frame()
for (para in para_names){
  param <- as.data.frame(post_kid[[para]])
  if (ncol(param) == 1){
    para_post <- rbind(para_post,
                       data.frame(value = param[,1], para = para, group = NA))
  }
  else {
    for (i in 1:ncol(param)){
      para_post <- rbind(para_post, 
                         data.frame(value = param[,i], para = para, group = i))
    }
  }
}

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

p1 <- ggplot(para_post[para_post$para %in% c("alpha", "sigma", "rho"),], aes(x = para, y = value)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_wrap(~ para, ncol = 3, scales = "free", 
             labeller = labeller(para = c(
               "alpha" = expression(alpha),
               "rho" = expression(rho),
               "sigma" = expression(sigma),
               "g_alpha" = "Group α",
               "g_rho" = "Group ρ",
               "mu_j" = expression(mu[j])),
               .default = label_parsed)) +
  scale_x_discrete(labels = c("alpha" = expression(alpha), "rho" = expression(rho), "sigma" = expression(sigma))) +
  theme_classic_box +
  labs(title = element_blank(), y = "Posterior Values", x = element_blank())

p2 <- ggplot(para_post[!(para_post$para %in% c("alpha", "sigma", "rho")), ],
       aes(x = factor(group), y = value)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_wrap(~ para, ncol = 3, scales = "free", 
             labeller = labeller(para = c(
               "g_alpha" = "alpha[j]",
               "g_rho"   = "rho[j]",
               "mu_j"    = "mu[j]"
             ),
               .default = label_parsed)) +
  scale_x_discrete(labels = c("1" = "j = 1 (Highschool False)", "2" = "j = 2 (Highschool True)")) +
  theme_classic_box +
  labs(title = element_blank(), y = "Posterior Values" , x = element_blank())
p1 / p2
# p1 | p2
# # Also works, for learning purpose
# ggplot(para_post[!(para_post$para %in% c("alpha", "sigma", "rho")),], aes(x = interaction(para, group), y = value)) +
#   geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
#   scale_x_discrete(labels = c("g_alpha.1" = "no highschool", "g_alpha.2" = "highschool", 
#                               "g_rho.1" = "no highschool", "g_rho.2" = "highschool", 
#                               "mu_j.1" = "no highschool", "mu_j.2" = "highschool")) +
#   facet_wrap( ~ para, scales = "free") +
#   theme_classic_box +
#   labs(title = "Posterior Boxplots", x = "Posterior Value", y = "Parameter")


## training point for plot
obs_df <- data.frame(
  x = kidiq_stan$X1,
  y = kidiq_stan$Y1,
  group = factor(kidiq_stan$group1,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True")))

obs_df[, 1:3] <- mapply(function(x, mean, sd) x*sd + mean,
                        obs_df[, 1:3], # x
                        std_stats[1,], # mean
                        std_stats[2,],
                        SIMPLIFY = TRUE)  # sd
X2_revert <- X2
X2_revert[,1:2] <- mapply(function(x, mean, sd) x*sd + mean,
                          as.data.frame(X2_revert), # x
                  std_stats[1,1:2], # mean
                  std_stats[2,1:2],
                  SIMPLIFY = TRUE)  # sd
## testing prediction (function values & predictive samples) for plot
f_test <- post_kid$f2 # function value
f_test <- f_test * std_stats["sd", "kid_score"] + std_stats["mean", "kid_score"]
# f_test <- sapply(f_test, function(x) std_stats[1,3]+std_stats[2,3]*f_test) # this is not good, very inefficient
y_test <- post_kid$y2 # sampleing value
y_test <- y_test * std_stats["sd", "kid_score"] + std_stats["mean", "kid_score"]
f_glob <- post_kid$f2glob # sampleing value
f_glob <- f_glob * std_stats["sd", "kid_score"] + std_stats["mean", "kid_score"]


n_draws <- 6 # number of functions drawn (irrespective of groups. They will be splitted ultimately in ggplot by  "interaction", and thus each group 20 drwas )
draws <- sample(1:nrow(f_test), n_draws)

mean_df <- data.frame(       # mean function/samples to be plot
  x = X2_revert,
  f = colMeans(f_test),
  y = colMeans(y_test),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True")))

f_test_ci <- apply(f_test, 2, quantile, probs = c(0.025, 0.975))  # 2 x N2 matrix
mean_df$lower <- f_test_ci[1, ]
mean_df$upper <- f_test_ci[2, ]

mean_glob_df <- data.frame(       # mean function/samples to be plot
  x = X2_revert,
  y = colMeans(f_glob),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True"))
)


plot_df <- data.frame(       ## functions to be plot (only plot function values, because sample values are wiggly and ugly)
  x = rep(X2_revert[,1], n_draws),
  f = as.vector(t(f_test[draws,])),
  group = factor(kidiq_stan$group2,
                 levels = c(1,2),
                 labels = c("Highschool False", "Highschool True")),
  sample_f = rep(1:n_draws, each = kidiq_stan$N2)
)

plot_kidiq <- ggplot() +
  # observed points
  geom_point(data = obs_df, aes(x = x.1, y = y, color = group),
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_ribbon(data = mean_df, aes(x = x.1, ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
  geom_line(data = mean_df, aes(x = x.1, y = f, color = group),
            linewidth = 1.5, alpha = 0.9) +
  geom_line(data = plot_df, aes(x = x, y = f, group = interaction(sample_f, group), color = group),
            linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
  geom_line(data = mean_glob_df,
            aes(x = x.1, y = y, color = "Global Effect"),
            linewidth = 0.8, alpha = 0.7, linetype = 1) +
  labs(title = element_blank(),
       x = "Mother's IQ", y = "Child Test Scores",
       color = "Predictive Function Values",
       fill = "95% Confidence Interval") +
  scale_color_manual(
    values = c("Highschool False" = "#c23726",
               "Highschool True" = "#1d336c",
               "Global Effect" = "#fac901")
  ) +
  scale_fill_manual(
    values = c("Highschool False" = "#c23726",
               "Highschool True" = "#1d336c")
  ) +
  theme_classic_box +
  theme(
    legend.position = "inside",
    legend.key.size = unit(0.2, 'cm'),
    legend.title = element_text(size = 7),
    legend.key.width = unit(0.4, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.text = element_text(size = 6),
    legend.position.inside = c(0.85, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.box.background = element_rect(color = "transparent")
  ) + 
  guides(
    color = guide_legend(override.aes = list(
      linewidth = 1,        # consistent line width
      linetype = c("solid", "solid", "solid")
    )),
    fill = guide_legend(override.aes = list(alpha = 0.4))
  ) +
  theme(
    legend.spacing.y = unit(0.5, "cm"),              # more vertical space between items
    legend.key.height = unit(0.4, "cm"),             # taller keys
    legend.key.width = unit(1.2, "cm")               # longer legend lines
  )
print(plot_kidiq)
