library(tidyverse)

theme_classic_box <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "right"
  )

# Parameters
alpha_vals <- c(0, 1)
alpha_j_combinations <- expand.grid(alpha1 = c(0, 1), alpha2 = c(0, 1))
x <- seq(-3, 3, length.out = 300)

# Functions
f_global_linear <- function(x) 0.5 * x
f_global_nonlinear <- function(x) 0.5 * (sin(pi * x) + 0.3 * x^2)
f_group1 <- function(x) 0.5 * (cos(pi * x / 2) + 0.3 * x^2)
f_group2 <- function(x) 0.5 * (sin(pi * x / 2) + 0.1 * x^2)

# Build 8 scenarios
all_scenarios <- bind_rows(lapply(alpha_vals, function(alpha) {
  w_global <- alpha
  global_func <- (1 - w_global) * f_global_linear(x) + w_global * f_global_nonlinear(x)
  
  map2_dfr(alpha_j_combinations$alpha1, alpha_j_combinations$alpha2, function(a1, a2) {
    g_shift <- c(-2, 2)
    group_funcs <- bind_rows(
      tibble(
        x = x,
        y = global_func + a1 * f_group1(x) + g_shift[1],
        group = "G1"
      ),
      tibble(
        x = x,
        y = global_func + a2 * f_group2(x) + g_shift[2],
        group = "G2"
      )
    ) %>%
      mutate(
        alpha = alpha,
        alpha1 = a1,
        alpha2 = a2
      )
    
    global_df <- tibble(
      x = x,
      y = global_func,
      group = "global",
      alpha = alpha,
      alpha1 = a1,
      alpha2 = a2
    )
    
    bind_rows(group_funcs, global_df)
  })
}), .id = NULL) %>%
    mutate(
    group = factor(group, levels = c("G1", "G2", "global")),
    alpha_label = factor(alpha, levels = c(0, 1), labels = c(
      "alpha[global]^2 == 0",
      "alpha[global]^2 != 0"
    )),
    alpha12_label = case_when(
      alpha1 == 0 & alpha2 == 0 ~ "alpha[group1]^2 == 0~','~alpha[group2]^2 == 0",
      alpha1 != 0 & alpha2 == 0 ~ "alpha[group1]^2 != 0~','~alpha[group2]^2 == 0",
      alpha1 == 0 & alpha2 != 0 ~ "alpha[group1]^2 == 0~','~alpha[group2]^2 != 0",
      alpha1 != 0 & alpha2 != 0 ~ "alpha[group1]^2 != 0~','~alpha[group2]^2 != 0"
    ) %>% factor(levels = c(
      "alpha[group1]^2 == 0~','~alpha[group2]^2 == 0",
      "alpha[group1]^2 != 0~','~alpha[group2]^2 == 0",
      "alpha[group1]^2 == 0~','~alpha[group2]^2 != 0",
      "alpha[group1]^2 != 0~','~alpha[group2]^2 != 0"
    ))
  )

# Plot
ggplot(all_scenarios, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(linewidth = 0.4) +
  facet_grid(rows = vars(alpha_label), cols = vars(alpha12_label), labeller = label_parsed) +
  scale_color_manual(
    values = c("G1" = "#c23726", "G2" = "#1d336c", "global" = "black"),
    labels = c("G1" = "Group 1", "G2" = "Group 2", "global" = "Global")
  ) +
  scale_linetype_manual(
    values = c("G1" = "solid", "G2" = "solid", "global" = "twodash"),
    labels = c("G1" = "Group 1", "G2" = "Group 2", "global" = "Global")
  ) +
  labs(x = "x", y = "y", color = "Function", linetype = "Function") +
  theme_classic_box


