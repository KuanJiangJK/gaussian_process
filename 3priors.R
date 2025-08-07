# Compare 3 priors

par(mfrow = c(1,1))
# Generate the distributions
set.seed(123)
sig1 <- abs(rcauchy(1e5, 0, 1))
sig2 <- abs(rnorm(1e5, 0, 1))
sig3 <- invgamma::rinvgamma(1e5, 1, 1)

beta1 <- rnorm(1e5, 0, sig1)
beta2 <- rnorm(1e5, 10, sig2)
beta3 <- rnorm(1e5, 5, sqrt(sig3))

# Plot the first histogram
hist(beta1, breaks = 1e6, freq = FALSE,
     xlim = c(-5, 15), ylim = c(0,2), 
     col  = "red",
     xlab = expression(beta))

# Overlay the others
hist(beta2, breaks = 1e3, freq = FALSE, col = rgb(0, 0, 1, 0.4), xlim = c(-10, 10), ylim = c(0,2), add = TRUE)
hist(beta3, breaks = 1e5, freq = FALSE, col = rgb(1, 1, 0, 0.4) ,xlim = c(-10, 10), ylim = c(0,2), add = TRUE)
# Add legend
legend("topright", legend = c("Half-Cauchy", "Half-Normal", "Inv-Gamma"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4), rgb(1, 1, 0, 0.4)),
       border = NA)
