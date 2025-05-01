# Gaussian Process

This project uses squared-Exponential Kernel for demonstration, where parameters include magnitude (coded as sigmaf or alpha), length scale (l, or 1/l^2^), and measurement error (coded as sigman or sigma).

Parameters could be determined through optimization or Bayesian hierarchical modelling. This project takes the second path.

By modelling the magnitude (sigmaf or alpha) through student-t with 1 degree of freedom, sigmaf^2^ or alpha^2^ would be following an F(1, 1) distribution. Thus forming a horseshoe type of prior useful for modelling sparsity of mean function of normal distribution.

### Some useful reference

Carvalho, Carlos M, Nicholas G Polson, and James G Scott. “Handling Sparsity via the Horseshoe,” n.d.

Mulder, Joris, and Luis Raúl Pericchi. “The Matrix-F Prior for Estimating and Testing Covariance Matrices.” *Bayesian Analysis* 13, no. 4 (December 1, 2018). <https://doi.org/10.1214/17-BA1092>.

Rasmussen, Carl Edward, and Christopher K. I. Williams. *Gaussian Processes for Machine Learning*. 3. print. Adaptive Computation and Machine Learning. Cambridge, Mass.: MIT Press, 2008.

Schulz, Eric, Maarten Speekenbrink, and Andreas Krause. “A Tutorial on Gaussian Process Regression: Modelling, Exploring, and Exploiting Functions.” *Journal of Mathematical Psychology* 85 (August 2018): 1–16. <https://doi.org/10.1016/j.jmp.2018.03.001>.
