# Gaussian Process

This project uses squared-Exponential Kernel for demonstration, whose hyperparameters include marginal standard deviation (also called amplitude, magnitude, coded alpha, it shows as the vertical flucation of GP around the mean function), length scale (coded as rho, shown as the horizontal wiggly in the path). The standard deviation of measurement error is coded as sigma.

Parameters could be determined through optimization or Bayesian hierarchical modelling. This project takes the second approach because this way we can encode our own belief in inference (In the paper, this belief is sparsity in nonlinearity).

By modelling the amplitude alpha through a half-Cauchy, we can have an induced horseshoe-type prior in alpha^2, leading to the sparsity in nonlinearity a posteriori. That is to say, when there is no sufficient evidence to justify the nonlinear term, the alpha^2 would be effectively pulled to close to 0, and the whole model degrades to its mean function, i.e. a linear model with a Gaussian noise.

# Code used in the paper (Stored in Thesis code folder)
## Stan code for models of the paper
GP_ANCOVA_constant.stan corresponds to the model of GP ANCOVA with constant group shifts. This is an nonlinear extension of classical linear ANCOVA. The difference is we use GP to model the reponse-covariates relationship, instead of linear model.

GP_ANCOVA_7.stan, GP_ANCOVA_finnish.stan, corresponds to the model of GP ANCOVA with varying group shift. The former file uses an original horse prior setting, and the latter one uses regularized horseshoe prior (also called Finnish horseshoe) with an inverse-gamma(2,1) placed on c^2 regularizing the large signals.
This model enables group varying shift by constructing a group GP, whose elements are determined by group kernel parameter. This paper also uses a zero-sum kernel, motivated by increasing identifiability. But the creator of these codes recently realized it does nothing to the group kernel because we have assumed independence between groups. The identifiability of this model demands more discussion and thoughts. The output of Stan is positive. I think there is no huge non-identifiability issues now. The reason might be, 1. the group kernel setting is correct, and we shouldn't have worried about identifiability issues; 2. Thanks to the horseshoe prior placed on group alpha, which pulls the group GP to 0, leading to disentanglement between group constant shift mu and group GP when GP is flat. 

## Simulation of priors on simple GP regression
This section aims to show effect of different priors placed on GP magnitude.

parallel_sim_compare_prior.R compares the effect of GP under different priors placed on alpha hyperparameter using the following stan codes "GPR_horseshoe.stan", "GPR_invgamma.stan" (inverse-gamma(1,1) on alpha^2), "GPR_invgamma105.stan"(inverse-gamma(1,0.5) on alpha^2), "GPR_invgamma101.stan" (inverse-gamma(1,0.1) on alpha^2).

## artifical dataset 
This section aims to show effect of GPANCOVA models in faced with different data scenarios. It indeed demonstrates the model's adaptability.

ANCOVA_sim.R simulate 8 artifical datasets with different nonlinearity levels in both global trend and group trend. And then used the GP_ANCOVA_constant.stan, GP_ANCOVA_7.stan, GP_ANCOVA_finnish.stan to fit the datasets.

## kidiq dataset
The thesis uses kidiq.csv for a practical implementation. Subsequently using GP_ANCOVA_linear.stan, GP_ANCOVA_constant.stan, GP_ANCOVA_7.stan to fit the data. 

This section aims to show the GPANCOVA, together with a linear ANCOVA works in a real dataset. It indeed could model the potential nonlinearity in both global and group trend, though, in this case, the linear ANCOVA could have give a group mu_j that reflects group difference. While the GPANCOVA with varying shift could add some more subtle info regarding group difference through its length scale, amplitude. 

GP_ANCOVA_linear.stan is a stan implementation of linear ANCOVA. 

The R code implements these analysis are: nieuw_kidiq_linear.R for data manipulation, calling stan and visualization in linear ANCOVA setting. nieuw_kidiq_constantGPANOVA.R for GPANOCOVA with constant shift, and nieuw_kidiq.R for GPANCOVA with varying shift.



### Useful references

Carvalho, Carlos M, Nicholas G Polson, and James G Scott. “Handling Sparsity via the Horseshoe,” n.d.

Mulder, Joris, and Luis Raúl Pericchi. “The Matrix-F Prior for Estimating and Testing Covariance Matrices.” *Bayesian Analysis* 13, no. 4 (December 1, 2018). <https://doi.org/10.1214/17-BA1092>.

Rasmussen, Carl Edward, and Christopher K. I. Williams. *Gaussian Processes for Machine Learning*. 3. print. Adaptive Computation and Machine Learning. Cambridge, Mass.: MIT Press, 2008.

Schulz, Eric, Maarten Speekenbrink, and Andreas Krause. “A Tutorial on Gaussian Process Regression: Modelling, Exploring, and Exploiting Functions.” *Journal of Mathematical Psychology* 85 (August 2018): 1–16. <https://doi.org/10.1016/j.jmp.2018.03.001>.
