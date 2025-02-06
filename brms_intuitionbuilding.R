data <- data.frame(
  x = seq(1, 20, by = 0.5),
  y = rpois(length(x), lambda = 5 * x + 10) * rbinom(length(x), size = 1, prob = 0.8)
)

formula <- bf(
  y ~ x,
  hu ~ 1 
)

family <- hurdle_poisson()

fit <- brm(
  formula, data = data, family = family,
  chains = 4, iter = 4000, warmup = 2000, seed = 1,
  core = 4, control = list(max_treedepth = 10, adapt_delta = 0.9),
  prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
  silent = 0
)

summary(fit)

conditional_effects(fit, "x", method = "posterior_epred")
#epred combines the mu and hu parts of the model and takes into account the link functions
#By default, this marginalizes across mean values for numeric data and the reference point for categorical data
#The expected y values (with the caveat that it's marginalized over other parameters)

conditional_effects(fit, "x", method = "posterior_linpred")
#linpred only shows the mu part of the model and doesn't take into account the link function
#(i.e. this is on the log scale)
#If you want to report the constraint on how much DNA is in the water, this is what you use
#The y-axis here has caveats - think about what is being marginalized over. 
#If we attributed all zeros to e.g. PCR error then this is what we'd want to report - how much DNA we predict to be in the water

conditional_effects(fit, "x", method = "posterior_predict")
#If you go measure more data tomorrow, this is what you will see. i.e. this is the falsifyable part
#This is if you want to reproduce this work
