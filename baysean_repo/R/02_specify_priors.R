# =============================================================================
# 02_specify_priors.R
# Bayesian Moderated Mediation Tutorial
# Step 3: Specify Priors
# =============================================================================

library(brms)
library(tidyverse)

# =============================================================================
# Prior specification for multivariate model
# =============================================================================

# Priors for multivariate model with different response families
priors <- c(
  # Mediator equation (Gaussian): TTM ~ BNR * FP
  prior(normal(0, 1.5), class = "b", resp = "TTMz"),
  prior(exponential(1), class = "sigma", resp = "TTMz"),


  # Outcome equation (can be Gaussian or Bernoulli)
  prior(normal(0, 1.5), class = "b", resp = "FS")
  # Note: No sigma prior for Bernoulli family
)

# --- View prior specifications ---
cat("Prior specifications:\n")
print(priors)

# =============================================================================
# Visualize what these priors imply
# =============================================================================

# For standardized variables, Normal(0, 1) means:
# - 68% of prior mass: effect between -1 and +1 SD
# - 95% of prior mass: effect between -2 and +2 SD
# This is "weakly informative" - allows large effects but skeptical of extreme

library(ggplot2)

# Visualize Normal(0, 1) prior
tibble(x = seq(-4, 4, 0.01)) %>%
  mutate(density = dnorm(x, 0, 1)) %>%
  ggplot(aes(x, density)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", alpha = 0.5) +
  annotate("text", x = 0, y = 0.2, label = "68%", size = 4) +
  annotate("text", x = 0, y = 0.05, label = "95%", size = 4) +
  theme_minimal() +
  labs(title = "Normal(0, 1) Prior for Standardized Coefficients",
       subtitle = "Weakly informative: allows moderate to large effects",
       x = "Standardized Effect Size", y = "Density")

# =============================================================================
# Prior predictive simulation
# =============================================================================

# Simulate from priors only (no data influence)
# This helps check if priors generate sensible predictions

prior_predictive <- brm(
  bf(TTMz ~ BNRz * FPz) +
  bf(FS ~ BNRz + TTMz, family = bernoulli()),
  data = model_data,
  prior = priors,
  sample_prior = "only",  # KEY: Only sample from priors
  chains = 2, iter = 1000,
  seed = 42
)

# Check if prior predictions are sensible
pp_check(prior_predictive, resp = "TTMz", ndraws = 50) +
  labs(title = "Prior Predictive Check: Mediator (TTM)")

# =============================================================================
# COMMON PITFALLS
# =============================================================================

# PROBLEMATIC: Flat priors can cause computational issues
# prior(uniform(-Inf, Inf), class = "b")  # DON'T DO THIS

# BETTER: Weakly informative
# prior(normal(0, 1), class = "b")

# -----------------------------------------------------------------------------

# WRONG: brms doesn't know which equation this applies to
# prior(normal(0, 1), class = "b")  # Ambiguous in multivariate models!

# CORRECT: Specify response variable
# prior(normal(0, 1), class = "b", resp = "TTMz")
# prior(normal(0, 1), class = "b", resp = "FS")

# -----------------------------------------------------------------------------

# For logistic regression (binary Y), consider tighter priors:
logistic_priors <- c(
  prior(normal(0, 1), class = "b", resp = "TTMz"),
  prior(exponential(1), class = "sigma", resp = "TTMz"),
  # Tighter for logistic: Normal(0, 1.5) on log-odds scale
  # This says "effects larger than ~3 log-odds are very unlikely"
  prior(normal(0, 1.5), class = "b", resp = "FS")
)

# =============================================================================
# CRITICAL: The standardization-prior interaction
# =============================================================================

# DISASTER: Applying Normal(0, 1.5) to unstandardized income
# If income is measured in dollars (mean = 50,000, SD = 20,000):
# beta_income ~ Normal(0, 1.5)
# This prior says "I expect a $1 increase in income to change Y by
# at most ~3 units" -- essentially forcing the coefficient to zero!

# WHAT HAPPENED: You accidentally imposed an incredibly informative
# prior that overwhelms your data completely.

# THE FIX: Either standardize (Step 2) or adjust your prior:
# prior(normal(0, 0.0001), class = "b", coef = "income")
# OR just standardize your predictors -- it's easier and safer.

# =============================================================================
# Save priors for later use
# =============================================================================

saveRDS(priors, "priors_specified.rds")
cat("Priors saved to priors_specified.rds\n")
