# =============================================================================
# prior_sensitivity.R
# Prior Sensitivity Analysis for Bayesian Moderated Mediation
# =============================================================================

library(brms)
library(tidyverse)
library(posterior)

# =============================================================================
# Why sensitivity analysis matters
# =============================================================================

# If conclusions change across reasonable priors, your data aren't strong
# enough to overcome prior assumptions. This is valuable information!

# =============================================================================
# Define three prior specifications
# =============================================================================

# Default (used throughout this tutorial): moderately informative
priors_default <- c(
  prior(normal(0, 1), class = "b", resp = "M"),
  prior(normal(0, 1), class = "b", resp = "Y"),
  prior(exponential(1), class = "sigma")
)

# Tighter: stronger regularization (more skeptical)
priors_tight <- c(
  prior(normal(0, 0.5), class = "b", resp = "M"),
  prior(normal(0, 0.5), class = "b", resp = "Y"),
  prior(exponential(2), class = "sigma")
)

# Wider: more diffuse (less informative)
priors_wide <- c(
  prior(normal(0, 2), class = "b", resp = "M"),
  prior(normal(0, 2), class = "b", resp = "Y"),
  prior(exponential(0.5), class = "sigma")
)

# =============================================================================
# Example data
# =============================================================================

set.seed(42)
n <- 200

data <- tibble(
  X = rnorm(n),
  W = rnorm(n),
  M = 0.3*X + 0.2*W + 0.15*X*W + rnorm(n, 0, 0.85),
  Y = 0.4*M + 0.1*X + rnorm(n, 0, 0.8)
) %>%
  mutate(across(everything(), ~ scale(.x)[,1]))

# =============================================================================
# Fit all three models
# =============================================================================

cat("Fitting models with three different prior specifications...\n\n")

fit_default <- brm(
  bf(M ~ X * W) + bf(Y ~ X + M),
  data = data, prior = priors_default,
  chains = 4, iter = 3000, warmup = 1500,
  cores = 4, seed = 42, silent = 2, refresh = 0
)

fit_tight <- brm(
  bf(M ~ X * W) + bf(Y ~ X + M),
  data = data, prior = priors_tight,
  chains = 4, iter = 3000, warmup = 1500,
  cores = 4, seed = 42, silent = 2, refresh = 0
)

fit_wide <- brm(
  bf(M ~ X * W) + bf(Y ~ X + M),
  data = data, prior = priors_wide,
  chains = 4, iter = 3000, warmup = 1500,
  cores = 4, seed = 42, silent = 2, refresh = 0
)

# =============================================================================
# Extract IMM from each model
# =============================================================================

extract_IMM <- function(fit) {
  posts <- as_draws_df(fit)
  a3 <- posts$`b_M_X:W`
  b <- posts$b_Y_M
  IMM <- a3 * b
  return(IMM)
}

IMM_default <- extract_IMM(fit_default)
IMM_tight <- extract_IMM(fit_tight)
IMM_wide <- extract_IMM(fit_wide)

# =============================================================================
# Compare results
# =============================================================================

cat("=== Prior Sensitivity Analysis ===\n\n")

sensitivity <- tibble(
  Prior = c("Tight (N(0, 0.5))", "Default (N(0, 1))", "Wide (N(0, 2))"),
  IMM_median = c(median(IMM_tight), median(IMM_default), median(IMM_wide)),
  CI_lower = c(quantile(IMM_tight, 0.025), quantile(IMM_default, 0.025),
               quantile(IMM_wide, 0.025)),
  CI_upper = c(quantile(IMM_tight, 0.975), quantile(IMM_default, 0.975),
               quantile(IMM_wide, 0.975)),
  Excludes_zero = c(
    (quantile(IMM_tight, 0.025) > 0) | (quantile(IMM_tight, 0.975) < 0),
    (quantile(IMM_default, 0.025) > 0) | (quantile(IMM_default, 0.975) < 0),
    (quantile(IMM_wide, 0.025) > 0) | (quantile(IMM_wide, 0.975) < 0)
  )
)

print(sensitivity, digits = 3)

# =============================================================================
# Visualize
# =============================================================================

tibble(
  Prior = rep(c("Tight", "Default", "Wide"), each = length(IMM_default)),
  IMM = c(IMM_tight, IMM_default, IMM_wide)
) %>%
  mutate(Prior = factor(Prior, levels = c("Tight", "Default", "Wide"))) %>%
  ggplot(aes(x = IMM, fill = Prior)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Prior Sensitivity: IMM Posterior Distributions",
       subtitle = "Conclusions should be similar across reasonable priors",
       x = "Index of Moderated Mediation", y = "Density")

# =============================================================================
# Interpretation
# =============================================================================

cat("\n=== INTERPRETATION ===\n\n")

# Check if conclusions are consistent
all_same_conclusion <- length(unique(sensitivity$Excludes_zero)) == 1

if (all_same_conclusion) {
  cat("ROBUST: All prior specifications lead to the same conclusion.\n")
  cat("Your results are not sensitive to reasonable prior choices.\n")
} else {
  cat("*** WARNING: CONCLUSIONS DIFFER ACROSS PRIORS ***\n\n")
  cat("This means your data aren't strong enough to overcome prior assumptions.\n")
  cat("Report the sensitivity analysis explicitly and interpret cautiously.\n\n")

  # Example of what to report:
  cat("Example reporting:\n")
  cat("'With default priors (N(0,1)), IMM =", round(median(IMM_default), 3),
      "[", round(quantile(IMM_default, 0.025), 3), ",",
      round(quantile(IMM_default, 0.975), 3), "].\n")
  cat("However, conclusions were sensitive to prior specification:\n")
  cat("tighter priors yielded", ifelse(sensitivity$Excludes_zero[1], "significant", "non-significant"),
      "results.\n")
  cat("These findings should be interpreted with caution.'\n")
}

# =============================================================================
# Best practice recommendation
# =============================================================================

cat("\n=== BEST PRACTICE ===\n")
cat("Always run sensitivity analysis with at least two prior specifications.\n")
cat("Report whether conclusions are robust.\n")
cat("If they aren't, this is valuable information - not a failure.\n")

# Save results
saveRDS(list(
  sensitivity_table = sensitivity,
  IMM_default = IMM_default,
  IMM_tight = IMM_tight,
  IMM_wide = IMM_wide
), "prior_sensitivity_results.rds")
