# =============================================================================
# frequentist_comparison.R
# Frequentist Bootstrap Approach Using lavaan
# For comparison with Bayesian brms results
# =============================================================================

library(lavaan)
library(tidyverse)

# =============================================================================
# Example Data (same as Bayesian tutorial)
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

# Create interaction term (lavaan requires explicit interaction variable)
data$XW <- data$X * data$W

# =============================================================================
# Specify lavaan Model
# =============================================================================

# Model 8: First-stage moderated mediation
# M ~ a1*X + a2*W + a3*X:W
# Y ~ c'*X + b*M

model_syntax <- '
  # Mediator equation
  M ~ a1*X + a2*W + a3*XW

  # Outcome equation
  Y ~ cprime*X + b*M

  # Define indirect effects at different W values
  # Low W (-1 SD)
  indirect_low := (a1 + a3*(-1)) * b

  # Mean W (0)
  indirect_mean := a1 * b

  # High W (+1 SD)
  indirect_high := (a1 + a3*(1)) * b

  # Index of Moderated Mediation
  IMM := a3 * b
'

# =============================================================================
# Fit Model with Bootstrap
# =============================================================================

cat("Fitting lavaan model with 5000 bootstrap samples...\n")
cat("This may take a few minutes.\n\n")

fit_lavaan <- sem(
  model_syntax,
  data = data,
  se = "bootstrap",
  bootstrap = 5000,
  iseed = 42
)

# =============================================================================
# Extract Results
# =============================================================================

cat("=== FREQUENTIST RESULTS ===\n\n")

# Get parameter estimates with bootstrap CIs
params <- parameterEstimates(fit_lavaan, boot.ci.type = "perc")

# Filter to key effects
key_effects <- params %>%
  filter(label %in% c("a1", "a2", "a3", "b", "cprime",
                      "indirect_low", "indirect_mean", "indirect_high", "IMM"))

print(key_effects %>%
        select(label, est, se, ci.lower, ci.upper, pvalue) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# =============================================================================
# Focus on IMM
# =============================================================================

cat("\n=== INDEX OF MODERATED MEDIATION ===\n\n")

IMM_result <- params %>% filter(label == "IMM")

cat("IMM =", round(IMM_result$est, 4), "\n")
cat("SE =", round(IMM_result$se, 4), "\n")
cat("95% Bootstrap CI: [", round(IMM_result$ci.lower, 4), ",",
    round(IMM_result$ci.upper, 4), "]\n")
cat("p-value =", round(IMM_result$pvalue, 4), "\n\n")

IMM_sig <- (IMM_result$ci.lower > 0) | (IMM_result$ci.upper < 0)
cat("95% CI excludes zero:", IMM_sig, "\n")

# =============================================================================
# Compare with Bayesian (if available)
# =============================================================================

cat("\n=== COMPARISON NOTES ===\n")
cat("1. Bootstrap CI uses percentile method (recommended)\n")
cat("2. For small N (< 200), Bayesian credible intervals often have better coverage\n")
cat("3. Large N (> 500): Results should be nearly identical\n")
cat("4. Key difference: Bayesian provides P(IMM > 0), frequentist provides p-value\n")

# =============================================================================
# Model Fit Indices
# =============================================================================

cat("\n=== MODEL FIT ===\n")
fit_measures <- fitMeasures(fit_lavaan, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
print(round(fit_measures, 4))

# Note: For just-identified models (like this one), fit indices may not be informative

# =============================================================================
# Effect Size: Partially Standardized Indirect Effects
# =============================================================================

cat("\n=== STANDARDIZED EFFECTS ===\n")
std_params <- standardizedSolution(fit_lavaan)

std_indirect <- std_params %>%
  filter(label %in% c("indirect_low", "indirect_mean", "indirect_high", "IMM"))

print(std_indirect %>%
        select(label, est.std, se, ci.lower, ci.upper) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# =============================================================================
# Johnson-Neyman Regions (Manual Calculation)
# =============================================================================

cat("\n=== JOHNSON-NEYMAN ANALYSIS ===\n")
cat("Finding regions where indirect effect is significant...\n\n")

# Extract coefficients
a1 <- coef(fit_lavaan)["a1"]
a3 <- coef(fit_lavaan)["a3"]
b <- coef(fit_lavaan)["b"]

# Get variance-covariance matrix for derived parameters
# This is approximate for the product a3*b
vcov_mat <- vcov(fit_lavaan)

# Simplified JN calculation
# For exact JN regions with products, use the interactions package or custom bootstrap

w_range <- seq(-3, 3, by = 0.1)
indirect_at_w <- (a1 + a3 * w_range) * b

# Bootstrap SE at each W (simplified - use interactions package for exact)
cat("Indirect effect at different moderator values:\n")
for (w in c(-1, 0, 1)) {
  indirect <- params %>% filter(label == paste0("indirect_",
                                                 ifelse(w == -1, "low",
                                                        ifelse(w == 0, "mean", "high"))))
  cat(sprintf("  W = %d: %.4f [%.4f, %.4f]\n",
              w, indirect$est, indirect$ci.lower, indirect$ci.upper))
}

# =============================================================================
# Save Results
# =============================================================================

saveRDS(list(
  fit = fit_lavaan,
  params = params,
  IMM = IMM_result
), "frequentist_results.rds")

cat("\n=== COMPLETE ===\n")
cat("Results saved to frequentist_results.rds\n")

# =============================================================================
# Side-by-Side Template
# =============================================================================

cat("\n=== REPORTING TEMPLATE ===\n")
cat("
For manuscript reporting, include both analyses:

'We analyzed moderated mediation using both frequentist (lavaan with 5000
bootstrap samples) and Bayesian (brms with default priors) approaches.
The Index of Moderated Mediation was [ESTIMATE] with bootstrap 95% CI
[LOWER, UPPER] and Bayesian 95% credible interval [LOWER, UPPER]. Both
methods [agreed/disagreed] that the CI [included/excluded] zero, suggesting
[no evidence for/evidence supporting] moderated mediation.'
")
