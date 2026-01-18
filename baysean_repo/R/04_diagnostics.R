# =============================================================================
# 04_diagnostics.R
# Bayesian Moderated Mediation Tutorial
# Step 5: Check Diagnostics
# =============================================================================

library(brms)
library(bayesplot)
library(posterior)
library(tidyverse)

# Load fitted model
fit <- readRDS("fitted_model.rds")

# =============================================================================
# Check Rhat and ESS for all parameters
# =============================================================================

cat("=== Convergence Diagnostics ===\n\n")

# Get summary with diagnostics
model_summary <- summary(fit)

# Check Rhat (should be < 1.01 for all parameters)
cat("Rhat values (should be < 1.01):\n")
print(model_summary$fixed[, c("Estimate", "Rhat", "Bulk_ESS", "Tail_ESS")])

# Quick diagnostic check
rhats <- rhat(fit)
cat("\nMax Rhat:", max(rhats, na.rm = TRUE), "\n")
cat("All Rhat < 1.01:", all(rhats < 1.01, na.rm = TRUE), "\n")

# =============================================================================
# Trace plots for key parameters
# =============================================================================

# Visual inspection of chain mixing
mcmc_trace(fit, pars = c("b_TTMz_BNR_z", "b_TTMz_FP_z", "b_TTMz_BNR_z:FP_z",
                          "b_FS_TTM_z")) +
  labs(title = "Trace Plots: Key Parameters",
       subtitle = "Chains should look like 'fuzzy caterpillars' - well-mixed")

# =============================================================================
# Posterior predictive check for mediator
# =============================================================================

pp_check(fit, resp = "TTMz", ndraws = 100) +
  labs(title = "Posterior Predictive Check: Mediator (TTM)",
       subtitle = "Dark line = observed data; light lines = model predictions")

# =============================================================================
# What good vs bad traces look like
# =============================================================================

# GOOD trace: "fuzzy caterpillar"
# - Chains overlap completely
# - No trends or drifts
# - Random fluctuation around stable mean

# BAD trace patterns:
# - Chains not overlapping = different chains found different modes
# - Trending up/down = hasn't converged yet
# - Stuck periods = sampler got stuck

# =============================================================================
# ESS warnings
# =============================================================================

# Check effective sample size
ess_bulk <- model_summary$fixed[, "Bulk_ESS"]
ess_tail <- model_summary$fixed[, "Tail_ESS"]

if (any(ess_bulk < 400) | any(ess_tail < 400)) {
  cat("\nWARNING: Some ESS values are low (<400)\n")
  cat("Consider running more iterations.\n")

  # SOLUTION: Run more iterations
  # fit <- update(fit, iter = 6000, warmup = 3000)

  # OR reduce model complexity
  # OR improve parameterization
}

# =============================================================================
# Complete diagnostic check
# =============================================================================

# BAD: Only checking Rhat
# if (max(rhat(fit)) < 1.01) { print("Converged!") }

# GOOD: Check both Rhat AND ESS
check_convergence <- function(fit) {
  rhats <- rhat(fit)
  ess <- neff_ratio(fit)

  passed <- TRUE

  if (max(rhats, na.rm = TRUE) >= 1.01) {
    cat("FAIL: Rhat >= 1.01 for some parameters\n")
    passed <- FALSE
  }

  if (min(ess, na.rm = TRUE) < 0.1) {
    cat("FAIL: ESS ratio < 0.1 for some parameters\n")
    passed <- FALSE
  }

  if (passed) {
    cat("PASS: All convergence diagnostics OK\n")
  }

  return(passed)
}

check_convergence(fit)

# =============================================================================
# Additional diagnostics
# =============================================================================

# Pairs plot to check for problematic correlations
mcmc_pairs(fit, pars = c("b_TTMz_BNR_z", "b_TTMz_BNR_z:FP_z"))

# Energy diagnostic (Stan-specific)
mcmc_nuts_energy(nuts_params(fit))

cat("\nDiagnostics complete.\n")
