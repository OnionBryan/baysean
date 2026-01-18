# =============================================================================
# 03_fit_model.R
# Bayesian Moderated Mediation Tutorial
# Step 4: Fit Models
# =============================================================================

library(brms)
library(tidyverse)

# Load prepared data and priors
model_data <- readRDS("model_data_prepared.rds")
priors <- readRDS("priors_specified.rds")

# =============================================================================
# Fit the multivariate mediation model
# =============================================================================

# This may take 2-5 minutes depending on your computer
cat("Fitting Bayesian moderated mediation model...\n")
cat("This may take 2-5 minutes.\n\n")

fit <- brm(
  # Multivariate formula:
  # 1. Mediator equation: TTM predicted by BNR, FP, and their interaction
  # 2. Outcome equation: FS predicted by BNR, FP, and TTM
  bf(TTM_z ~ BNR_z * FP_z) +
  bf(FS ~ BNR_z + FP_z + TTM_z, family = bernoulli()),

  data = model_data,
  prior = priors,

  # MCMC settings
  chains = 4,        # Run 4 independent chains

  iter = 4000,       # Total iterations per chain
  warmup = 2000,     # Discard first 2000 as warmup
  cores = 4,         # Parallel processing
  seed = 42,         # Reproducibility

  # Sampler settings
  control = list(
    adapt_delta = 0.95,  # Increase if divergent transitions
    max_treedepth = 12   # Increase if hitting max treedepth
  ),

  # For cleaner output during fitting
  silent = 2,
  refresh = 500  # Print progress every 500 iterations
)

cat("\nModel fitting complete!\n\n")

# =============================================================================
# View model summary
# =============================================================================

cat("=== Model Summary ===\n\n")
summary(fit)

# --- Extract fixed effects ---
cat("\n=== Fixed Effects ===\n\n")
fixef(fit)

# =============================================================================
# COMMON PITFALLS
# =============================================================================

# --- Iterations ---

# TOO FEW: May not converge
# iter = 1000, warmup = 500  # Risky

# ADEQUATE: Standard recommendation
# iter = 2000, warmup = 1000  # Minimum for publication

# FOR FINAL ANALYSIS: More samples for stable estimates
# iter = 4000, warmup = 2000  # Recommended

# --- Divergent transitions ---

# If you see: "X divergent transitions after warmup"
# This indicates the sampler had difficulty exploring the posterior

# SOLUTION 1: Increase adapt_delta
# control = list(adapt_delta = 0.99)

# SOLUTION 2: Reparameterize or simplify model
# SOLUTION 3: Check data for outliers

# --- Residual correlation (set_rescor) ---

# WRONG: Can't estimate residual correlation between
# Gaussian (TTM) and Bernoulli (FS) families
# set_rescor(TRUE)  # Will error with mixed families

# CORRECT: Disable residual correlation for mixed families
# brms automatically handles this, but be aware

# =============================================================================
# Save fitted model
# =============================================================================

saveRDS(fit, "fitted_model.rds")
cat("Model saved to fitted_model.rds\n")
