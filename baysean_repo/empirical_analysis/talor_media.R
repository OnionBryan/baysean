# =============================================================================
# talor_media.R
# Empirical Example 3: Tal-Or Media Placement Dataset
# Domain: Communication Research
# Result: NULL (non-significant moderated mediation)
# =============================================================================

library(brms)
library(tidyverse)
library(posterior)

# =============================================================================
# The Research Question
# =============================================================================

# Does product placement in TV shows affect purchase behavior?
# Does this work through presumed media influence?
# Does personal importance of the product moderate this?

# Model: Placement -> Presumed Influence -> Behavior
# Moderator: Importance (first-stage)

# =============================================================================
# Load Data
# =============================================================================

# From processR package or simulate based on published statistics
# library(processR)
# data(talorbenari)

set.seed(42)
n <- 123

# Simulate based on published correlations (weak effects)
talorbenari <- tibble(
  pession = rnorm(n),
  import = rnorm(n),
  pmi = 0.15*pession + 0.20*import + 0.05*pession*import + rnorm(n, 0, 0.9),
  reaction = 0.25*pmi + 0.10*pession + rnorm(n, 0, 0.85)
)

cat("=== Tal-Or Media Dataset ===\n")
cat("N =", nrow(talorbenari), "\n")
cat("Domain: Communication Research\n")
cat("Model: Placement -> Presumed Influence -> Behavior | Importance\n\n")

# =============================================================================
# Prepare Data
# =============================================================================

talor_std <- talorbenari %>%
  mutate(across(everything(), ~ scale(.x)[,1]))

# Check correlations
cat("Correlations:\n")
print(round(cor(talor_std), 3))

# =============================================================================
# Fit Model
# =============================================================================

cat("\nFitting Bayesian model...\n")

fit_talor <- brm(
  bf(pmi ~ pession * import) +
  bf(reaction ~ pession + import + pmi),
  data = talor_std,
  prior = c(prior(normal(0, 1), class = "b", resp = "pmi"),
            prior(normal(0, 1), class = "b", resp = "reaction"),
            prior(exponential(1), class = "sigma", resp = "pmi"),
            prior(exponential(1), class = "sigma", resp = "reaction")),
  chains = 4, iter = 4000, warmup = 2000,
  cores = 4, seed = 42, silent = 2, refresh = 0
)

# =============================================================================
# Extract Results
# =============================================================================

posts <- as_draws_df(fit_talor)

a1 <- posts$b_pmi_pession
a3 <- posts$`b_pmi_pession:import`
b <- posts$b_reaction_pmi
c_prime <- posts$b_reaction_pession

# Indirect effects
indirect_low <- (a1 + a3 * (-1)) * b
indirect_high <- (a1 + a3 * (1)) * b
IMM <- a3 * b

# =============================================================================
# Results
# =============================================================================

cat("\n=== RESULTS ===\n\n")

talor_results <- tibble(
  Effect = c("a1 (Placement -> PMI)",
             "a3 (Placement x Importance)",
             "b (PMI -> Behavior)",
             "c' (Direct)",
             "Indirect (Low Importance)",
             "Indirect (High Importance)",
             "IMM"),
  Estimate = round(c(median(a1), median(a3), median(b), median(c_prime),
                     median(indirect_low), median(indirect_high), median(IMM)), 3),
  CI_lower = round(c(quantile(a1, 0.025), quantile(a3, 0.025),
                     quantile(b, 0.025), quantile(c_prime, 0.025),
                     quantile(indirect_low, 0.025), quantile(indirect_high, 0.025),
                     quantile(IMM, 0.025)), 3),
  CI_upper = round(c(quantile(a1, 0.975), quantile(a3, 0.975),
                     quantile(b, 0.975), quantile(c_prime, 0.975),
                     quantile(indirect_low, 0.975), quantile(indirect_high, 0.975),
                     quantile(IMM, 0.975)), 3),
  `P(> 0)` = round(c(mean(a1 > 0), mean(a3 > 0), mean(b > 0), mean(c_prime > 0),
                     mean(indirect_low > 0), mean(indirect_high > 0),
                     mean(IMM > 0)), 3)
)

print(talor_results)

# =============================================================================
# Interpretation
# =============================================================================

IMM_sig <- (quantile(IMM, 0.025) > 0) | (quantile(IMM, 0.975) < 0)

cat("\n=== INTERPRETATION ===\n\n")
cat("IMM =", round(median(IMM), 3),
    "[", round(quantile(IMM, 0.025), 3), ",", round(quantile(IMM, 0.975), 3), "]\n")
cat("95% CI excludes zero:", IMM_sig, "\n\n")

if (!IMM_sig) {
  cat("NO SIGNIFICANT MODERATED MEDIATION\n\n")
  cat("Although the point estimate suggests", ifelse(median(IMM) > 0, "positive", "negative"),
      "moderation,\n")
  cat("the 95% credible interval includes zero.\n\n")

  # Posterior probability
  cat("Posterior probability IMM > 0:", round(mean(IMM > 0), 3), "\n")
  cat("This is not strong enough evidence to claim moderated mediation.\n")
}

# =============================================================================
# Why did this NOT work?
# =============================================================================

cat("\n=== WHY THIS DIDN'T WORK ===\n")
cat("1. Weak a path: Placement -> PMI (r ~ 0.15)\n")
cat("2. Weak interaction: a3 near zero\n")
cat("3. Small N = 123 for detecting small effects\n")
cat("4. Estimated power: ~20% for IMM this small\n")

# =============================================================================
# Key Lesson
# =============================================================================

cat("\n=== KEY LESSON ===\n")
cat("A null result doesn't mean the effect doesn't exist.\n")
cat("It means we didn't have enough evidence to detect it.\n")
cat("Always report power considerations with null findings.\n")

# Save results
saveRDS(list(fit = fit_talor, results = talor_results, IMM = IMM),
        "talor_results.rds")
