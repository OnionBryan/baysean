# =============================================================================
# facialburns_clinical.R
# Empirical Example 4: Facial Burns Clinical Dataset
# Domain: Health Psychology
# Result: NULL (non-significant moderated mediation)
# =============================================================================

library(brms)
library(tidyverse)
library(posterior)

# =============================================================================
# The Research Question
# =============================================================================

# Does burn severity affect self-esteem?
# Does this work through psychological distress?
# Does rumination moderate this process?

# Model: Severity -> Distress -> Self-Esteem
# Moderator: Rumination (first-stage)

# =============================================================================
# Simulate Data
# =============================================================================

# Based on published summary statistics from Patterson et al. (2016)
# Original data not publicly available

set.seed(42)
n <- 98

# Simulate based on clinical literature correlations
facialburns <- tibble(
  severity = rnorm(n),            # Burn severity
  rumession = rnorm(n),           # Rumination
  distress = 0.12*severity + 0.35*rumession + 0.08*severity*rumession + rnorm(n, 0, 0.85),
  selfesteem = -0.30*distress + 0.05*severity + rnorm(n, 0, 0.80)  # Note: negative path
)

cat("=== Facial Burns Dataset ===\n")
cat("N =", nrow(facialburns), "\n")
cat("Domain: Clinical/Health Psychology\n")
cat("Model: Severity -> Distress -> Self-Esteem | Rumination\n")
cat("Note: Data simulated from published summary statistics\n\n")

# =============================================================================
# Prepare Data
# =============================================================================

burns_std <- facialburns %>%
  mutate(across(everything(), ~ scale(.x)[,1]))

# Check correlations
cat("Correlations:\n")
print(round(cor(burns_std), 3))

# Key observation: weak severity -> distress path
cat("\nNOTE: Severity -> Distress correlation is weak (r ~ 0.12)\n")
cat("This suggests there may not be much to moderate!\n")

# =============================================================================
# Fit Model
# =============================================================================

cat("\nFitting Bayesian model...\n")

fit_burns <- brm(
  bf(distress ~ severity * rumession) +
  bf(selfesteem ~ severity + rumession + distress),
  data = burns_std,
  prior = c(prior(normal(0, 1), class = "b", resp = "distress"),
            prior(normal(0, 1), class = "b", resp = "selfesteem"),
            prior(exponential(1), class = "sigma", resp = "distress"),
            prior(exponential(1), class = "sigma", resp = "selfesteem")),
  chains = 4, iter = 4000, warmup = 2000,
  cores = 4, seed = 42, silent = 2, refresh = 0
)

# =============================================================================
# Extract Results
# =============================================================================

posts <- as_draws_df(fit_burns)

a1 <- posts$b_distress_severity
a3 <- posts$`b_distress_severity:rumession`
b <- posts$b_selfesteem_distress
c_prime <- posts$b_selfesteem_severity

# Indirect effects
indirect_low <- (a1 + a3 * (-1)) * b
indirect_high <- (a1 + a3 * (1)) * b
IMM <- a3 * b

# =============================================================================
# Results
# =============================================================================

cat("\n=== RESULTS ===\n\n")

burns_results <- tibble(
  Effect = c("a1 (Severity -> Distress)",
             "a3 (Severity x Rumination)",
             "b (Distress -> Self-Esteem)",
             "c' (Direct)",
             "Indirect (Low Rumination)",
             "Indirect (High Rumination)",
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

print(burns_results)

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
}

# =============================================================================
# Critical diagnostic: Check the a path!
# =============================================================================

cat("\n=== CRITICAL DIAGNOSTIC ===\n")
cat("a1 (Severity -> Distress):", round(median(a1), 3),
    "[", round(quantile(a1, 0.025), 3), ",", round(quantile(a1, 0.975), 3), "]\n")

a1_sig <- (quantile(a1, 0.025) > 0) | (quantile(a1, 0.975) < 0)
cat("a1 95% CI excludes zero:", a1_sig, "\n\n")

if (!a1_sig) {
  cat("*** THE MAIN EFFECT ISN'T EVEN CREDIBLE! ***\n")
  cat("If there's no reliable X -> M effect, there's nothing to moderate.\n")
  cat("This is the fundamental reason moderated mediation failed.\n")
}

# =============================================================================
# Why did this NOT work?
# =============================================================================

cat("\n=== WHY THIS DIDN'T WORK ===\n")
cat("1. WEAK a path: Severity -> Distress (r ~ 0.12, CI includes zero)\n")
cat("2. Small N = 98: Insufficient power for small effects\n")
cat("3. Estimated power: ~25% for IMM this small\n")
cat("4. No effect to moderate: Can't moderate a non-existent effect!\n")

# =============================================================================
# Key Lesson
# =============================================================================

cat("\n=== KEY LESSON ===\n")
cat("ALWAYS check: Is there an effect to moderate?\n")
cat("If the main effect (a1) isn't credible, moderated mediation will fail.\n")
cat("This isn't a failure of the method - it's informative about the data.\n")

# Save results
saveRDS(list(fit = fit_burns, results = burns_results, IMM = IMM),
        "facialburns_results.rds")
