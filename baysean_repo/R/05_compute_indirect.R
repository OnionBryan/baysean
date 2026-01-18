# =============================================================================
# 05_compute_indirect.R
# Bayesian Moderated Mediation Tutorial
# Step 6: Compute Indirect Effects
# =============================================================================

library(brms)
library(posterior)
library(tidyverse)

# Load fitted model
fit <- readRDS("fitted_model.rds")

# =============================================================================
# Extract posterior samples
# =============================================================================

posts <- as_draws_df(fit)

# --- Extract path coefficients ---
a1 <- posts$b_TTMz_BNR_z           # X -> M (main effect)
a2 <- posts$b_TTMz_FP_z            # W -> M (moderator main effect)
a3 <- posts$`b_TTMz_BNR_z:FP_z`    # X*W -> M (interaction)
b <- posts$b_FS_TTM_z              # M -> Y
c_prime <- posts$b_FS_BNR_z        # X -> Y (direct effect)

# =============================================================================
# Compute indirect effects at different levels of moderator
# =============================================================================

# For first-stage moderation: indirect = (a1 + a3*W) * b
# At different values of W (standardized):
W_low <- -1   # 1 SD below mean
W_mean <- 0   # At mean
W_high <- 1   # 1 SD above mean

# Conditional indirect effects
indirect_low <- (a1 + a3 * W_low) * b
indirect_mean <- (a1 + a3 * W_mean) * b
indirect_high <- (a1 + a3 * W_high) * b

# =============================================================================
# Index of Moderated Mediation (IMM)
# =============================================================================

# IMM = a3 * b
# This is the KEY quantity: how much does the indirect effect change
# for a 1-unit increase in the moderator?

IMM <- a3 * b

# =============================================================================
# Total effect
# =============================================================================

total_effect <- (a1 + a3 * W_mean) * b + c_prime

# =============================================================================
# Create summary table
# =============================================================================

mediation_summary <- tibble(
  Effect = c("a₁ (X → M)", "a₃ (X×W → M)", "b (M → Y)", "c' (direct)",
             "Indirect (W = -1 SD)", "Indirect (W = 0)",
             "Indirect (W = +1 SD)", "IMM"),
  Estimate = c(median(a1), median(a3), median(b), median(c_prime),
               median(indirect_low), median(indirect_mean),
               median(indirect_high), median(IMM)),
  CI_lower = c(quantile(a1, 0.025), quantile(a3, 0.025),
               quantile(b, 0.025), quantile(c_prime, 0.025),
               quantile(indirect_low, 0.025), quantile(indirect_mean, 0.025),
               quantile(indirect_high, 0.025), quantile(IMM, 0.025)),
  CI_upper = c(quantile(a1, 0.975), quantile(a3, 0.975),
               quantile(b, 0.975), quantile(c_prime, 0.975),
               quantile(indirect_low, 0.975), quantile(indirect_mean, 0.975),
               quantile(indirect_high, 0.975), quantile(IMM, 0.975)),
  Prob_positive = c(mean(a1 > 0), mean(a3 > 0), mean(b > 0), mean(c_prime > 0),
                    mean(indirect_low > 0), mean(indirect_mean > 0),
                    mean(indirect_high > 0), mean(IMM > 0))
)

cat("=== Mediation Summary ===\n\n")
print(mediation_summary, digits = 3)

# Check if IMM is significant (95% CI excludes zero)
IMM_significant <- (quantile(IMM, 0.025) > 0) | (quantile(IMM, 0.975) < 0)
cat("\n*** IMM 95% CI excludes zero:", IMM_significant, "***\n")

# =============================================================================
# Visualize indirect effect distribution
# =============================================================================

tibble(
  `W = -1 SD` = indirect_low,
  `W = 0` = indirect_mean,
  `W = +1 SD` = indirect_high
) %>%
  pivot_longer(everything(), names_to = "Moderator Level", values_to = "Indirect") %>%
  ggplot(aes(x = Indirect, fill = `Moderator Level`)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Posterior Distribution of Conditional Indirect Effects",
       subtitle = "Vertical line at zero; density away from zero = evidence of effect",
       x = "Indirect Effect", y = "Density")

# =============================================================================
# IMPORTANT: First-stage vs Second-stage formulas
# =============================================================================

# FIRST-STAGE moderation (W moderates a path):
# indirect(W) = (a1 + a3*W) * b
# IMM = a3 * b

# SECOND-STAGE moderation (W moderates b path):
# indirect(W) = a * (b1 + b3*W)
# IMM = a * b3

# =============================================================================
# For binary outcomes: interpretation note
# =============================================================================

# b coefficient is on LOG-ODDS scale
# b = -0.5 means: 1 unit increase in M decreases log-odds of Y by 0.5

# To get odds ratio:
OR <- exp(median(b))
cat("\nOdds ratio for b path:", round(OR, 3), "\n")

# To get probability change: depends on baseline probability
# Use marginal effects for interpretable quantities

# =============================================================================
# COMMON PITFALL: Not propagating uncertainty
# =============================================================================

# WRONG: Using point estimates (ignores uncertainty)
# IMM_wrong <- median(a3) * median(b)

# CORRECT: Multiply full posteriors (preserves uncertainty)
# IMM_correct <- a3 * b  # This is a vector of 4000+ values

# =============================================================================
# Save results
# =============================================================================

results <- list(
  summary = mediation_summary,
  IMM_posterior = IMM,
  indirect_low = indirect_low,
  indirect_mean = indirect_mean,
  indirect_high = indirect_high,
  IMM_significant = IMM_significant
)

saveRDS(results, "indirect_effects_results.rds")
cat("\nResults saved to indirect_effects_results.rds\n")
