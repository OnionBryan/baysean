# =============================================================================
# 06_test_moderation.R
# Bayesian Moderated Mediation Tutorial
# Step 7: Test Moderation (Simple Slopes & Johnson-Neyman)
# =============================================================================

library(brms)
library(posterior)
library(tidyverse)

# Load fitted model
fit <- readRDS("fitted_model.rds")
posts <- as_draws_df(fit)

# =============================================================================
# Test interaction effect using hypothesis()
# =============================================================================

cat("=== Testing Moderation ===\n\n")

# Bayesian hypothesis test for interaction
hyp_test <- hypothesis(fit, "TTMz_BNR_z:FP_z = 0")
print(hyp_test)

# Manual calculation with posterior
a3 <- posts$`b_TTMz_BNR_z:FP_z`
cat("\nInteraction (a3):\n")
cat("  Median:", round(median(a3), 3), "\n")
cat("  95% CI: [", round(quantile(a3, 0.025), 3), ",",
    round(quantile(a3, 0.975), 3), "]\n")
cat("  P(a3 > 0):", round(mean(a3 > 0), 3), "\n")

# =============================================================================
# Simple slopes: Effect of X on M at different W levels
# =============================================================================

a1 <- posts$b_TTMz_BNR_z

# Effect of BNR on TTM at different FP levels
W_values <- c(-1, 0, 1)  # Low, Mean, High (standardized)

simple_slopes <- tibble(
  W_level = c("Low (-1 SD)", "Mean (0)", "High (+1 SD)"),
  W_value = W_values,
  slope = map_dbl(W_values, ~ median(a1 + a3 * .x)),
  CI_lower = map_dbl(W_values, ~ quantile(a1 + a3 * .x, 0.025)),
  CI_upper = map_dbl(W_values, ~ quantile(a1 + a3 * .x, 0.975)),
  Prob_positive = map_dbl(W_values, ~ mean((a1 + a3 * .x) > 0))
)

cat("\n=== Simple Slopes ===\n")
cat("Effect of X on M at different levels of W:\n\n")
print(simple_slopes, digits = 3)

# =============================================================================
# Visualize simple slopes
# =============================================================================

# Create data for plotting
plot_data <- tibble(
  W = rep(W_values, each = length(a1)),
  W_label = rep(c("Low (-1 SD)", "Mean", "High (+1 SD)"), each = length(a1)),
  slope = c(a1 + a3 * (-1), a1 + a3 * 0, a1 + a3 * 1)
)

ggplot(plot_data, aes(x = slope, fill = W_label)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Simple Slopes: Effect of X on M at Different W Levels",
       subtitle = "Distribution of slopes; vertical line at zero",
       x = "Slope (effect of X on M)", y = "Density",
       fill = "Moderator (W)")

# =============================================================================
# Johnson-Neyman: Region of significance
# =============================================================================

# At what value of W does the X effect become significant?

# Calculate effect of X at many W values
W_range <- seq(-3, 3, by = 0.1)

jn_results <- tibble(
  W = W_range,
  slope = map_dbl(W_range, ~ median(a1 + a3 * .x)),
  CI_lower = map_dbl(W_range, ~ quantile(a1 + a3 * .x, 0.025)),
  CI_upper = map_dbl(W_range, ~ quantile(a1 + a3 * .x, 0.975)),
  significant = CI_lower > 0 | CI_upper < 0
)

# Find transition points
jn_transitions <- jn_results %>%
  mutate(sig_change = significant != lag(significant)) %>%
  filter(sig_change == TRUE) %>%
  pull(W)

cat("\n=== Johnson-Neyman Analysis ===\n")
cat("Transition points (where effect becomes significant):\n")
print(jn_transitions)

# Plot Johnson-Neyman
ggplot(jn_results, aes(x = W, y = slope)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.3, fill = "steelblue") +
  geom_line(size = 1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = jn_transitions, linetype = "dotted", color = "red") +
  theme_minimal() +
  labs(title = "Johnson-Neyman Plot",
       subtitle = "Shaded region = 95% CI; red lines = significance transitions",
       x = "Moderator (W)", y = "Effect of X on M")

# =============================================================================
# COMMON PITFALLS
# =============================================================================

# If a3 CI includes zero:
# "The effect was stronger at high FP than low FP"  # WRONG

# CORRECT interpretation:
# "We did not find evidence of moderation (95% CI includes zero)"

# -----------------------------------------------------------------------------

# QUESTIONABLE: Probing when a3 CI includes zero
# The differences in slopes could be due to chance

# RECOMMENDED: Only probe if:
# 1. a3 CI excludes zero, OR
# 2. You have theoretical reasons and note exploratory nature

# -----------------------------------------------------------------------------

# If FP was measured 0-100 (not standardized):
# "Low" = 0, "High" = 100 may not be meaningful

# SOLUTION: Always standardize or use meaningful values
# e.g., 25th percentile, median, 75th percentile

# =============================================================================
# Save results
# =============================================================================

moderation_results <- list(
  hypothesis_test = hyp_test,
  simple_slopes = simple_slopes,
  jn_results = jn_results,
  jn_transitions = jn_transitions
)

saveRDS(moderation_results, "moderation_results.rds")
cat("\nModeration results saved to moderation_results.rds\n")
