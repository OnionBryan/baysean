# =============================================================================
# 07_reporting.R
# Bayesian Moderated Mediation Tutorial
# Step 8: Report Results (Publication-Ready Tables)
# =============================================================================

library(tidyverse)
library(knitr)
library(kableExtra)

# Load results
results <- readRDS("indirect_effects_results.rds")
moderation <- readRDS("moderation_results.rds")

# =============================================================================
# Create publication-ready table
# =============================================================================

# Format results for publication
pub_table <- results$summary %>%
  mutate(
    # Format estimate with CI
    `Estimate [95% CI]` = sprintf("%.3f [%.3f, %.3f]",
                                   Estimate, CI_lower, CI_upper),
    # Format posterior probability
    `P(> 0)` = sprintf("%.3f", Prob_positive)
  ) %>%
  select(Effect, `Estimate [95% CI]`, `P(> 0)`)

cat("=== Publication Table ===\n\n")
print(pub_table)

# For LaTeX output
kable(pub_table, format = "latex", booktabs = TRUE,
      caption = "Bayesian Moderated Mediation Results") %>%
  kable_styling(latex_options = "hold_position")

# =============================================================================
# Interpretation guidelines
# =============================================================================

# WRONG: Bayesian analysis doesn't produce p-values
# "The indirect effect was significant (p < .05)"

# CORRECT: Report credible intervals and posterior probabilities
# "The indirect effect at high W was 0.15, 95% CrI [0.03, 0.28],
#  with 98% posterior probability of being positive."

# -----------------------------------------------------------------------------

# FREQUENTIST CI: "If we repeated this study infinitely,
# 95% of CIs would contain the true value"

# BAYESIAN CrI: "Given the data and priors, there is a
# 95% probability the parameter lies in this interval"

# The Bayesian interpretation is more intuitive!

# =============================================================================
# What to report in methods section
# =============================================================================

methods_template <- "
## Statistical Analysis

We conducted Bayesian moderated mediation analysis using brms version X.X
(Bürkner, 2017) in R version X.X.X. The model estimated first-stage moderation
(Hayes Model 8) where W moderates the X → M path.

**Priors**: We specified weakly informative Normal(0, 1) priors for all
regression coefficients on standardized variables, and Exponential(1) priors
for residual standard deviations. These priors allow moderate to large effects
while providing regularization.
[Note: Also cite Gelman et al. (2008) for prior justification]

**MCMC Settings**: We ran 4 chains of 4,000 iterations each (2,000 warmup),
yielding 8,000 posterior samples. Convergence was assessed via Rhat (all < 1.01)
and effective sample size (all ESS > 400).
[Report actual values]

**Inference**: We computed the Index of Moderated Mediation (IMM = a₃ × b)
and conditional indirect effects at ±1 SD of the moderator. We report
posterior medians with 95% equal-tailed credible intervals. Effects were
considered credible if the 95% CI excluded zero.
"

cat("\n=== Methods Template ===\n")
cat(methods_template)

# =============================================================================
# Results paragraph template
# =============================================================================

# Fill in with actual values
IMM_median <- median(results$IMM_posterior)
IMM_lower <- quantile(results$IMM_posterior, 0.025)
IMM_upper <- quantile(results$IMM_posterior, 0.975)
prob_positive <- mean(results$IMM_posterior > 0)

results_template <- sprintf("
## Results

The Index of Moderated Mediation was %.3f, 95%% CrI [%.3f, %.3f], with
%.1f%% posterior probability of being positive. This indicates that the
indirect effect of X on Y through M %s as W increases.

[If significant:]
At low W (-1 SD), the indirect effect was [value], 95%% CrI [lower, upper].
At high W (+1 SD), the indirect effect was [value], 95%% CrI [lower, upper].
The difference in indirect effects was [value], 95%% CrI [lower, upper].

[If non-significant:]
Although the point estimate suggested [direction], the 95%% credible interval
included zero, indicating insufficient evidence for moderated mediation.
",
IMM_median, IMM_lower, IMM_upper, prob_positive * 100,
ifelse(IMM_median > 0, "strengthens", "weakens"))

cat("\n=== Results Template ===\n")
cat(results_template)

# =============================================================================
# APA-style table function
# =============================================================================

create_apa_table <- function(results_df, caption = "Results") {
  results_df %>%
    kable(format = "latex", booktabs = TRUE, caption = caption,
          align = c("l", rep("c", ncol(results_df) - 1))) %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    footnote(general = "CrI = Credible Interval; P(> 0) = Posterior probability of positive effect.",
             general_title = "Note.",
             footnote_as_chunk = TRUE)
}

# Example usage
# create_apa_table(pub_table, "Bayesian Moderated Mediation Results")

cat("\n\nReporting complete.\n")
