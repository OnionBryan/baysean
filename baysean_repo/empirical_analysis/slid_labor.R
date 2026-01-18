# =============================================================================
# slid_labor.R
# Empirical Example 2: SLID-inspired Labor Economics Dataset
# Domain: Labor Economics
# Result: SIGNIFICANT moderated mediation
# =============================================================================

library(brms)
library(tidyverse)
library(posterior)

# =============================================================================
# The Research Question
# =============================================================================

# Does age affect wages through education?
# And does this indirect effect differ by sex?

# Model: Age -> Education -> Wages
# Moderator: Sex (first-stage)

# =============================================================================
# Simulate Data
# =============================================================================

# Based on patterns from Survey of Labour and Income Dynamics
# Simulated to ensure detectable moderated mediation

set.seed(42)
n <- 500

# Sex moderates how age relates to education
# (e.g., older cohorts had different educational patterns by sex)
slid_data <- tibble(
  age = rnorm(n),
  sex = rbinom(n, 1, 0.5),  # 0 = Female, 1 = Male
  education = -0.20*age + 0.10*sex + 0.30*age*sex + rnorm(n, 0, 0.80),
  wages = 0.45*education + 0.20*age + 0.10*sex + rnorm(n, 0, 0.70)
)

cat("=== SLID-inspired Dataset ===\n")
cat("N =", nrow(slid_data), "\n")
cat("Domain: Labor Economics\n")
cat("Model: Age -> Education -> Wages | Sex\n")
cat("Note: Data simulated based on SLID patterns\n\n")

# =============================================================================
# Prepare Data
# =============================================================================

slid_std <- slid_data %>%
  mutate(across(c(age, education, wages), ~ scale(.x)[,1]))

# Check correlations
cat("Correlations:\n")
print(round(cor(slid_std[, c("age", "education", "wages")]), 3))

# Check moderation
mod_check <- lm(education ~ age * sex, data = slid_std)
cat("\nDoes sex moderate age -> education?\n")
cat("Interaction coefficient:", round(coef(mod_check)["age:sex"], 3), "\n")
cat("Interaction p-value:", round(summary(mod_check)$coefficients["age:sex", "Pr(>|t|)"], 4), "\n")

# =============================================================================
# Fit Model
# =============================================================================

cat("\nFitting Bayesian model...\n")

fit_slid <- brm(
  bf(education ~ age * sex) +
  bf(wages ~ age + sex + education),
  data = slid_std,
  prior = c(prior(normal(0, 1), class = "b", resp = "education"),
            prior(normal(0, 1), class = "b", resp = "wages"),
            prior(exponential(1), class = "sigma", resp = "education"),
            prior(exponential(1), class = "sigma", resp = "wages")),
  chains = 4, iter = 4000, warmup = 2000,
  cores = 4, seed = 42,
  control = list(adapt_delta = 0.95),
  silent = 2, refresh = 0
)

# =============================================================================
# Extract Results
# =============================================================================

posts <- as_draws_df(fit_slid)

a1 <- posts$b_education_age
a3 <- posts$`b_education_age:sex`
b <- posts$b_wages_education
c_prime <- posts$b_wages_age

# Conditional indirect effects
indirect_female <- a1 * b               # Female (sex = 0)
indirect_male <- (a1 + a3) * b          # Male (sex = 1)
IMM <- a3 * b

# =============================================================================
# Results
# =============================================================================

cat("\n=== RESULTS ===\n\n")

slid_results <- tibble(
  Effect = c("a1 (Age -> Education, Female)",
             "a3 (Age x Male -> Education)",
             "b (Education -> Wages)",
             "c' (Direct: Age -> Wages)",
             "Indirect (Female)",
             "Indirect (Male)",
             "IMM"),
  Estimate = round(c(median(a1), median(a3), median(b), median(c_prime),
                     median(indirect_female), median(indirect_male), median(IMM)), 4),
  CI_lower = round(c(quantile(a1, 0.025), quantile(a3, 0.025),
                     quantile(b, 0.025), quantile(c_prime, 0.025),
                     quantile(indirect_female, 0.025), quantile(indirect_male, 0.025),
                     quantile(IMM, 0.025)), 4),
  CI_upper = round(c(quantile(a1, 0.975), quantile(a3, 0.975),
                     quantile(b, 0.975), quantile(c_prime, 0.975),
                     quantile(indirect_female, 0.975), quantile(indirect_male, 0.975),
                     quantile(IMM, 0.975)), 4),
  `P(> 0)` = round(c(mean(a1 > 0), mean(a3 > 0), mean(b > 0), mean(c_prime > 0),
                     mean(indirect_female > 0), mean(indirect_male > 0),
                     mean(IMM > 0)), 4)
)

print(slid_results)

# =============================================================================
# Interpretation
# =============================================================================

IMM_sig <- (quantile(IMM, 0.025) > 0) | (quantile(IMM, 0.975) < 0)

cat("\n=== INTERPRETATION ===\n\n")
cat("IMM =", round(median(IMM), 4),
    "[", round(quantile(IMM, 0.025), 4), ",", round(quantile(IMM, 0.975), 4), "]\n")
cat("95% CI excludes zero:", IMM_sig, "\n\n")

if (IMM_sig) {
  cat("SIGNIFICANT MODERATED MEDIATION DETECTED!\n\n")
  cat("The indirect effect of age on wages through education\n")
  cat("differs significantly between males and females.\n\n")
  cat("Indirect effect for females:", round(median(indirect_female), 4), "\n")
  cat("Indirect effect for males:", round(median(indirect_male), 4), "\n")
}

# =============================================================================
# Convergence Check
# =============================================================================

cat("\n=== Convergence ===\n")
cat("Max Rhat:", round(max(rhat(fit_slid)), 4), "\n")
cat("Min ESS ratio:", round(min(neff_ratio(fit_slid)), 4), "\n")

# Save results
saveRDS(list(fit = fit_slid, results = slid_results, IMM = IMM),
        "slid_results.rds")
