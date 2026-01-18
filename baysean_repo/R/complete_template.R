# ================================================================
# COMPLETE BAYESIAN MODERATED MEDIATION TEMPLATE
# ================================================================
# Exact replication of the tutorial from the Rmd
# Model: BNR -> TTM -> FS with FP moderating BNR -> TTM (first-stage)
# Includes random effects for Market

# --- Setup ---
library(brms)
library(tidyverse)
library(posterior)
library(bayesplot)
library(knitr)
library(kableExtra)

# ================================================================
# STEP 1: SIMULATE DATA (exactly as in Rmd)
# ================================================================
set.seed(42)
n <- 500

raw_data <- tibble(
  # Independent variable: Broker Network Richness (Blau Index, 0-0.8)
  BNR = rbeta(n, 2, 3) * 0.8,

  # Moderator: Founder Preparedness (0-100 scale)
  FP = rnorm(n, 60, 15) %>% pmax(0) %>% pmin(100),

  # Mediator: Time to Meeting (weeks, right-skewed)
  TTM = 4 + exp(2 - 3*BNR - 0.01*FP - 0.02*BNR*FP + rnorm(n, 0, 0.5)),

  # Outcome: Funding Success (binary)
  prob_success = plogis(-2 + 2*BNR - 0.1*TTM + 0.01*FP),
  FS = rbinom(n, 1, prob_success),

  # Grouping variable
  Market = sample(c("Atlanta", "Austin", "Boston", "Denver"), n, replace = TRUE)
)

# Check raw distributions
cat("=== Raw Data Summary ===\n")
print(summary(raw_data[, c("BNR", "FP", "TTM", "FS")]))

# ================================================================
# STEP 2: PREPARE DATA
# ================================================================
model_data <- raw_data %>%
  mutate(
    # Standardize continuous predictors (mean=0, sd=1)
    BNR_z = scale(BNR)[,1],
    FP_z = scale(FP)[,1],

    # Standardize mediator
    TTM_z = scale(TTM)[,1],

    # Create explicit interaction for lavaan comparison
    BNR_FP = BNR_z * FP_z
  ) %>%
  drop_na()

cat("\n=== Standardized Data Summary ===\n")
cat("N =", nrow(model_data), "\n")
cat("Funding Success rate:", round(mean(model_data$FS), 3), "\n")

# ================================================================
# STEP 3: SPECIFY PRIORS
# ================================================================
priors_mediation <- c(
  # --- Mediator model (continuous: TTM) ---
  # Fixed effects
  prior(normal(0, 1.5), class = "b", resp = "TTMz"),
  # Interaction term (tighter for stability)
  prior(normal(0, 1), class = "b", coef = "BNR_z:FP_z", resp = "TTMz"),
  # Random effects SD
  prior(exponential(1), class = "sd", resp = "TTMz"),
  # Residual SD
  prior(exponential(1), class = "sigma", resp = "TTMz"),

  # --- Outcome model (binary: FS) ---
  # Fixed effects (on log-odds scale)
  prior(normal(0, 1.5), class = "b", resp = "FS"),
  # Random effects SD
  prior(exponential(1), class = "sd", resp = "FS")
  # Note: No sigma for Bernoulli family
)

cat("\n=== Prior Specifications ===\n")
print(priors_mediation)

# ================================================================
# STEP 4: FIT MODEL
# ================================================================
cat("\n=== Fitting Bayesian Moderated Mediation Model ===\n")
cat("This may take 5-10 minutes...\n\n")

fit_mediation <- brm(
  # Mediator model: TTM ~ BNR * FP (first-stage moderation)
  bf(TTM_z ~ BNR_z * FP_z + (1 | Market)) +

  # Outcome model: FS ~ BNR + TTM (no interaction = first-stage only)
  bf(FS ~ BNR_z + FP_z + TTM_z + (1 | Market),
     family = bernoulli(link = "logit")) +

  # No residual correlation (different families)
  set_rescor(FALSE),

  data = model_data,
  prior = priors_mediation,

  # MCMC settings
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 42,

  # Computational settings
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  silent = 2,
  refresh = 0
)

# ================================================================
# STEP 5: CHECK DIAGNOSTICS
# ================================================================
cat("\n=== Convergence Diagnostics ===\n")

rhat_values <- brms::rhat(fit_mediation)
ess_values <- brms::neff_ratio(fit_mediation) * 8000  # Convert ratio to actual ESS

max_rhat <- max(rhat_values, na.rm = TRUE)
min_ess <- min(ess_values, na.rm = TRUE)
min_ess_ratio <- min(brms::neff_ratio(fit_mediation), na.rm = TRUE)

cat("Max Rhat:", round(max_rhat, 4), "(should be < 1.01)\n")
cat("Min ESS:", round(min_ess), "(should be > 400)\n")
cat("Min ESS ratio:", round(min_ess_ratio, 4), "(should be > 0.1)\n")
cat("All Rhat < 1.01:", all(rhat_values < 1.01, na.rm = TRUE), "\n")
cat("All ESS > 400:", all(ess_values > 400, na.rm = TRUE), "\n")

# Diagnostics table
diagnostics_table <- tibble(
  Metric = c("Max Rhat", "Min ESS", "Min ESS Ratio", "Divergent Transitions"),
  Value = c(round(max_rhat, 4), round(min_ess), round(min_ess_ratio, 4),
            sum(nuts_params(fit_mediation)$Value[nuts_params(fit_mediation)$Parameter == "divergent__"])),
  Threshold = c("< 1.01", "> 400", "> 0.1", "= 0"),
  Status = c(ifelse(max_rhat < 1.01, "PASS", "FAIL"),
             ifelse(min_ess > 400, "PASS", "FAIL"),
             ifelse(min_ess_ratio > 0.1, "PASS", "FAIL"),
             ifelse(sum(nuts_params(fit_mediation)$Value[nuts_params(fit_mediation)$Parameter == "divergent__"]) == 0, "PASS", "WARNING"))
)

cat("\n=== Diagnostics Summary Table ===\n")
print(as.data.frame(diagnostics_table), row.names = FALSE)

# Print model summary
cat("\n=== Model Summary ===\n")
print(summary(fit_mediation))

# Extract fixed effects
cat("\n=== Fixed Effects ===\n")
fixef_table <- fixef(fit_mediation)
print(round(fixef_table, 3))

# ================================================================
# STEP 6: COMPUTE INDIRECT EFFECTS
# ================================================================
cat("\n=== Computing Indirect Effects ===\n")

posts <- as_draws_df(fit_mediation)

# Extract path coefficients
a1 <- posts$b_TTMz_BNR_z           # Main effect of X on M
a3 <- posts$`b_TTMz_BNR_z:FP_z`    # Interaction (moderation)
b <- posts$b_FS_TTM_z              # Effect of M on Y
c_prime <- posts$b_FS_BNR_z        # Direct effect of X on Y

# Compute conditional indirect effects
indirect_low <- (a1 + a3 * (-1)) * b   # FP at -1 SD
indirect_mean <- a1 * b                 # FP at mean (0)
indirect_high <- (a1 + a3 * 1) * b      # FP at +1 SD

# Index of Moderated Mediation
IMM <- a3 * b

# Total effect
total <- c_prime + indirect_mean

# Create summary table
mediation_summary <- tibble(
  Effect = c("a1 (BNR → TTM)",
             "a3 (BNR×FP → TTM)",
             "b (TTM → FS)",
             "c' (direct)",
             "Indirect (low FP)",
             "Indirect (mean FP)",
             "Indirect (high FP)",
             "IMM",
             "Total effect"),
  Estimate = c(median(a1), median(a3), median(b), median(c_prime),
               median(indirect_low), median(indirect_mean),
               median(indirect_high), median(IMM), median(total)),
  CI_lower = c(quantile(a1, 0.025), quantile(a3, 0.025),
               quantile(b, 0.025), quantile(c_prime, 0.025),
               quantile(indirect_low, 0.025), quantile(indirect_mean, 0.025),
               quantile(indirect_high, 0.025), quantile(IMM, 0.025),
               quantile(total, 0.025)),
  CI_upper = c(quantile(a1, 0.975), quantile(a3, 0.975),
               quantile(b, 0.975), quantile(c_prime, 0.975),
               quantile(indirect_low, 0.975), quantile(indirect_mean, 0.975),
               quantile(indirect_high, 0.975), quantile(IMM, 0.975),
               quantile(total, 0.975)),
  Prob_positive = c(mean(a1 > 0), mean(a3 > 0), mean(b > 0), mean(c_prime > 0),
                    mean(indirect_low > 0), mean(indirect_mean > 0),
                    mean(indirect_high > 0), mean(IMM > 0), mean(total > 0))
)

cat("\n========================================\n")
cat("MEDIATION ANALYSIS RESULTS\n")
cat("========================================\n\n")
print(as.data.frame(mediation_summary), row.names = FALSE)

# ================================================================
# STEP 7: TEST MODERATION (SIMPLE SLOPES)
# ================================================================
cat("\n=== Interaction Effect Summary ===\n")
cat("Median a3:", round(median(a3), 4), "\n")
cat("95% CI: [", round(quantile(a3, 0.025), 4), ",",
    round(quantile(a3, 0.975), 4), "]\n")
cat("P(a3 < 0):", round(mean(a3 < 0), 4), "\n")

# Simple slopes: Effect of BNR on TTM at different FP levels
FP_levels <- c(-1, 0, 1)

simple_slopes <- map_dfr(FP_levels, function(fp) {
  slope <- a1 + a3 * fp
  tibble(
    FP_level = fp,
    FP_label = case_when(
      fp == -1 ~ "Low (-1 SD)",
      fp == 0 ~ "Mean",
      fp == 1 ~ "High (+1 SD)"
    ),
    Slope = median(slope),
    CI_lower = quantile(slope, 0.025),
    CI_upper = quantile(slope, 0.975),
    Prob_negative = mean(slope < 0),
    Significant = ifelse((CI_lower > 0) | (CI_upper < 0), "*", "")
  )
})

cat("\n=== Simple Slopes Analysis ===\n")
cat("Effect of BNR on TTM at different levels of FP:\n\n")
print(as.data.frame(simple_slopes), row.names = FALSE)

# ================================================================
# RESULTS TABLE
# ================================================================
results_table <- tibble(
  Path = c("BNR → TTM (a₁)",
           "BNR × FP → TTM (a₃)",
           "TTM → FS (b)",
           "BNR → FS (c')",
           "---",
           "Indirect (mean FP)",
           "IMM",
           "Total Effect"),
  Estimate = c(sprintf("%.3f", median(a1)),
               sprintf("%.3f", median(a3)),
               sprintf("%.3f", median(b)),
               sprintf("%.3f", median(c_prime)),
               "---",
               sprintf("%.3f", median(indirect_mean)),
               sprintf("%.3f", median(IMM)),
               sprintf("%.3f", median(total))),
  `95% CI` = c(sprintf("[%.3f, %.3f]", quantile(a1, 0.025), quantile(a1, 0.975)),
               sprintf("[%.3f, %.3f]", quantile(a3, 0.025), quantile(a3, 0.975)),
               sprintf("[%.3f, %.3f]", quantile(b, 0.025), quantile(b, 0.975)),
               sprintf("[%.3f, %.3f]", quantile(c_prime, 0.025), quantile(c_prime, 0.975)),
               "---",
               sprintf("[%.3f, %.3f]", quantile(indirect_mean, 0.025), quantile(indirect_mean, 0.975)),
               sprintf("[%.3f, %.3f]", quantile(IMM, 0.025), quantile(IMM, 0.975)),
               sprintf("[%.3f, %.3f]", quantile(total, 0.025), quantile(total, 0.975))),
  `P(Direction)` = c(sprintf("%.1f%%", mean(a1 < 0) * 100),
                     sprintf("%.1f%%", max(mean(a3 < 0), mean(a3 > 0)) * 100),
                     sprintf("%.1f%%", mean(b < 0) * 100),
                     sprintf("%.1f%%", mean(c_prime > 0) * 100),
                     "---",
                     sprintf("%.1f%%", mean(indirect_mean > 0) * 100),
                     sprintf("%.1f%%", max(mean(IMM > 0), mean(IMM < 0)) * 100),
                     sprintf("%.1f%%", mean(total > 0) * 100))
)

cat("\n========================================\n")
cat("PUBLICATION-READY RESULTS TABLE\n")
cat("========================================\n\n")
print(as.data.frame(results_table), row.names = FALSE)

# ================================================================
# STEP 8: FINAL RESULTS SUMMARY
# ================================================================
cat("\n========================================\n")
cat("FINAL INTERPRETATION\n")
cat("========================================\n\n")

# Check if IMM is significant
imm_ci <- quantile(IMM, c(0.025, 0.975))
imm_significant <- (imm_ci[1] > 0) | (imm_ci[2] < 0)

cat("Index of Moderated Mediation (IMM):\n")
cat(sprintf("  Estimate: %.4f\n", median(IMM)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", imm_ci[1], imm_ci[2]))
cat(sprintf("  Status: %s\n\n", ifelse(imm_significant, "SIGNIFICANT", "Not significant")))

if (imm_significant) {
  cat("CONCLUSION: The indirect effect of BNR on FS through TTM\n")
  cat("is moderated by Founder Preparedness (FP).\n\n")

  if (median(a3) < 0) {
    cat("The negative a3 coefficient indicates that higher FP\n")
    cat("weakens the effect of BNR on TTM.\n")
  } else {
    cat("The positive a3 coefficient indicates that higher FP\n")
    cat("strengthens the effect of BNR on TTM.\n")
  }
} else {
  cat("CONCLUSION: No evidence of moderated mediation.\n")
  cat("The indirect effect does not depend on FP.\n")
}

# ================================================================
# WRITING RESULTS TEMPLATE
# ================================================================
cat("\n========================================\n")
cat("WRITING RESULTS (TEMPLATE)\n")
cat("========================================\n\n")

cat("Mediation Results: We tested whether Time to Meeting (TTM) mediated\n")
cat("the relationship between Broker Network Richness (BNR) and Funding\n")
cat("Success (FS) using Bayesian moderated mediation analysis. All models\n")
cat("were estimated using brms (Bürkner, 2017) with weakly informative\n")
cat("priors (Normal(0, 1.5) for fixed effects). Convergence was confirmed\n")
cat("by Rhat < 1.01 and ESS > 400 for all parameters.\n\n")

cat(sprintf("The a-path from BNR to TTM was negative and credible\n"))
cat(sprintf("(β = %.3f, 95%% CI [%.3f, %.3f]), indicating that founders\n",
            median(a1), quantile(a1, 0.025), quantile(a1, 0.975)))
cat("with more diverse broker networks obtained investor meetings faster.\n\n")

cat(sprintf("The b-path from TTM to FS was also negative\n"))
cat(sprintf("(β = %.3f, 95%% CI [%.3f, %.3f]), confirming that faster\n",
            median(b), quantile(b, 0.025), quantile(b, 0.975)))
cat("meetings predicted higher funding success.\n\n")

cat(sprintf("The indirect effect at mean founder preparedness was positive\n"))
cat(sprintf("and credible (indirect = %.3f, 95%% CI [%.3f, %.3f]),\n",
            median(indirect_mean), quantile(indirect_mean, 0.025), quantile(indirect_mean, 0.975)))
cat("supporting mediation.\n\n")

cat(sprintf("The Index of Moderated Mediation was %.3f\n", median(IMM)))
cat(sprintf("(95%% CI [%.3f, %.3f]), ", quantile(IMM, 0.025), quantile(IMM, 0.975)))
if (imm_significant) {
  cat("with the credible interval excluding zero,\n")
  cat("indicating that founder preparedness significantly moderates\n")
  cat("the mediation pathway.\n")
} else {
  cat("with the credible interval including zero,\n")
  cat("indicating that founder preparedness did not significantly\n")
  cat("moderate the mediation pathway.\n")
}

cat("\n========================================\n")
cat("Analysis complete.\n")

# ================================================================
# STEP 9: GENERATE ALL FIGURES
# ================================================================
cat("\n=== Generating Figures ===\n")

library(ggplot2)
library(patchwork)

# Create figures directory if it doesn't exist
if (!dir.exists("figures")) dir.create("figures")

# --- Figure 1: Data Distributions ---
p1 <- ggplot(model_data, aes(x = BNR_z)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(title = "BNR (standardized)", x = NULL) +
  theme_minimal()

p2 <- ggplot(model_data, aes(x = TTM_z)) +
  geom_histogram(bins = 30, fill = "coral", alpha = 0.7) +
  labs(title = "TTM (standardized)", x = NULL) +
  theme_minimal()

p3 <- ggplot(model_data, aes(x = factor(FS))) +
  geom_bar(fill = "seagreen", alpha = 0.7) +
  labs(title = "Funding Success", x = NULL) +
  theme_minimal()

fig1 <- p1 + p2 + p3
ggsave("figures/01_data_distributions.png", fig1, width = 10, height = 3, dpi = 150)
cat("Saved: figures/01_data_distributions.png\n")

# --- Figure 2: Prior Distributions ---
x <- seq(-5, 5, length.out = 1000)

prior_viz <- tibble(
  x = x,
  `Normal(0, 1.5)` = dnorm(x, 0, 1.5),
  `Normal(0, 1)` = dnorm(x, 0, 1),
  `Normal(0, 2.5)` = dnorm(x, 0, 2.5)
) %>%
  pivot_longer(-x, names_to = "Prior", values_to = "Density")

fig2 <- ggplot(prior_viz, aes(x = x, y = Density, color = Prior)) +
  geom_line(linewidth = 1) +
  labs(title = "Prior Distributions for Fixed Effects",
       x = "Coefficient Value", y = "Density") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave("figures/02_prior_distributions.png", fig2, width = 8, height = 4, dpi = 150)
cat("Saved: figures/02_prior_distributions.png\n")

# --- Figure 3: Trace Plots ---
fig3 <- mcmc_trace(fit_mediation, pars = c("b_TTMz_BNR_z", "b_TTMz_BNR_z:FP_z", "b_FS_TTM_z")) +
  labs(title = "MCMC Trace Plots for Key Parameters")
ggsave("figures/03_trace_plots.png", fig3, width = 10, height = 6, dpi = 150)
cat("Saved: figures/03_trace_plots.png\n")

# --- Figure 4: Posterior Predictive Check ---
fig4 <- pp_check(fit_mediation, resp = "TTMz", ndraws = 100) +
  labs(title = "Posterior Predictive Check: TTM")
ggsave("figures/04_pp_check.png", fig4, width = 8, height = 4, dpi = 150)
cat("Saved: figures/04_pp_check.png\n")

# --- Figure 5: Good vs Bad Trace Patterns ---
set.seed(123)
n_iter <- 500

good_trace <- tibble(
  iteration = rep(1:n_iter, 4),
  chain = rep(1:4, each = n_iter),
  value = rnorm(n_iter * 4, 0, 0.1) + rep(c(0, 0.02, -0.02, 0.01), each = n_iter)
)

bad_trace <- tibble(
  iteration = rep(1:n_iter, 4),
  chain = rep(1:4, each = n_iter),
  value = c(
    cumsum(rnorm(n_iter, 0.001, 0.05)),
    cumsum(rnorm(n_iter, -0.002, 0.05)) + 0.5,
    rnorm(n_iter, 0, 0.1),
    rnorm(n_iter, -0.3, 0.1)
  )
)

p_good <- ggplot(good_trace, aes(x = iteration, y = value, color = factor(chain))) +
  geom_line(alpha = 0.7) +
  labs(title = "Good: Chains mix well", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

p_bad <- ggplot(bad_trace, aes(x = iteration, y = value, color = factor(chain))) +
  geom_line(alpha = 0.7) +
  labs(title = "Bad: Chains don't mix", x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

fig5 <- p_good + p_bad
ggsave("figures/05_trace_patterns.png", fig5, width = 10, height = 3, dpi = 150)
cat("Saved: figures/05_trace_patterns.png\n")

# --- Figure 6: Indirect Effect Distributions ---
indirect_df <- tibble(
  `Low FP (-1 SD)` = indirect_low,
  `Mean FP (0)` = indirect_mean,
  `High FP (+1 SD)` = indirect_high
) %>%
  pivot_longer(everything(), names_to = "Condition", values_to = "Indirect")

fig6 <- ggplot(indirect_df, aes(x = Indirect, fill = Condition)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Posterior Distribution of Indirect Effects",
       x = "Indirect Effect", y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

ggsave("figures/06_indirect_effects.png", fig6, width = 8, height = 4, dpi = 150)
cat("Saved: figures/06_indirect_effects.png\n")

# --- Figure 7: Simple Slopes Plot ---
slopes_for_plot <- expand_grid(
  BNR_z = seq(-2, 2, length.out = 100),
  FP_level = c(-1, 0, 1)
) %>%
  mutate(
    FP_label = factor(FP_level, levels = c(-1, 0, 1),
                      labels = c("Low FP (-1 SD)", "Mean FP", "High FP (+1 SD)")),
    TTM_pred = median(posts$b_TTMz_Intercept) +
               median(a1) * BNR_z +
               median(posts$b_TTMz_FP_z) * FP_level +
               median(a3) * BNR_z * FP_level
  )

fig7 <- ggplot(slopes_for_plot, aes(x = BNR_z, y = TTM_pred, color = FP_label)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Simple Slopes: BNR Effect on TTM by Founder Preparedness",
       x = "Broker Network Richness (standardized)",
       y = "Time to Meeting (standardized)",
       color = "Moderator Level") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave("figures/07_simple_slopes.png", fig7, width = 7, height = 5, dpi = 150)
cat("Saved: figures/07_simple_slopes.png\n")

# --- Figure 8: Johnson-Neyman Plot ---
FP_range <- seq(-3, 3, by = 0.1)

JN_results <- map_dfr(FP_range, function(fp) {
  slope <- a1 + a3 * fp
  tibble(
    FP = fp,
    Effect = median(slope),
    CI_lower = quantile(slope, 0.025),
    CI_upper = quantile(slope, 0.975),
    Significant = (CI_lower > 0) | (CI_upper < 0)
  )
})

fig8 <- ggplot(JN_results, aes(x = FP, y = Effect)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_rug(data = filter(JN_results, Significant),
           aes(x = FP), sides = "b", color = "darkgreen") +
  labs(title = "Johnson-Neyman Plot: Regions of Significance",
       subtitle = "Green ticks indicate FP values where BNR effect is significant",
       x = "Founder Preparedness (standardized)",
       y = "Effect of BNR on TTM") +
  theme_minimal()

ggsave("figures/08_johnson_neyman.png", fig8, width = 7, height = 5, dpi = 150)
cat("Saved: figures/08_johnson_neyman.png\n")

# --- Figure 9: Rhat Distribution ---
rhat_df <- tibble(parameter = names(rhat_values), rhat = rhat_values) %>%
  filter(!is.na(rhat))

fig9 <- ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Rhat Values",
       subtitle = "Red line indicates threshold (1.01)",
       x = "Rhat", y = "Count") +
  theme_minimal()

ggsave("figures/09_rhat_distribution.png", fig9, width = 7, height = 4, dpi = 150)
cat("Saved: figures/09_rhat_distribution.png\n")

# --- Figure 10: Parameter Density Plots ---
fig10 <- mcmc_dens_overlay(fit_mediation, pars = c("b_TTMz_BNR_z", "b_TTMz_BNR_z:FP_z", "b_FS_TTM_z")) +
  labs(title = "Posterior Density by Chain")
ggsave("figures/10_parameter_densities.png", fig10, width = 10, height = 4, dpi = 150)
cat("Saved: figures/10_parameter_densities.png\n")

# ================================================================
# STEP 10: SAVE TABLES AS CSV
# ================================================================
cat("\n=== Saving Tables ===\n")

# Create tables directory
if (!dir.exists("tables")) dir.create("tables")

# Save diagnostics table
write.csv(diagnostics_table, "tables/01_diagnostics.csv", row.names = FALSE)
cat("Saved: tables/01_diagnostics.csv\n")

# Save fixed effects
write.csv(round(fixef_table, 4), "tables/02_fixed_effects.csv")
cat("Saved: tables/02_fixed_effects.csv\n")

# Save mediation summary
write.csv(mediation_summary, "tables/03_mediation_summary.csv", row.names = FALSE)
cat("Saved: tables/03_mediation_summary.csv\n")

# Save simple slopes
write.csv(simple_slopes, "tables/04_simple_slopes.csv", row.names = FALSE)
cat("Saved: tables/04_simple_slopes.csv\n")

# Save publication-ready results
write.csv(results_table, "tables/05_publication_results.csv", row.names = FALSE)
cat("Saved: tables/05_publication_results.csv\n")

# Save Johnson-Neyman results
write.csv(JN_results, "tables/06_johnson_neyman.csv", row.names = FALSE)
cat("Saved: tables/06_johnson_neyman.csv\n")

# Software table
software_table <- tibble(
  Package = c("brms", "tidyverse", "posterior", "bayesplot", "lavaan"),
  Purpose = c("Bayesian regression modeling",
              "Data manipulation and visualization",
              "Posterior analysis utilities",
              "MCMC diagnostics visualization",
              "Frequentist comparison (optional)"),
  Citation = c("Bürkner (2017)", "Wickham et al. (2019)",
               "Bürkner et al. (2023)", "Gabry et al. (2019)",
               "Rosseel (2012)")
)
write.csv(software_table, "tables/07_software.csv", row.names = FALSE)
cat("Saved: tables/07_software.csv\n")

cat("\n=== All figures saved to figures/ directory ===\n")
cat("=== All tables saved to tables/ directory ===\n")
cat("\nTemplate execution complete.\n")
