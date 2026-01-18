# =============================================================================
# 01_data_preparation.R
# Bayesian Moderated Mediation Tutorial
# Step 2: Prepare Data
# =============================================================================

library(tidyverse)

# --- Simulate example data ---
set.seed(42)
n <- 500

# Generate correlated predictors
BNR <- rnorm(n, 0, 1)  # Broker Network Richness
FP <- rnorm(n, 0, 1)   # Founder Preparedness

# Generate mediator with interaction
TTM <- 0.4 * BNR + 0.3 * FP + 0.25 * BNR * FP + rnorm(n, 0, 0.8)

# Generate outcome
FS_latent <- -0.3 * TTM + 0.2 * BNR + rnorm(n, 0, 1)
FS <- rbinom(n, 1, plogis(FS_latent))

# Create data frame
model_data <- tibble(
  BNR = BNR,
  FP = FP,
  TTM = TTM,
  FS = FS
)

# --- Check raw distributions ---
model_data %>%
  pivot_longer(cols = c(BNR, FP, TTM), names_to = "variable") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(title = "Raw Variable Distributions")

# --- Standardize predictors ---
# IMPORTANT: Standardize BEFORE creating interaction terms
model_data <- model_data %>%
  mutate(
    BNR_z = scale(BNR)[,1],
    FP_z = scale(FP)[,1],
    TTM_z = scale(TTM)[,1],
    # Create interaction AFTER standardizing
    BNR_FP = BNR_z * FP_z
  )

# --- Check distributions ---
summary(model_data)

# Verify standardization
cat("Mean BNR_z:", mean(model_data$BNR_z), "\n")
cat("SD BNR_z:", sd(model_data$BNR_z), "\n")

# =============================================================================
# COMMON PITFALLS
# =============================================================================

# WRONG: Don't standardize binary outcomes
# FS_z <- scale(model_data$FS)  # NO!

# CORRECT: Keep binary outcomes as 0/1
# Binary outcomes stay as-is; brms handles them with family = bernoulli()

# -----------------------------------------------------------------------------

# WRONG: Creates correlation between main effect and interaction
# BAD_interaction <- model_data$BNR * model_data$FP
# BAD_BNR_z <- scale(model_data$BNR)
# Correlation between BAD_BNR_z and BAD_interaction can exceed 0.5!

# CORRECT: Standardize first, then create interaction
model_data <- model_data %>%
  mutate(
    BNR_z = scale(BNR)[,1],
    FP_z = scale(FP)[,1],
    BNR_FP = BNR_z * FP_z  # Interaction of standardized variables
  )

# Check correlation (should be near 0)
cat("Correlation between BNR_z and BNR_FP:",
    cor(model_data$BNR_z, model_data$BNR_FP), "\n")

# =============================================================================
# Save prepared data
# =============================================================================

saveRDS(model_data, "model_data_prepared.rds")
cat("Data preparation complete. Saved to model_data_prepared.rds\n")
