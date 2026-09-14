# Check:
# Error structure consistency, error_it = nu_i + eta_it --> sometimes inside mutate before full panel alignment
# ATT(g,t) is no longer well-defined without conditioning on entire history. Estimand becomes: ATT(g,t,history)
# Simulation 5 models "dynamic treatment response function"
# ToDo:
# Tweak violation of parallel trends more


# ==============================================================================
# Master Data Simulation: Complex DiD Environments
# ==============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(did)
library(ggplot2)
library(gridExtra)

set.seed(20260608) # Grounded seed

# --- Global Parameter Parameters ---
N <- 1000
T <- 10

# 1. Generate Static Unit-Level Characteristics
units <- tibble(
  id = 1:N,
  alpha_i = rnorm(N, mean = 0, sd = 1),     # Unit fixed effect
  nu_i    = rnorm(N, mean = 0, sd = 0.5),   # Unit-level error cluster component
  X1      = rnorm(N, mean = 2, sd = 1),     # Confounder driving conditional trends
  X2      = rnorm(N, mean = 0, sd = 1),     # Noise covariate 1
  X3      = rbinom(N, 1, 0.5)               # Noise covariate 2
) %>%
  mutate(
    # Base cohort assignment for staggered settings (50% end up treated)
    cohort_prob = exp(X1) / (1 + exp(X1)), # Confounding entry likelihood
    G = sample(c(2:10, 0), size = N, replace = TRUE, 
               prob = c(rep(0.5/9, 9), 0.5))
  )

# 2. Build Structural Panel Shell
panel_base <- expand_grid(id = 1:N, period = 1:T) %>%
  left_join(units, by = "id") %>%
  mutate(
    theta_t  = sin(period / 2),              # Time fixed effect
    eta_it   = rnorm(n(), mean = 0, sd = 0.5),# Idiosyncratic shock
    error_it = nu_i + eta_it                 # Clustered error term
  )

# ==============================================================================
# SETTING 1: Staggered Binary, Heterogeneous Effects, Unconditional Parallel Trends
# ==============================================================================
df_setting1 <- panel_base %>%
  mutate(
    treatment_status = if_else(G != 0 & period >= G, 1, 0),
    # Causal Effect: Dynamically grows over time and varies by cohort group
    true_att = if_else(treatment_status == 1, 1.5 * (period - G + 1) + (G * 0.1), 0),
    # Potential Outcomes
    Y0 = alpha_i + theta_t + 0.3*X1 + 0.2*X2 + error_it,
    Y1 = Y0 + true_att,
    observed_outcome = if_else(treatment_status == 1, Y1, Y0),
    treatment_dose = NA_real_,
    cumulative_dose = NA_real_
  )

# ==============================================================================
# SETTING 2: Staggered Binary, Heterogeneous Effects, Conditional Parallel Trends
# ==============================================================================
df_setting2 <- panel_base %>%
  mutate(
    treatment_status = if_else(G != 0 & period >= G, 1, 0),
    true_att = if_else(treatment_status == 1, 1.5 * (period - G + 1) + (G * 0.1), 0),
    # Potential Outcome Y0 contains a non-linear covariate-by-time interaction
    Y0 = alpha_i + theta_t + 0.3*X1 + 0.5*(X1 * period) + 0.2*X2 + error_it,
    Y1 = Y0 + true_att,
    observed_outcome = if_else(treatment_status == 1, Y1, Y0),
    treatment_dose = NA_real_,
    cumulative_dose = NA_real_
  )

# ==============================================================================
# SETTING 3: Continuous Doses & Conditional Parallel Trends (Callaway et al., 2024)
# ==============================================================================
df_setting3 <- panel_base %>%
  mutate(
    treatment_status = if_else(G != 0 & period >= G, 1, 0),
    # Assign permanent intensity dosage if treated
    treatment_dose = if_else(G != Inf, sample(c(1, 2, 3), n(), replace = TRUE), 0),
    treatment_dose = if_else(treatment_status == 1, treatment_dose, 0),
    # ATT is scaled non-linearly by treatment dose intensity
    true_att = if_else(treatment_status == 1, (1.2 * (period - G + 1)) * (treatment_dose^1.2), 0),
    Y0 = alpha_i + theta_t + 0.3*X1 + 0.5*(X1 * period) + 0.2*X2 + error_it,
    Y1 = Y0 + true_att,
    observed_outcome = if_else(treatment_status == 1, Y1, Y0)
  ) %>%
  group_by(id) %>%
  mutate(cumulative_dose = cumsum(treatment_dose)) %>%
  ungroup()

# ==============================================================================
# SETTING 4: Reversible Treatments (On/Off) with Conditional Parallel Trends
# ==============================================================================
df_setting4 <- panel_base %>%
  mutate(
    # Generate an arbitrary alternating on-off mechanism tied to unit attributes
    treatment_status = if_else(G != 0 & period >= G & (period %% 3 != 0), 1, 0),
    # Shock effect manifests during active periods
    true_att = if_else(treatment_status == 1, 2.0 + 0.5 * (period - G), 0),
    Y0 = alpha_i + theta_t + 0.3*X1 + 0.5*(X1 * period) + 0.2*X2 + error_it,
    Y1 = Y0 + true_att,
    observed_outcome = if_else(treatment_status == 1, Y1, Y0),
    treatment_dose = NA_real_,
    cumulative_dose = NA_real_
  )

# ==============================================================================
# SETTING 5: Reversible Multi-Dose Treatment Matrix (Fully Complex)
# ==============================================================================
df_setting5 <- panel_base %>%
  mutate(
    # Multi-dose assignment can hit zero or shift intensity fluidly at any time
    treatment_dose = case_when(
      G == 0 ~ 0,
      period < G ~ 0,
      (period + id) %% 4 == 0 ~ 1,
      (period + id) %% 4 == 1 ~ 2,
      (period + id) %% 4 == 2 ~ 3,
      TRUE ~ 0
    ),
    treatment_status = if_else(treatment_dose > 0, 1, 0)
  ) %>%
  group_by(id) %>%
  mutate(cumulative_dose = cumsum(treatment_dose)) %>%
  ungroup() %>%
  mutate(
    att_gt  = 1.5 * (period - G + 1),
    att_d   = 0.8 * treatment_dose,
    att_cum = 0.3 * cumulative_dose,
    # ATT depends on both immediate intensity injection and historical usage
    true_att = if_else(treatment_status == 1, 
                       att_gt + att_d + att_cum, 0),
    Y0 = alpha_i + theta_t + 0.3*X1 + 0.5*(X1 * period) + 0.2*X2 + error_it,
    Y1 = Y0 + true_att,
    observed_outcome = if_else(treatment_status == 1, Y1, Y0)
  )

# ==============================================================================
# Output and Verification Check
# ==============================================================================
# Sanity check validation of dimensions
cat("Dimensions Check (Rows x Columns):\n")
cat("Setting 1:", dim(df_setting1), "\n")
cat("Setting 2:", dim(df_setting2), "\n")
cat("Setting 3:", dim(df_setting3), "\n")
cat("Setting 4:", dim(df_setting4), "\n")
cat("Setting 5:", dim(df_setting5), "\n")

# View structure snippet of the final complex model
head(df_setting5 %>% select(id, period, X1, treatment_status, treatment_dose, cumulative_dose, Y0, Y1, observed_outcome))





# ==============================================================================
# Master Data Simulation: Complex DiD Environments
# ==============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(did)        # For att_gt and aggte
library(ggplot2)    # For plotting
library(gridExtra)  # For grid.arrange


# ==============================================================================
# CS ESTIMATOR ADAPTATION
# ==============================================================================

# 1. Define your variables based on the simulated dataset
outcome_vars <- c("observed_outcome")

# Optional: If you want conditional parallel trends, use X1 and X2.
# For unconditional parallel trends as specified in Setting 1, leave this NULL or pass NULL.
covariates_formula <- NULL  # Alternative: ~ X1 + X2

# 3. Create a function to run DiD and Event Study for a single outcome
run_cs_did <- function(outcome, data, covariates_formula = NULL) {
  message(paste("Running CS-DiD for:", outcome))
  
  # Estimate group-time average treatment effects (ATT(g,t))
  atts <- att_gt(
    yname = outcome,
    tname = "period",                # Simulated time variable
    idname = "id",                    # Simulated unit ID
    gname = "G",                     # Simulated treatment cohort variable
    xformla = covariates_formula,
    data = data,
    panel = TRUE,
    allow_unbalanced_panel = FALSE,  # Data is a balanced panel shell
    clustervars = "id",              # Clustering at the unit level
    control_group = "notyettreated", # Alternatively "nevertreated"
    anticipation = 0,
    bstrap      = TRUE,
    biters      = 1000               # Reduced from 100k for faster execution
  )
  
  # Aggregate into Event Study (Dynamic Effects)
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  
  return(list(atts = atts, es = es))
}

# 4. Loop through all outcomes using purrr::map
did_results <- map(outcome_vars, ~run_cs_did(.x, df_setting5, covariates_formula = covariates_formula))
names(did_results) <- outcome_vars

# 5. Results

# Summary for Observed Outcome
summary(did_results$observed_outcome$es)

# Plot the Event Study for Observed Outcome
ggdid(did_results$observed_outcome$es) + 
  ggtitle("Event Study: Effect of Treatment on Observed Outcome") +
  theme_minimal()


# Extract ATT(g,t) estimates for Observed Outcome
att_observed <- did_results$observed_outcome$atts

# Build a tidy dataframe from the att_gt object
att_df <- data.frame(
  group    = att_observed$group,   # treatment cohort (G value)
  time     = att_observed$t,       # calendar time (period value)
  att      = att_observed$att,
  se       = att_observed$se
) %>%
  mutate(
    event_time = time - group,    # relative time to treatment
    ci_low  = att - 1.96 * se,
    ci_high = att + 1.96 * se,
    cohort  = factor(paste("Cohort", group))
  ) %>%
  filter(group > 0) # Remove never-treated group placeholder if present in rows

# Plot: one panel per cohort, x-axis = event time
ggplot(att_df, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(
    title    = "ATT(g,t) by Treatment Cohort — Observed Outcome",
    subtitle = "Each panel shows one treatment cohort; red line = treatment onset",
    x        = "Event Time (periods relative to treatment)",
    y        = "ATT"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
