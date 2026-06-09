###################### Descriptive Analysis of Panel Data ######################

###################### Load packages & data ######################

library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)

library(fixest)
library(bacondecomp)
library(did)
library(HonestDiD)
library(contdid)
library(fect)

###################### Set the scene ######################

#### Load data ####
df_panel_cov <- readRDS('_data/df_panel_cov.rds')


#### Re-construct treatment variables for different packages & add additional subgroup dummies ####

# Create sorted vector of election years
unique_years <- sort(unique(df_panel_cov$election_year))

# Add variables
df_did_ready <- df_panel_cov %>%
  mutate(
    # Create sequential time index (1 to 10) for did package
    seq_time = match(election_year, unique_years),
    # Create sequential group index (never treated == 0) for did package
    seq_group = match(treat_absorbing_cs, unique_years),
    seq_group = ifelse(is.na(seq_group), 0, seq_group),
    # Create time-to-treatment variable for fixest package
    time_to_treatment = ifelse(seq_group > 0, seq_time - seq_group, -1000),
    # Clean lagged wind count variable of NAs
    cat_cum_lag_wind_count_3km = ifelse(election_year == 1990 & is.na(cat_cum_lag_wind_count_3km),
                                        0, cat_cum_lag_wind_count_3km),
    # Create current incumbent election results
    current_incumbent = case_when(
      election_year %in% c(1990, 1994, 1998, 2013) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(fdp, 0),
      election_year %in% c(2002, 2005) ~ coalesce(spd, 0) + coalesce(gruene, 0),
      election_year %in% c(2009, 2017, 2021) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(spd, 0),
      election_year == 2025 ~ coalesce(gruene, 0) + coalesce(spd, 0) + coalesce(fdp, 0),
      TRUE ~ NA_real_),
    # Create "other parties" election results
    others = 1 - (coalesce(cdu, 0) + coalesce(afd, 0) + coalesce(csu, 0)
                  + coalesce(spd, 0) + coalesce(fdp, 0) + coalesce(linke_pds, 0)
                  + coalesce(gruene, 0)),
    # Add east Germany dummy
    east_ger = ifelse(substr(ags, 1, 2) %in% c("12", "13", "14", "15", "16"), 1, 0)) %>%
  # Filter for certain states
  # filter(str_starts(ags, "12")) %>%
  # Filter for certain time period
  # filter(!election_year < 1990)
  # Ensure ID variable is numeric for the did package
  mutate(ags = as.numeric(ags)) %>%
  as.data.frame()



###################### Outcomes ######################

outcome_vars <- c("turnout", "cdu", "csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")


###################### TWFE Event Study ######################

# Controls
covariates_str <- "+ pop_density + cum_lag_wind_count_3km"


# Function for TWFE event study estimation
run_twfe_event_study <- function(outcome, data, fixed_effects = "ags + election_year", covariates = "") {
  # Formula construction: outcome ~ treatment + covariates | fixed effects
  formula_str <- paste0(outcome, " ~ i(time_to_treatment, ref = c(-1, -1000)) ", 
    covariates, " | ", fixed_effects)
  # Run the TWFE model
  model <- feols(as.formula(formula_str), data = data, cluster = ~ags)
  return(list(model = model))
}


# Estimate TWFE event study
twfe_results <- map(outcome_vars, ~run_twfe_event_study(
  outcome = .x, 
  data = df_did_ready, 
  covariates = covariates_str
))
names(twfe_results) <- outcome_vars


# Summary of outcomes
etable(map(twfe_results, ~ .x$model), keep = "time_to_treatment")


# Plot results
par(mfrow = c(3, 3))
for (i in 1:9) {
  iplot(twfe_results[[i]], 
        main = paste(outcome_vars[i], "(Ref: 1 Year Pre-Treatment)"))
}
par(mfrow = c(1, 1))



###################### Goodman-Bacon Decomposition (Goodman-Bacon, 2021) ######################
# Running decomposition on all outcomes takes some time. Possibly parallelise code here with future.

# Covariates
covariate_str <- "~ treat_absorbing + pop_density + cum_lag_wind_count_3km"

did_short <- df_did_ready %>%
  filter(election_year > 2013, str_starts(ags, "12"))

run_bacon_decomposition <- function(outcome, data, covariates) {
  # Drop NAs
  df_clean <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(treat_absorbing), !is.na(ags), !is.na(seq_time), !is.na(pop_density)
    )
  # Find no. of time periods
  expected_periods <- n_distinct(df_clean$seq_time)
  # Balance panel
  df_balanced <- df_clean %>%
    group_by(ags) %>%
    filter(n() == expected_periods) %>%  
    ungroup()
  # Check if any data left after balancing
  if (nrow(df_balanced) == 0) {
    warning(paste("No balanced data left for outcome:", outcome))
    return(NULL)
  }
  # Run decomposition
  fml <- as.formula(paste(outcome, covariates))
  bgd <- bacon(fml, 
               data = df_balanced, 
               id_var = "ags", 
               time_var = "seq_time")
  return(bgd)
}

# Run decomposition
bacon_results <- map(outcome_vars, ~run_bacon_decomposition(.x, did_short, covariate_str))
names(bacon_results) <- outcome_vars

# Combine 2x2 decompositions into df
bacon_plot_data <- map_dfr(names(bacon_results), function(outcome) {
  df_2x2 <- bacon_results[[outcome]]$two_by_twos
  if (!is.null(df_2x2)) {
    # Add column to identify outcome
    df_2x2$outcome_var <- outcome
    return(df_2x2)
  }
  return(NULL)
})

# Calculate overall TWFE (weighted mean) for each outcome group
twfe_means <- bacon_plot_data %>%
  group_by(outcome_var) %>%
  summarize(bgd_wm = weighted.mean(estimate, weight), .groups = "drop")

# Plot results
ggplot(bacon_plot_data, aes(x = weight, y = estimate, shape = type, col = type)) +
  # Add outcome-specific TWFE horizontal line
  geom_hline(data = twfe_means, aes(yintercept = bgd_wm), lty = 2, color = "darkgray") +
  # Add  points
  geom_point(size = 2.5, alpha = 0.8) +
  # Create grid
  facet_wrap(~ outcome_var, scales = "free", ncol = 3) + 
  # Styling
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 11), # Style grid titles
    legend.position = "bottom"
  ) +
  labs(
    x = "Weight", 
    y = "Estimate", 
    shape = "Type", 
    col = "Type",
    title = "Bacon-Goodman Decomposition by Outcome",
    subtitle = "Dotted lines depict the full TWFE estimate for each outcome."
  )



###################### Staggered DiD with binary, absorbing treatment (Callaway & Sant'Anna, 2021) ######################
# Only conditional parallel trends need to hold
# DR is default
# ToDo: Check if time-varying covariates such as pop_density evolve bc of treatment, e.g. through predicting pop_density or use net migration


#### Formula for Conditional Parallel Trends ####
covariates_formula <- ~ pop_density + east_ger

# Keeps periods from 1994:
#covariates_formula <- ~ pop_density + share_fem + tax_rev
# Keeps periods from 1998:
#covariates_formula <- ~ pop_density + share_fem + tax_rev + hinc


#### Check Propensity Score Distribution between Treated and Control for Conditional Parallel Trends ####
ps_model <- glm(ifelse(seq_group > 0, 1, 0) ~ pop_density + east_ger,
                data = df_did_ready, family = binomial(), na.action = na.exclude)

df_did_ready$pscore <- predict(ps_model, type = "response")

ggplot(df_did_ready, aes(x = pscore, fill = factor(ifelse(seq_group > 0, 1, 0)))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("steelblue", "tomato"),
                    labels = c("Control", "Treated"),
                    name = "") +
  labs(title = "Overlap Check: Propensity Score Distribution",
       x = "Propensity Score", y = "Density") +
  theme_minimal()


#### Function for CS-DiD estimation ####
run_cs_did <- function(outcome, data, covariates_formula = NULL) {
  message(paste("Running CS-DiD for:", outcome))
  # Estimate (ATT(g,t))
  atts <- att_gt(
    yname = outcome,
    tname = "seq_time",
    idname = "ags",
    gname = "seq_group",
    xformla = covariates_formula,
    data = data,
    panel = TRUE,
    allow_unbalanced_panel = TRUE,
    clustervars = "ags",
    control_group = "notyettreated",
    anticipation = 0,
    bstrap      = TRUE,
    biters      = 100000
  )
  # Aggregate into event study
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}


#### Estimation for all outcomes ####
did_results <- map(outcome_vars, ~run_cs_did(.x, df_did_ready, covariates_formula = covariates_formula))
names(did_results) <- outcome_vars


#### Create event study plots for all outcomes ####
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_results[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  
  return(p)
})

grid.arrange(grobs = plot_list, ncol = 3)


#### Create cohort-specific plots for one outcome ####

# Extract ATT(g,t) estimates
att_turnout <- did_results$gruene$atts

# Build df from att_gt object
att_df <- data.frame(
  group    = att_turnout$group,
  time     = att_turnout$t,
  att      = att_turnout$att,
  se       = att_turnout$se
) %>%
  mutate(
    event_time = time - group,    # relative time to treatment
    ci_low  = att - 1.96 * se,
    ci_high = att + 1.96 * se,
    cohort  = factor(paste("Cohort", group))
  )

# Plot: one panel per cohort, x-axis = event time
ggplot(att_df, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(
    title    = "ATT(g,t) by Treatment Cohort",
    subtitle = "Each panel shows one treatment cohort; red line = treatment onset",
    x        = "Event Time (periods relative to treatment)",
    y        = "ATT"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))



###################### (Conditional) parallel trends for TWFE and staggered adoption ######################
# Honest DiD
# Check CS estimates: How strong are violations of pre-trends?
# Could add interacted linear trends as robustness check in CS-DiD or residualise
# Report breakdown value M


###################### Staggered DiD with continuous, absorbing treatment (Callaway, Goodman-Bacon & Sant'Anna, 2025) ######################

# Problem 1: Stronger parallel trends assumption (which cannot be tested) OR bias term
# Problem 2: No change in dose once treated
# Potentially address problem 2 through binning?
# Also: Not ATT but Average Causal Response

###################### Staggered DiD with Spatial Spillover (Butts, 2021) ######################
# Let's see if this rabbit hole is worth it to go down

###################### Staggered DiD with non-absorbing treatment ######################
# Is it non-absorbing?
# Or rather:
# - Is treatment multiple times? This specification would lead to a growth of potential outcomes
# - The many different potential outcomes make weighting, comparisons and estimates very hard to interpret
# - One could make assumption "treatment effect fades after 5 years"
# - Or "first treatment is only treatment" and then put dummy for subsequent treatment (see Bailey & Goodman-Bacon)
# - Or put unit multiple times in their data set 

# Possible additional tests with fect package: "no pretrend test" and "no carryover effect test"
# Need to test/assess exogeneity assumption, control for lagged vote share (potentially)

#### Controls ####

control_vars <- " ~ treat_nonabsorbing + pop_density + share_fem + hinc + share_foreign + tax_rev"

#### fect estimation function ####
run_fect_nonabs <- function(outcome, controls, data) {
  message(paste("Running fect (non-absorbing) for:", outcome))
  formula <- as.formula(paste0(outcome, controls))
  fit <- fect(
    formula,
    data    = data,
    index   = c("ags", "election_year"),
    method  = "ife",
    force   = "two-way",
    se      = TRUE,
    nboots  = 1000,
    min.T0  = 1
  )
  return(fit)
}

#### Run for all outcomes ####
fect_nonabs_results <- map(outcome_vars, ~run_fect_nonabs(.x, control_vars, df_did_ready))
names(fect_nonabs_results) <- outcome_vars


#### Plot results ####
grid.arrange(
  plot(fect_nonabs_results[[outcome_vars[1]]], main = outcome_vars[1]),
  plot(fect_nonabs_results[[outcome_vars[2]]], main = outcome_vars[2]),
  plot(fect_nonabs_results[[outcome_vars[3]]], main = outcome_vars[3]),
  plot(fect_nonabs_results[[outcome_vars[4]]], main = outcome_vars[4]),
  plot(fect_nonabs_results[[outcome_vars[5]]], main = outcome_vars[5]),
  plot(fect_nonabs_results[[outcome_vars[6]]], main = outcome_vars[6]),
  plot(fect_nonabs_results[[outcome_vars[7]]], main = outcome_vars[7]),
  plot(fect_nonabs_results[[outcome_vars[8]]], main = outcome_vars[8]),
  plot(fect_nonabs_results[[outcome_vars[9]]], main = outcome_vars[9]),
  ncol = 3,
  top = textGrob("FECT Absorbing DiD - All Outcomes",
                 gp = gpar(fontsize = 16, fontface = "bold"))
)


