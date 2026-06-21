###################### Main Analysis of Panel Data ######################

###################### Set the scene  ######################

#### Packages ####
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

library(fixest)
library(bacondecomp)
library(did)
library(HonestDiD)
library(contdid)
library(fect)


#### Helper functions ####
source("_support_functions/twfe_did_mult_outcome_wrapper.R")
source("_support_functions/honest_did.R")
source("_support_functions/HonestDiD_mult_outcome_wrapper.R")
source("_support_functions/mult_stage_loop.R")


#### Data ####
df_did_ready <- readRDS("_data/df_did_ready.rds")


#### Outcomes ####
outcome_vars <- c("turnout", "cdu", "csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")


###################### TWFE Event Study ######################

# Controls
covariates_str <- "+ pop_density + cat_cum_lag_wind_count_3km"

# Estimate TWFE event study
twfe_results <- map(outcome_vars, ~run_twfe_event_study(outcome = .x, data = df_did_ready, covariates = covariates_str))
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
covariate_str <- "~ treat_absorbing"

# Shorten data bc. vector memory was reached
did_short <- df_did_ready %>%
  filter(election_year < 2026)

# Run decomposition
bacon_results <- map(outcome_vars, ~run_bacon_decomposition(.x, did_short, "~ treat_absorbing")) %>% set_names(outcome_vars)

# Combine 2x2 decompositions into df
bacon_plot_data <- map_dfr(names(bacon_results), function(outcome) {
  df_2x2 <- bacon_results[[outcome]]
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
    strip.text = element_text(face = "bold", size = 11),
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
# Only conditional parallel trends need to hold, DR is default
# Checked ps overlap in descriptive_analyses script
# ToDo: Check if time-varying covariates such as pop_density evolve bc of treatment, e.g. through predicting pop_density or use net migration

#### Set CS options ####
# Set vars for conditional parallel trends with xformla
att_options_base <- list(
  data = df_did_ready,
  tname = "seq_time",
  idname = "ags",
  gname = "seq_group",
  xformla = ~ east_ger + pop_density,
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)

# Set vars for conditional parallel trends with xformla
# Keeps all periods: ~ east_ger + pop_density,
# Keeps periods from 1994: ~ pop_density + share_fem + tax_rev
# Keeps periods from 1998: ~ pop_density + share_fem + tax_rev + hinc

#### Estimate CS-DiD for all outcomes ####
did_results <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)


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
  group = att_turnout$group, time = att_turnout$t, att = att_turnout$att, se = att_turnout$se
  ) %>%
  mutate(event_time = time - group, ci_low  = att - 1.96 * se, ci_high = att + 1.96 * se,
         cohort  = factor(paste("Cohort", group))
         )

# Plot one panel per cohort, x-axis = event time
ggplot(att_df, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(title = "ATT(g,t) by Treatment Cohort", subtitle = "Each panel shows one treatment cohort; red line = treatment onset",
       x = "Event Time (periods relative to treatment)", y = "ATT") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))


###################### (Conditional) parallel trends for TWFE and staggered adoption ######################
# Check CS estimates: How strong are violations of pre-trends?
# Could add interacted linear trends as robustness check in CS-DiD or residualise
# Report breakdown value M

# Estimate CS-DiD for all outcomes using base_period = uniform
did_results_uni <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, base_period = "universal"), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)


#### Run HonestDiD estimations ####
# Smoothness
honest_smooth_results <- map(outcome_vars, ~run_honest_smoothness(.x, did_results_uni)) %>% set_names(outcome_vars)

# Relative Magnitude
honest_rm_results <- map(outcome_vars, ~run_honest_rm(.x, did_results_uni, mbar_seq = seq(0, 0.5, by = 0.05))) %>% set_names(outcome_vars)

#### Generate plots ####
# Smoothness Grid
generate_sensitivity_grid_plot(
  outcome_vars = outcome_vars, 
  results_list = honest_smooth_results, 
  type = "smooth", 
  ncol = 3)

# Relative Magnitudes Grid
generate_sensitivity_grid_plot(
  outcome_vars = outcome_vars, 
  results_list = honest_rm_results, 
  type = "rm", 
  ncol = 3)


###################### Staggered DiD with continuous, absorbing treatment (Callaway, Goodman-Bacon & Sant'Anna, 2025) ######################

# Notes (theoretical):
# - Continuous treatment required either strong parallel trends assumption (which cannot be tested) OR estimate is biased
# - No change in dose once treated (Potentially address this through binning?)
# Two treatment effects: Level vs. slope
# - Level treatment effect (ATT): Difference between untreated and treated under dose d
# - Causal response (ACRT): Difference in a units potential outcome under marginal increase of dose d
# --> Comparison between adjacent dose groups ≠ global effect (only with strong parallel trends assumption)

# Notes (practical):
# - contdid has mayor bugs described in this report: https://github.com/bcallaway11/contdid/issues/11
# - Can only implement constellations "slope + eventstudy" and "level + eventstudy" (under transformations)
# - Actually more interesting: "slope + dose" and "level + dose" --> But they don't work

# Clean data
df_clean <- df_did_ready %>%
  # Drop periods and group without within-period-dose variation and filter all NA
  filter(!is.na(treat_dose), !is.na(ags), !is.na(seq_time), !is.na(pop_density), election_year > 1994, !seq_group == 2)
# Balance panel
expected_periods <- n_distinct(df_clean$seq_time)
df_balanced <- df_clean %>%
  group_by(ags) %>%
  filter(n() == expected_periods) %>%
  ungroup()

# Slope + eventstudy
did_results_se <- map(outcome_vars, function(var) {
  df_temp <- df_balanced %>% filter(!is.na(.data[[var]]))
  res_cont_did <- cont_did(
    yname = var,
    tname = "seq_time",
    idname = "ags",
    dname = "treat_dose",
    gname = "seq_group",
    data = df_temp,
    target_parameter = "slope",
    aggregation     = "eventstudy",
    treatment_type  = "continuous",
    control_group   = "nevertreated",
    biters          = 1000,
    cband           = TRUE,
    num_knots       = 2,
    degree          = 5
  )
  return(list(res = res_cont_did))
}) %>% set_names(outcome_vars)


#### Create event study plots for all outcomes ####
plot_list_se <- map(outcome_vars, function(var) {
  ggcont_did(did_results_se[[var]]$res, type = "slope") +
    ggtitle(var) +
    theme_minimal()
})

grid.arrange(grobs = plot_list_se, ncol = 3)


# Level - eventstudy
# Needs treat_dose to be between 0 and 1
df_balanced <- df_balanced %>% mutate (treat_dose_squish = treat_dose / (1 + treat_dose))

did_results_le <- map(outcome_vars, function(var) {
  df_temp <- df_balanced %>% filter(!is.na(.data[[var]]))
  res_cont_did <- cont_did(
    yname = var,
    tname = "seq_time",
    idname = "ags",
    dname = "treat_dose_squish",
    gname = "seq_group",
    data = df_temp,
    target_parameter = "level",
    aggregation     = "eventstudy",
    treatment_type  = "continuous",
    control_group   = "nevertreated",
    biters          = 1000,
    cband           = TRUE,
    num_knots       = 2,
    degree          = 5
  )
  return(list(res = res_cont_did))
}) %>% set_names(outcome_vars)


#### Create event study plots for all outcomes ####
plot_list_le <- map(outcome_vars, function(var) {
  ggcont_did(did_results_le[[var]]$res, type = "slope") +
    ggtitle(var) +
    theme_minimal()
})

grid.arrange(grobs = plot_list_le, ncol = 3)



###################### New Specification of Treatment ######################
# Treatment is neither absorbing nor non-absorbing, but multiple times with different doses
# Idea: Investigate this multiple treatment effect through new specification
# 1.) Estimate effect of first WT in period X
# 2.) Filter for units treated in period X
# 3.) Use this sample to estimate effect of new WT in period X+1 (or X+2)
# ToDo:
# - How to define treatment? Every unit receiving one turbine or should it be binned, somehow accounting for dose?

# Outcomes (without AfD):
outcome_vars <- c("turnout", "cdu", "csu", "spd", "fdp", 
                  "linke_pds", "gruene", "current_incumbent")

# Set stages
stages <- list(
  list(year = 2002, suffix = "s1"), # seq_time = 4
  list(year = 2005, suffix = "s2"),
  list(year = 2009, suffix = "s3"),
  list(year = 2013, suffix = "s4"),
  list(year = 2017, suffix = "s5")
)

# Set CS options
att_opt_seq <- list(
  tname = "seq_time",
  idname = "ags",
  xformla = ~ east_ger + pop_density,
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)

# Display table of treated and untreated units by stage
generate_pipeline_summary(df_did_ready, stages, initial_filter_groups = c(4, 0))

# Run loop across all outcomes
seq_results <- map(outcome_vars, ~run_sequential_stages(.x, att_opt_seq, stages)) %>% set_names(outcome_vars)

# Create plot list
seq_plots <- generate_sequential_plots(seq_results)

# Combine into a 3x3 plot grid
final_grid <- wrap_plots(seq_plots, ncol = 3, nrow = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

final_grid


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


