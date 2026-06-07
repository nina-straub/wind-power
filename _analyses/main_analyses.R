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

#### Transform Treatment Variable ####

# Pre-process with sequential time and group indices
# Create sorted vector of election years
unique_years <- sort(unique(df_panel_cov$election_year))
# Create df mapping year to an index (1 to 10)
year_lookup <- data.frame(
  calendar_year = unique_years,
  seq_time = 1:length(unique_years)
)

# Reconstruct treat-var for did package, add incumbent variable and east Germany dummy
df_did_ready <- df_panel_cov %>%
  # Create sequential time index (1-10)
  left_join(year_lookup, by = c("election_year" = "calendar_year")) %>%
  # Create sequential group index for treatment year
  left_join(year_lookup, by = c("treat_absorbing_cs" = "calendar_year"), suffix = c("", "_group")) %>%
  mutate(
    # Create treatment group with sequential group index. If never treated, it stays 0.
    seq_group = ifelse(is.na(treat_absorbing_cs) | treat_absorbing_cs == 0, 0, seq_time_group),
    # Create current incumbent outcome
    current_incumbent = case_when(
      election_year %in% c(1990, 1994, 1998, 2013) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(fdp, 0),
      election_year %in% c(2002, 2005) ~ coalesce(spd, 0) + coalesce(gruene, 0),
      election_year %in% c(2009, 2017, 2021) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(spd, 0),
      election_year == 2025 ~ coalesce(gruene, 0) + coalesce(spd, 0) + coalesce(fdp, 0),
      TRUE ~ NA_real_),
    # Create "other parties" election outcome
    others = 1 - (coalesce(cdu, 0) + coalesce(afd, 0) + coalesce(csu, 0)
                  + coalesce(spd, 0) + coalesce(fdp, 0) + coalesce(linke_pds, 0)
                  + coalesce(gruene, 0)),
    # Add east Germany dummy
    east_ger = ifelse(substr(ags, 1, 2) %in% c("12", "13", "14", "15", "16"), 1, 0),
    # Ensure ID variable is numeric for the did package
    ags = as.numeric(ags)
  ) %>%
  # Filtering for certain states
  # filter(str_starts(ags, "12")) %>%
  # Filter for certain time period
  # filter(!election_year < 1990)
  as.data.frame()


###################### Vanilla TWFE ######################

# Filter outcomes
ols_outcomes <- c("spd", "cdu", "gruene", "afd")

# Create a clean post-treatment dummy for 2025
df_did_ready <- df_did_ready %>%
  mutate(post_2025 = ifelse(election_year == 2025, 1, 0))

# Run the TWFE models simultaneously for all 4 outcomes
# - Unit FE: ags
# - Time FE: election_year
# - Clustered SEs: clustered at the unit level (ags) by default in fixest when using FEs
twfe_models <- feols(
  c(spd, cdu, gruene, afd) ~ treat_nonabsorbing : post_2025 + 
    pop_density + 
    cat_cum_lag_wind_count_3km 
  | ags + election_year, 
  data = df_did_ready, cluster = ~ags
)

# Display a clean summary table of the results
etable(twfe_models, keep = "treat_nonabsorbing")

# 1. Estimate the Event Study model 
# i(election_year, treat_nonabsorbing, ref = 2021) creates leads and lags automatically
event_study_models <- feols(spd ~ i(election_year, treat_nonabsorbing, ref = 2021) + 
                              pop_density + 
                              cum_lag_wind_count_3km 
                            | ags + election_year, 
                            data = df_did_ready
)

# 2. Visualize the pre-trends and dynamic effects for all 4 outcomes
# Look for pre-treatment coefficients (years before 2025) that are tightly clustered around 0
iplot(event_study_models, main = "Parallel Trends Check (Ref: 2021)")


###################### Goodman-Bacon Decomposition (Goodman-Bacon, 2021) ######################
# bacondecomp


###################### Staggered DiD with binary, absorbing treatment (Callaway & Sant'Anna, 2021) ######################
# Only conditional parallel trends need to hold
# DR is default
# Need to check if time-varying covariates such as pop_density evolve bc of treatment, e.g. through predicting pop_density or use net migration

# 1. Define your variables
# List of all the vote share and turnout outcomes you want to analyze
outcome_vars <- c("turnout", "cdu", "csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent", "others")
# Pre-treatment covariates formula
# Keeps all periods:
covariates_formula <- ~ pop_density + east_ger
# Keeps periods from 1994:
covariates_formula <- ~ pop_density + share_fem + tax_rev
# Keeps periods from 1998:
covariates_formula <- ~ pop_density + share_fem + tax_rev + hinc


# 3. Create a function to run DiD and Event Study for a single outcome
run_cs_did <- function(outcome, data, covariates_formula = NULL) {
  message(paste("Running CS-DiD for:", outcome))
  
  # Estimate group-time average treatment effects (ATT(g,t))
  # Note: If panel is unbalanced, change allow_unbalanced_panel = TRUE
  atts <- att_gt(
    yname = outcome,
    tname = "seq_time",
    idname = "ags",
    gname = "seq_group",
    xformla = covariates_formula,
    data = data,
    panel = TRUE,
    allow_unbalanced_panel = TRUE, # Recommended for real-world election panels
    control_group = "notyettreated", # Alternatively "nevertreated"
    anticipation = 0
  )
  
  # Aggregate into an Event Study (Dynamic Effects)
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  
  return(list(atts = atts, es = es))
}

# 4. Loop through all outcomes using purrr::map
did_results <- map(outcome_vars, ~run_cs_did(.x, df_did_ready, covariates_formula = covariates_formula))
names(did_results) <- outcome_vars

# 5. Results

# Summary for Turnout
summary(did_results$turnout$es)

# Create event study plots for all outcomes
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_results[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  
  return(p)
})

grid.arrange(grobs = plot_list, ncol = 3)

# Plot the Event Study for Turnout
ggdid(did_results$turnout$es) + 
  ggtitle("Event Study: Effect of Treatment on Voter Turnout") +
  theme_minimal()


# Extract ATT(g,t) estimates for turnout
att_turnout <- did_results$spd$atts

# Build a tidy dataframe from the att_gt object
att_df <- data.frame(
  group    = att_turnout$group,   # treatment cohort (seq_group value)
  time     = att_turnout$t,       # calendar time (seq_time value)
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
    title    = "ATT(g,t) by Treatment Cohort — Turnout",
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


###################### Staggered DiD with contiunuous, absorbing treatment (Callaway, Goodman-Bacon & Sant'Anna, 2025) ######################

# Problem 1: Stronger parallel trends assumption (which cannot be tested) OR bias term
# Problem 2: No change in dose once treated
# Potentially address problem 2 through binning?


###################### Staggered DiD with non-absorbing treatment ######################
# Is it non-abosbring?
# Or rather:
# - Is treatment multiple times? This specification would lead to a growth of potential outcomes
# - The many different potetnial outcome make weighting, comparisons and estimates very hard to interpret
# - One could make asusmption "treatment effect fades after 5 years"
# - Or "first treatment is only treatment" and then put dummy for subsequent treatment (see Bailey & Goodman-Bacon)
# - Or put unit multiple times in their data set 

outcome_vars <- c("turnout", "cdu", "spd", "fdp", 
                  "linke_pds", "gruene", "cdu_csu")

# ========================
# FECT Non-Absorbing Setup
# ========================

# Prepare non-absorbing treatment data
main_data_nonabs <- df_panel_cov %>%
  filter(!is.na(treat_nonabsorbing), !is.na(ags), !is.na(election_year)) %>%
  mutate(
    cat_cum_lag_wind_count_3km = ifelse(
      election_year == 1990 & is.na(cat_cum_lag_wind_count_3km),
      0,
      cat_cum_lag_wind_count_3km
    ),
    cum_lag_wind_count_3km = ifelse(
      election_year == 1990 & is.na(cum_lag_wind_count_3km),
      0,
      cum_lag_wind_count_3km
    )
  )
  

# Remove always-treated units
treatment_share <- aggregate(treat_nonabsorbing ~ ags, data = main_data_nonabs, mean)
always_treated_ids <- treatment_share$ags[treatment_share$treat_nonabsorbing == 1]
main_data_nonabs <- subset(main_data_nonabs, !(ags %in% always_treated_ids))


# ========================
# fect runner function
# ========================

run_fect_nonabs <- function(outcome, data) {
  message(paste("Running fect (non-absorbing) for:", outcome))
  
  formula <- as.formula(paste0(
    outcome, " ~ treat_absorbing + pop_density + cat_cum_lag_wind_count_3km"
  ))
  
  fit <- fect(
    formula,
    data    = data,
    index   = c("ags", "election_year"),
    method  = "ife",
    force   = "two-way",
    se      = TRUE,
    nboots  = 100,
    min.T0  = 1
  )
  
  return(fit)
}

# ========================
# Run for all outcomes
# ========================

fect_nonabs_results <- map(outcome_vars, ~run_fect_nonabs(.x, main_data_nonabs))
names(fect_nonabs_results) <- outcome_vars

# ========================
# Plot all results
# ========================

grid.arrange(
  plot(fect_nonabs_results[[outcome_vars[1]]], main = outcome_vars[1]),
  plot(fect_nonabs_results[[outcome_vars[2]]], main = outcome_vars[2]),
  plot(fect_nonabs_results[[outcome_vars[3]]], main = outcome_vars[3]),
  plot(fect_nonabs_results[[outcome_vars[4]]], main = outcome_vars[4]),
  plot(fect_nonabs_results[[outcome_vars[5]]], main = outcome_vars[5]),
  plot(fect_nonabs_results[[outcome_vars[6]]], main = outcome_vars[6]),
  ncol = 3,
  top = textGrob("FECT Non-Absorbing DiD - All Outcomes",
                 gp = gpar(fontsize = 16, fontface = "bold"))
)


