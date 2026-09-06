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
source("_support_functions/honest_did_cs_bridge.R")
source("_support_functions/HonestDiD_mult_outcome_wrapper.R")
source("_support_functions/mult_stage_loop.R")


###################### Data ######################
df_did_ready <- readRDS("_data/df_did_ready.rds")

# Add categorical variable "state" (bundesland)
df_did_ready <- df_did_ready %>%
  mutate(
    ags_char   = str_pad(as.character(ags), width = 8, side = "left", pad = "0"),
    bundesland = substr(ags_char, 1, 2)
  )

# Create static baseline strata for population density
df_did_ready <- df_did_ready %>%
  group_by(ags) %>%
  mutate(pop_density_1990 = pop_density[seq_time == min(seq_time)]) %>%
  ungroup() %>%
  mutate(
    pop_quartile_1990 = ntile(pop_density_1990, 4),
    state_pop_strata  = paste(bundesland, pop_quartile_1990, sep = "_")
  )

# Units treated once
# Discard all units treated more than once
df_treat_once <- df_did_ready %>%
  group_by(ags) %>%
  filter(sum(treat_nonabsorbing, na.rm = TRUE) <= 1) %>%
  ungroup()


#### Outcomes ####
outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "afd", "current_incumbent")


###################### Sensitivity Analysis with HonestDiD ######################
# Check CS estimates: How strong are violations of pre-trends?
# Could add interacted linear trends as robustness check in CS-DiD or residualise
# Report breakdown value M

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

# saveRDS(honest_smooth_results, '_results/hyp1_main_analyses/res_honest_smooth.rds')

# Relative Magnitude
honest_rm_results <- map(outcome_vars, ~run_honest_rm(.x, did_results_uni, mbar_seq = seq(0, 07., by = 0.07))) %>% set_names(outcome_vars)

# saveRDS(honest_rm_results, '_results/hyp1_main_analyses/res_honest_rm.rds')


#### Generate plots ####
# Smoothness Grid
generate_sensitivity_grid_plot(
  outcome_vars = outcome_vars, 
  results_list = honest_smooth_results, 
  type = "smooth", 
  ncol = 2)

# Relative Magnitudes Grid
generate_sensitivity_grid_plot(
  outcome_vars = outcome_vars, 
  results_list = honest_rm_results, 
  type = "rm", 
  ncol = 2)

# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_honest_smooth.png", plot = plot_honest_smooth, width = 12, height = 8, dpi = 300)
# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_honest_rm.png", plot = plot_honest_rm, width = 12, height = 8, dpi = 300)



###################### Units treated once - CS Estimator ######################

# Set CS options
# Set vars for conditional parallel trends with xformla
att_options_base <- list(
  tname = "seq_time",
  idname = "ags",
  gname = "seq_group",
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)

# Check distribution of units treated once, twice, etc.
df_did_ready %>%
  group_by(ags) %>%
  summarise(total_treat = sum(treat_nonabsorbing, na.rm = TRUE)) %>%
  count(total_treat)

# Check if 0 - 1 holds in filtered df
df_treat_once %>%
  group_by(ags) %>%
  summarise(total_treat = sum(treat_nonabsorbing, na.rm = TRUE)) %>%
  count(total_treat)

# Estimation loop
did_cs_once <- run_csdid_pipeline(df_treat_once, outcome_vars, label = "One-time Treatment")

df_cs_once <- extract_cs_df(did_cs_once, "Treated Once")
create_overlay_plot(df_cs_once)


###################### Units treated once - dCDH Estimator ######################
# Set dCDH option
dcdh_options_base <- list(
  treatment = "cum_wind_count_3km",
  group = "ags",
  time = "seq_time",
  trends_nonparam = "state_pop_strata",
  cluster = "ags",
  normalized = TRUE,
  only_never_switchers = FALSE,
  same_switchers = FALSE,
  same_switchers_pl = FALSE,
  graph_off = TRUE
)

did_dcdh_once <- run_dcdh_pipeline(df_treat_once, outcome_vars, label = "Treated Once",
                                   effects_default = 6,
                                   placebo_default = 3,
                                   placebo_afd = 1,
                                   effects_afd = 3)

df_once_dcdh <- bind_rows(did_dcdh_once)
create_overlay_plot(df_once_dcdh)



###################### Placebo Outcomes ######################
# Does treatment affect net migration or agricultural land size? If yes --> Signs for sorting/negative effects of treatment
# Does it affect education or share of female/foreign population? If yes --> Signs for sorting
outcomes_alt <- c("net_migration", "educ", "share_fem", "share_foreign")

#### With CS ####
did_alt <- run_csdid_pipeline(df_did_ready, outcomes_alt, label = "Test Assumptions", formula = ~ pop_density + east_ger)
df_cs_alt <- extract_cs_df(did_alt, "Alternative Outcomes")
create_overlay_plot(df_cs_alt)

#### With dCDH ####
did_dcdh_once <- run_dcdh_pipeline(df_did_ready, outcomes_alt, label = "Alternative Outcomes",
                                   effects_default = 6,
                                   placebo_default = 3,
                                   placebo_afd = 1,
                                   effects_afd = 3)

df_once_dcdh <- bind_rows(did_dcdh_once)
create_overlay_plot(df_once_dcdh)


# How are socioeconomic indicators affected?
# Could give idea on mechanisms --> If positive effect on tax/employment etc., positive perception of WP
outcomes_se <- c("inc_tax", "busi_tax", "tax_rev", "commute_balance", "dens_work", "purch_pow", "tourism", "hinc", "share_unemp", "share_emp")
did_se <- run_csdid_pipeline(df_did_ready, outcomes_se, label = "Socioeconomic Indicators", formula = ~ pop_density + east_ger)


###################### Anticipation ######################

# Set CS options
# Set vars for conditional parallel trends with xformla
att_options_base_anti <- list(
  tname = "seq_time",
  idname = "ags",
  gname = "seq_group",
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  anticipation = 1,
  bstrap = TRUE,
  biters = 1000
)

# Estimation loop
did_cs_anti <- run_csdid_pipeline(df_did_ready, outcome_vars, options = att_options_base_anti, label = "Anticipation")

df_cs_anti <- extract_cs_df(did_cs_anti, "Anticipation")
create_overlay_plot(df_cs_anti)


###################### TWFE Event Study ######################

# Controls
covariates_str <- "+ pop_density"

# Estimate TWFE event study
twfe_results <- map(outcome_vars, ~run_twfe_event_study(outcome = .x, data = df_did_ready, covariates = covariates_str))
names(twfe_results) <- outcome_vars

# Summary of outcomes
etable(map(twfe_results, ~ .x$model), keep = "time_to_treatment")


# Plot results
plot_twfe <- function() {
  par(mfrow = c(3, 2))
  for (i in 1:6) {
    iplot(twfe_results[[i]], 
          main = paste(outcome_vars[i], "(Ref: 1 Year Pre-Treatment)"))
  }
  par(mfrow = c(1, 1))
}

# To save it to a file:
# png("_results/hyp1_main_analyses/_figures/plot_twfe.png", width = 1200, height = 800)
plot_twfe()
dev.off()



###################### Goodman-Bacon Decomposition (Goodman-Bacon, 2021) ######################
# Running decomposition on all outcomes takes some time. Possibly parallelise code here with future.

# Covariates
covariate_str <- "~ treat_absorbing"

# Shorten data bc. vector memory was reached
did_short <- df_did_ready %>%
  filter(election_year > 2012)

# Run decomposition
bacon_results <- map(outcome_vars, ~run_bacon_decomposition(.x, did_short, "~ treat_absorbing")) %>% set_names(outcome_vars)

# saveRDS(bacon_results, '_results/hyp1_main_analyses/res_bacon_decomp.rds')

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
bacon_plot <- ggplot(bacon_plot_data, aes(x = weight, y = estimate, shape = type, col = type)) +
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

# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_bacon_decomp.png", plot = bacon_plot, width = 12, height = 8, dpi = 300)



###################### New Specification of Treatment ######################
# Treatment is neither absorbing nor non-absorbing, but multiple times with different doses
# Idea: Investigate this multiple treatment effect through new specification
# 1.) Estimate effect of first WT in period X
# 2.) Filter for units treated in period X
# 3.) Use this sample to estimate effect of new WT in period X+1 (or X+2)
# ToDo:
# - How to define treatment? Every unit receiving one turbine or should it be binned, somehow accounting for dose?

# Outcomes (without AfD):
outcome_vars_wo_afd <- c("turnout", "cdu_csu", "spd", "gruene", "far_right", "current_incumbent")

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
seq_results <- map(outcome_vars_wo_afd, ~run_sequential_stages(.x, att_opt_seq, stages)) %>% set_names(outcome_vars_wo_afd)

# saveRDS(seq_results, '_results/hyp1_main_analyses/res_seq_treat.rds')

# Create plot list
seq_plots <- generate_sequential_plots(seq_results)

# Combine into a 3x3 plot grid
final_grid <- wrap_plots(seq_plots, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

final_grid

# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_seq_treat.png", plot = final_grid, width = 12, height = 8, dpi = 300)






