###################### Subgroup Analysis of Panel Data ######################

###################### Set the scene  ######################

#### Packages ####
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(did)


#### Helper functions ####
source("_support_functions/twfe_did_mult_outcome_wrapper.R")


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


#### Split data by subgroups ####

# East vs. West Germany (wo Bavaria)
df_east <- df_did_ready %>% filter(east_ger == 1)
df_west <- df_did_ready %>% filter(east_ger == 0)

# Early vs. Late Germany (polarisation hypothesis)
df_early <- df_did_ready %>% filter(election_year < 2012)
df_late <- df_did_ready %>% filter(election_year > 2012)


# Interaction: Late x West Germany
df_west_late <- df_did_ready %>% filter(east_ger == 0, election_year > 2012)

# Units treated once
# Discard all units treated more than once
df_treat_once <- df_did_ready %>%
  group_by(ags) %>%
  filter(sum(treat_nonabsorbing, na.rm = TRUE) <= 1) %>%
  ungroup()


###################### Outcomes & CS/dCDH Options ######################

outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "afd", "current_incumbent")

# Drop AfD in Outcomes
outcome_vars_wo_afd <- c("turnout", "cdu_csu", "spd", "gruene", "far_right", "current_incumbent")


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


###################### East vs. West Germany - CS Estimator ######################
did_east_cs <- run_csdid_pipeline(df_east, outcome_vars, label = "East Germany", formula = ~ pop_density)
did_west_cs <- run_csdid_pipeline(df_west, outcome_vars, label = "West Germany", formula = ~ pop_density)

df_cs_east <- extract_cs_df(did_east_cs, "East")
df_cs_west <- extract_cs_df(did_west_cs, "West")

#### Combine into a df for plotting ####
df_all_region_cs <- bind_rows(df_cs_east, df_cs_west)

create_overlay_plot(df_all_region_cs)


###################### East vs. West Germany - dCDH Estimator ######################
did_east_dcdh <- run_dcdh_pipeline(df_east, outcome_vars, label = "East", options = dcdh_options_base)
did_west_dcdh <- run_dcdh_pipeline(df_west, outcome_vars, label = "West", options = dcdh_options_base)

# Combine into a single dCDH data frame
df_dcdh_all_east <- bind_rows(did_east_dcdh)
df_dcdh_all_west <- bind_rows(did_west_dcdh)

#### Combine into a df for plotting ####
df_all_region_dcdh <- bind_rows(df_dcdh_all_east, df_dcdh_all_west)

create_overlay_plot(df_all_region_dcdh)


###################### Early vs. Late Germany (Polarisation Hypothesis) - CS Estimator ######################

did_cs_early <- run_csdid_pipeline(df_early, outcome_vars_wo_afd, label = "Early", formula = ~ pop_density + east_ger)
did_cs_late  <- run_csdid_pipeline(df_late, outcome_vars_wo_afd, label = "Late", formula = ~ pop_density + east_ger)

df_cs_early <- extract_cs_df(did_cs_early, "Early")
df_cs_late <- extract_cs_df(did_cs_late, "Late")

#### Combine into a df for plotting ####
df_all_time_cs <- bind_rows(df_cs_early, df_cs_late)

create_overlay_plot(df_all_time_cs)


###################### Early vs. Late Germany (Polarisation Hypothesis) - dCDH Estimator ######################

did_dcdh_early <- run_dcdh_pipeline(df_early, outcome_vars_wo_afd, label = "Early",
                                    effects_default = 3,
                                    placebo_default = 1)
did_dcdh_late  <- run_dcdh_pipeline(df_late, outcome_vars_wo_afd, label = "Late",
                                    effects_default = 3,
                                    placebo_default = 1)

# Combine into a single dCDH data frame
df_dcdh_all_early <- bind_rows(did_dcdh_early)
df_dcdh_all_late <- bind_rows(did_dcdh_late)

#### Combine into a df for plotting ####
df_all_time_dcdh <- bind_rows(df_dcdh_all_early, df_dcdh_all_late)

create_overlay_plot(df_all_time_dcdh)


###################### Units treated once - CS Estimator ######################
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
did_cs_treat_once <- run_csdid_pipeline(df_treat_once, outcome_vars, label = "One-time Treatment")


###################### Units treated once - dCDH Estimator ######################
did_dcdh_treat_once <- run_dcdh_pipeline(df_treat_once, outcome_vars, label = "Early",
                                    effects_default = 4,
                                    placebo_default = 3)

df_all_once_dcdh <- bind_rows(did_dcdh_treat_once)

create_overlay_plot(df_all_once_dcdh)


