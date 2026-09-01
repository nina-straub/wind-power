###################### Compare CS and dCDH from scratch ######################

# Step 1: Naive implementation of CS and dCDH
# Step 2: Create equivalence case described by dCDH (2025) as sanity check
# Step 3: Built up complexity

###################### Set the scene  ######################

#### Packages ####
# install.packages(c("rlang", "S7"))
# install.packages("polars", repos = "https://rpolars.r-universe.dev")
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(did)
library(DIDmultiplegtDYN)
library(polars)

#### Helper functions ####
source("_support_functions/twfe_did_mult_outcome_wrapper.R")


#### Data ####
df_did_ready <- readRDS("_data/df_did_ready.rds")


###################### Step 1: Naive implementation of both  ######################

#### Outcomes ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "afd", "current_incumbent")


#### CS Naive Implementation ####
att_options_naive <- list(
  data = df_did_ready,
  tname = "seq_time",
  idname = "ags",
  gname = "seq_group",
  xformla = ~ east_ger + pop_density,
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  base_period = "varying",
  control_group = "nevertreated",
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)

cs_results_list_naive <- map(outcome_vars, function(var) {
  message(paste("Running CS-DiD for:", var))
  
  atts <- do.call(att_gt, c(list(yname = var), att_options_naive))
  es   <- aggte(atts, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = 4)
  
  # Extract uniform critical value directly
  crit_val <- es$crit.val.egt
  
  data.frame(
    outcome   = var,
    e         = es$egt,
    estimate  = es$att.egt,
    se        = es$se.egt,
    estimator = "CS (did)"
  ) %>%
    mutate(
      conf.low  = estimate - crit_val * se,
      conf.high = estimate + crit_val * se
    )
}) %>% set_names(outcome_vars)

df_cs_all <- bind_rows(cs_results_list_naive)


#### dCDH Naive Implementation ####
dcdh_results_list_naive <- map(outcome_vars, function(var) {
  message(paste("Running dCDH for:", var))
  
  # Dynamic placebo/effects window: 3 for 'afd', standard for others
  plc_val <- if (var == "afd") 1 else 4
  eff_val <- if (var == "afd") 3 else 5
  
  res_dcdh <- did_multiplegt_dyn(
    df            = df_did_ready,
    outcome        = var,
    treatment       = "cum_wind_count_3km",
    group           = "ags",
    time            = "seq_time",
    effects         = eff_val,
    placebo         = plc_val,
    controls        = "pop_density",
    trends_nonparam = "east_ger",
    cluster         = "ags",
    normalized      = T,
    only_never_switchers = F,
    same_switchers = F,
    same_switchers_pl = F,
    graph_off = T
  )
  
  extract_dcdh_results(res_dcdh, var)
}) %>% set_names(outcome_vars)

# Combine into a single dCDH data frame
df_dcdh_all <- bind_rows(dcdh_results_list_naive)


#### Overlay Comparison Plots ####

# Combine CS and dCDH outputs
df_all_estimators <- bind_rows(df_cs_all, df_dcdh_all)

# Create Plot
create_overlay_plot(df_all_estimators)

# Create Plot for nAVSQ only
create_overlay_plot(df_dcdh_all)


###################### Step 2: Create baseline  ######################

#### Outcomes and Controls ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "afd", "current_incumbent")

#### CS DiD baseline ####
att_options_base <- list(
  data                   = df_did_ready,
  tname                  = "seq_time",
  idname                 = "ags",
  gname                  = "seq_group",
  xformla                = ~ 1,
  panel                  = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars            = "ags",
  control_group          = "nevertreated",
  base_period            = "universal",
  anticipation           = 0,
  bstrap                 = TRUE,
  biters                 = 1000,
  cband                  = FALSE
)

# Run CS DiD across all outcome variables
cs_results_list_baseline <- map(outcome_vars, function(var) {
  message(paste("Running CS-DiD for:", var))
  
  # Estimate group-time ATTs
  atts <- do.call(att_gt, c(list(yname = var), att_options_base))
  
  # Match event-study windows to dCDH (placebo & effects)
  min_e_val <- if (var == "afd") -1 else -5  # Matches dCDH placebo = 1 vs 4
  max_e_val <- if (var == "afd")  3 else  5  # Matches dCDH effects = 3 vs 5
  bal_val   <- if (var == "afd")  2 else  4  # Dynamic balance window
  
  # Aggregate to dynamic event study
  es <- aggte(
    atts, 
    type      = "dynamic", 
    na.rm     = TRUE,
    min_e     = min_e_val,
    max_e     = max_e_val,
    balance_e = bal_val  # Holds cohort composition constant across e
  )
  
  # Extract uniform critical value directly
  crit_val <- es$crit.val.egt
  
  # Return df
  data.frame(
    outcome   = var,
    e         = es$egt,
    estimate  = es$att,
    se        = es$se,
    estimator = "CS (did)"
  ) %>%
    mutate(
      conf.low  = estimate - crit_val * se,
      conf.high = estimate + crit_val * se
    )
}) %>% set_names(outcome_vars)

# Combine into a single CS data frame
df_cs_all <- bind_rows(cs_results_list_baseline)



#### dCDH baseline ####


# Run dCDH across all outcome variables
dcdh_results_list_baseline <- map(outcome_vars, function(var) {
  message(paste("Running dCDH for:", var))
  
  # Dynamic placebo/effects window: 3 for 'afd', standard for others
  plc_val <- if (var == "afd") 1 else 4
  eff_val <- if (var == "afd") 3 else 5
  
  res_dcdh <- did_multiplegt_dyn(
    df                   = df_did_ready,
    outcome              = var,
    group                = "ags",
    time                 = "seq_time",
    treatment            = "cum_wind_count_3km",
    effects              = eff_val,
    placebo              = plc_val,
    cluster              = "ags",
    only_never_switchers = T, # Matches control_group = "nevertreated"
    same_switchers = T,
    same_switchers_pl    = F, # Matches balance_e = 4
    normalized           = FALSE, # Matches raw ATT outcome units,
    graph_off = T
  )
  
  extract_dcdh_results(res_dcdh, var)
}) %>% set_names(outcome_vars)

# Combine into a single dCDH data frame
df_dcdh_all <- bind_rows(dcdh_results_list_baseline)


#### Overlay Comparison Plots ####

# Combine CS and dCDH outputs
df_all_estimators <- bind_rows(df_cs_all, df_dcdh_all)

# Create Plot
create_overlay_plot(df_all_estimators)


