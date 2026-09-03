###################### Hypothesis 2: Costs and Benefits ######################

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


#### Data ####
df_did_ready <- readRDS("_data/df_did_ready.rds")


#### Outcomes & CS Options ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "far_right", "current_incumbent")

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

att_options_base_2 <- list(
  tname = "seq_time",
  idname = "ags",
  gname = "second_treat_time",
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)


###################### Compare high vs. low benefit units' reactions to later treatment ######################

#### Construct data ####
# Identify treated units and winners/losers during first treatment (by inc_tax)
indic_cb <- df_did_ready %>%
  filter(seq_group > 0, seq_time >= seq_group - 1, seq_time <= seq_group + 1) %>%
  group_by(ags, treat_absorbing) %>%
  summarise(avg_inctax = mean(hinc, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = treat_absorbing, values_from = avg_inctax) %>%
  mutate(
    economic_benefit = `1` - `0`,
    # Define winners vs losers based on median growth or a 0-threshold
    shock_type = if_else(economic_benefit > median(economic_benefit, na.rm = TRUE), "winner", "loser"),
    # New shock type with 25th percentile breaks
    shock_type_4way = case_when(
      economic_benefit <= quantile(economic_benefit, 0.25, na.rm = TRUE) ~ "big_loser",
      economic_benefit <= quantile(economic_benefit, 0.50, na.rm = TRUE) ~ "loser",
      economic_benefit <= quantile(economic_benefit, 0.75, na.rm = TRUE) ~ "winner",
      TRUE                                                               ~ "big_winner"
    )
  )

df_cb <- df_did_ready %>%
  left_join(
    indic_cb %>% select(ags, economic_benefit, shock_type, shock_type_4way), by = "ags") %>%
  group_by(ags) %>%
  arrange(seq_time, .by_group = TRUE) %>% # Ensure strict chronological order
  mutate(
    # Get the time period for the first and second treatments
    n_treats = sum(treat_nonabsorbing, na.rm = TRUE),
    first_treat_time = if_else(n_treats >= 1, seq_time[treat_nonabsorbing == 1][1], NA_real_),
    # Set second_treat_time = 0 for units with 0 treatments (did package convention)
    second_treat_time = if_else(n_treats >= 2, seq_time[treat_nonabsorbing == 1][2], 0),
    # Calculate the gap. consecutive periods (e.g., 2015 and 2016) = 1.
    # We need a gap >= 2 (e.g., 2015 and 2017, leaving 2016 as a '0' break).
    treatment_gap = if_else(n_treats >= 2, second_treat_time - first_treat_time, NA_real_)
  ) %>%
  # Retain:
  # 1. Multi-treated units (n_treats >= 2) after their first treatment with gap > 1
  # 2. Truly never-treated units (n_treats == 0)
  # Exclude: Single-treated units (n_treats == 1)
  filter(
    (n_treats >= 2 & treatment_gap > 1 & seq_time > first_treat_time) |
      (n_treats == 0)
  ) %>%
  ungroup()

# Extract pure controls (units never receiving any treatment)
never_treated_pool <- df_cb %>% filter(n_treats == 0)

# Build subgroup datasets: Multi-treated units split by 1st treatment shock type + pure controls
df_cb_los  <- bind_rows(df_cb %>% filter(shock_type == "loser"), never_treated_pool) %>% distinct()
df_cb_win  <- bind_rows(df_cb %>% filter(shock_type == "winner"), never_treated_pool) %>% distinct()
df_cb_blos <- bind_rows(df_cb %>% filter(shock_type_4way == "big_loser"), never_treated_pool) %>% distinct()
df_cb_bwin <- bind_rows(df_cb %>% filter(shock_type_4way == "big_winner"), never_treated_pool) %>% distinct()




#### Estimation ####
did_los <- run_csdid_pipeline(df_cb_los, outcome_vars, options = att_options_base_2, 
                              label = "Loser", gname = "second_treat_time", formula = ~ 1)
did_win <- run_csdid_pipeline(df_cb_win, outcome_vars, options = att_options_base_2,
                              label = "Winner", gname = "second_treat_time", formula = ~ 1)

did_blos <- run_csdid_pipeline(df_cb_blos, outcome_vars, options = att_options_base_2,
                               label = "Big Loser", gname = "second_treat_time", formula = ~ pop_density)
did_bwin <- run_csdid_pipeline(df_cb_bwin, outcome_vars, options = att_options_base_2,
                               label = "Big Winner", gname = "second_treat_time", formula = ~ pop_density)


#### Extract Loser estimates from did_los pipeline output ####
df_cs_los <- map_dfr(outcome_vars, function(var) {
  es <- did_los[[var]]$es
  crit_val <- es$crit.val.egt
  
  data.frame(
    outcome   = var,
    e         = es$egt,
    estimate  = es$att.egt,
    se        = es$se.egt,
    estimator = "Loser"
  ) %>%
    mutate(
      conf.low  = estimate - crit_val * se,
      conf.high = estimate + crit_val * se
    )
})

# Extract Winner estimates from did_win pipeline output
df_cs_win <- map_dfr(outcome_vars, function(var) {
  es <- did_win[[var]]$es
  crit_val <- es$crit.val.egt
  
  data.frame(
    outcome   = var,
    e         = es$egt,
    estimate  = es$att.egt,
    se        = es$se.egt,
    estimator = "Winner"
  ) %>%
    mutate(
      conf.low  = estimate - crit_val * se,
      conf.high = estimate + crit_val * se
    )
})

#### Combine into a df for plotting ####
df_all <- bind_rows(df_cs_win, df_cs_los)

create_overlay_plot(df_all)


###################### Construct index ######################

# Omitted


###################### Simple Approach: Any changes in MeckPomm post 2016? ######################
# Filter df and split into early/late df's
list2env(
  df_did_ready %>%
    mutate(ags = str_pad(as.character(ags), width = 8, side = "left", pad = "0")) %>%
    filter(str_starts(ags, "16")) %>%
    mutate(ags = as.numeric(ags)) %>%
    split(.$election_year < 2016) %>%
    setNames(c("df_meckpomm_late", "df_meckpomm_early")),
  envir = .GlobalEnv
)

# Estimation
did_mp_early <- run_csdid_pipeline(df_meckpomm_early, outcome_vars, label = "MeckPomm Early", formula = ~ pop_density)
did_mp_late <- run_csdid_pipeline(df_meckpomm_late, outcome_vars, label = "MeckPomm Late", formula = ~ pop_density)

# Okay, no effects
# --> But honestly, in the East there were no effects in the first place, so what to expect?






# ============================================================
# Alternative:
# CDE via Sequential G-Estimation (Panel-Safe & Optimized)
# Outcome: vote shares | Mediator: hinc
# Use Blackwell et al.'s extension of DirectEffect
# ============================================================

