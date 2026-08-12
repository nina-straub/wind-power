###################### Hyp1 and Hyp2 using dCDH ######################

###################### Set the scene  ######################

#### Packages ####
# install.packages(c("rlang", "S7"))
# install.packages("polars", repos = "https://rpolars.r-universe.dev")
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(did)
library(DIDmultiplegtDYN)
library(polars)

#### Helper functions ####
source("_support_functions/twfe_did_mult_outcome_wrapper.R")

#### Data ####
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


#### Outcomes and Controls ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")



###################### Hyp 1 ######################

# Option 0: Use absorbing, non-normalised, same switchers and no covariates for CS comparison
res_dcdh_0 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 0 (CS Comparison)",
  treatment       = "treat_absorbing",
  group           = "ags",
  time            = "seq_time",
  effects         = 4,
  placebo         = c(default = 4, afd = 0),
  cluster         = "ags",
  normalized      = F,
  same_switchers = T,
  same_switchers_pl = F
)

# Option 1: Most simple setup, no controls
res_dcdh_1 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 1 (Baseline)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = c(default = 5, afd = 3),
  placebo         = c(default = 4, afd = 1),
  cluster         = "ags",
  normalized      = TRUE
)

# Individual outcome estimates, e.g. gruene
summary(res_dcdh_1[["gruene"]])


# Option 2: Controls as in CS DiD
res_dcdh_2 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 2 (With Controls)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = c(default = 5, afd = 3),
  placebo         = c(default = 4, afd = 1),
  controls        = "pop_density",
  trends_nonparam = "bundesland",
  by              = "east_ger",
  cluster         = "ags",
  normalized      = TRUE
)


# Option 3: Only units with same trajectory
res_dcdh_3 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 3 (Only Complete Observations)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = 3,
  placebo         = c(default = 3, afd = 0),
  controls        = "pop_density",
  trends_nonparam = "east_ger",
  cluster         = "ags",
  normalized      = TRUE,
  same_switchers = TRUE,
  same_switchers_pl = TRUE
)


# Open Question:
# Normalisation yes/no
# Which controls and how (trends_nonparam vs. controls)?
# Same switchers/same_witchers_pl to ensure fully balanced panel across pre- and post periods? 

# Investigate further:
# Normalised weights
# Treatment paths


###################### Compare CS and dCDH from scratch ######################

# Absolute baseline: Try to force CS and dCDH to use sample as similar as possible
# I.e. create equivalence case described by dCDH (2025)


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
  control_group          = "notyettreated",
  base_period            = "universal",
  anticipation           = 0,
  bstrap                 = TRUE,
  biters                 = 1000,
  cband                  = FALSE
)

# Run CS DiD across all outcome variables
cs_results_list <- map(outcome_vars, function(var) {
  message(paste("Running CS-DiD for:", var))
  
  # Estimate group-time ATTs
  atts <- do.call(att_gt, c(list(yname = var), att_options_base))
  
  # Aggregate to dynamic event study
  es <- aggte(
    atts, 
    type      = "dynamic", 
    na.rm     = TRUE,
    min_e     = -5,
    max_e     = 5,
    balance_e = 4  # Holds cohort composition constant across e
  )
  
  # Return df
  data.frame(
    outcome   = var,
    e         = es$egt,
    estimate  = es$att,
    se        = es$se,
    estimator = "CS (did)"
  ) %>%
    mutate(
      conf.low  = estimate - 1.96 * se,
      conf.high = estimate + 1.96 * se
    )
}) %>% set_names(outcome_vars)

# Combine into a single CS data frame
df_cs_all <- bind_rows(cs_results_list)


#### dCDH baseline ####

# Extraction function to apply index shift (e = x - 1)
extract_dcdh_results <- function(dcdh_res, var_name) {
  # 1.) Placebos (x = -1, -2, ...)
  p_df <- if (!is.null(dcdh_res$results$Placebos)) {
    p_mat <- dcdh_res$results$Placebos
    data.frame(
      x        = -seq_len(nrow(p_mat)),
      estimate = p_mat[, "Estimate"],
      se       = p_mat[, "SE"]
    )
  } else { data.frame() }
  
  # 2.) Fixed Baseline Anchor (x = 0, normalized to 0)
  b_df <- data.frame(x = 0, estimate = 0, se = 0)
  
  # 3.) Post-Treatment Effects (x = 1, 2, ...)
  e_df <- if (!is.null(dcdh_res$results$Effects)) {
    e_mat <- dcdh_res$results$Effects
    data.frame(
      x        = seq_len(nrow(e_mat)),
      estimate = e_mat[, "Estimate"],
      se       = e_mat[, "SE"]
    )
  } else { data.frame() }
  
  # 4.) Combine and apply time shift: e = x - 1
  bind_rows(p_df, b_df, e_df) %>%
    mutate(
      e         = x - 1, # Maps x = 0 (baseline) to e = -1, matching CS timeline
      outcome   = var_name,
      estimator = "dCDH (did_multiplegt_dyn)",
      conf.low  = estimate - 1.96 * se,
      conf.high = estimate + 1.96 * se
    ) %>%
    select(outcome, e, estimate, se, conf.low, conf.high, estimator)
}

# Run dCDH across all outcome variables
dcdh_results_list <- map(outcome_vars, function(var) {
  message(paste("Running dCDH for:", var))
  
  res_dcdh <- did_multiplegt_dyn(
    df                   = df_did_ready,
    outcome              = var,
    group                = "ags",
    time                 = "seq_time",
    treatment            = "treat_absorbing",
    effects              = 5,
    placebo              = 4,
    cluster              = "ags",
    only_never_switchers = F, # Matches control_group = "notyettreated"
    same_switchers = T,
    same_switchers_pl    = F, # Matches balance_e = 4
    normalized           = FALSE, # Matches raw ATT outcome units,
    graph_off = T
  )
  
  extract_dcdh_results(res_dcdh, var)
}) %>% set_names(outcome_vars)

# Combine into a single dCDH data frame
df_dcdh_all <- bind_rows(dcdh_results_list)


#### Overlay Comparison Plots ####

# Combine CS and dCDH outputs
df_all_estimators <- bind_rows(df_cs_all, df_dcdh_all)

# Clean NA confidence bounds at the baseline period (e = -1)
df_all_estimators <- df_all_estimators %>%
  mutate(
    conf.low  = ifelse(is.na(conf.low), estimate, conf.low),
    conf.high = ifelse(is.na(conf.high), estimate, conf.high)
  )

# Generate grid of comparison plots
plot_list <- map(outcome_vars, function(var) {
  df_plot <- df_all_estimators %>% filter(outcome == var)
  
  ggplot(df_plot, aes(x = e, y = estimate, color = estimator, shape = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = -1, linetype = "dotted", color = "gray60") +
    geom_pointrange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = 0.3)
    ) +
    geom_line(aes(group = estimator), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = seq(-4, 4, 1)) +
    scale_color_manual(values = c("CS (did)" = "#2b5c8f", "dCDH (did_multiplegt_dyn)" = "#d95f02")) +
    labs(
      title = paste("Outcome:", var),
      x = "Event Time (e = 0 is first treated period; e = -1 is baseline)",
      y = "Estimate",
      color = "Estimator",
      shape = "Estimator"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 11)
    )
})

# Outcome grid plot
grid.arrange(grobs = plot_list, ncol = 3)


