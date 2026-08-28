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

# Overlay Plot
create_overlay_plot <- function(
    df_all, 
    outcome_vars = unique(df_all$outcome),
    outcome_labels = c(
      "turnout"           = "Voter Turnout",
      "cdu_csu"           = "CDU/CSU",
      "spd"               = "SPD",
      "gruene"            = "Greens",
      "afd"               = "AfD",
      "current_incumbent" = "Current Incumbent"
    ),
    est_cols   = c("CS (did)" = "#2b5c8f", "dCDH (did_multiplegt_dyn)" = "#fd7107"),
    est_shapes = c("CS (did)" = 16,        "dCDH (did_multiplegt_dyn)" = 17),
    est_lines  = c("CS (did)" = "solid",   "dCDH (did_multiplegt_dyn)" = "solid"),
    ncol = 2
) {
  # Build subplots
  plot_list <- map(outcome_vars, function(var) {
    ggplot(df_all %>% filter(outcome == var), 
           aes(x = e, y = estimate, color = estimator, shape = estimator, linetype = estimator)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
      geom_vline(xintercept = -1, linetype = "dotted", color = "gray70", linewidth = 0.4) +
      geom_pointrange(aes(ymin = coalesce(conf.low, estimate), ymax = coalesce(conf.high, estimate)),
                      position = position_dodge(0.35), size = 0.6) +
      geom_line(aes(group = estimator), position = position_dodge(0.35), linewidth = 0.7) +
      scale_x_continuous(breaks = seq(-6, 6, 1)) +
      scale_y_continuous(labels = function(x) format(x, scientific = FALSE, trim = TRUE)) +
      scale_color_manual(values = est_cols) +
      scale_shape_manual(values = est_shapes) +
      scale_linetype_manual(values = est_lines) +
      labs(title = outcome_labels[[var]] %||% var, x = NULL, y = NULL) +
      theme_classic() +
      theme(
        plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.text        = element_text(size = 15, color = "black"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 14, face = "bold"),
        legend.key.width = unit(1.4, "cm")
      )
  })
  
  # Extract shared legend safely via cowplot
  shared_legend <- cowplot::get_legend(plot_list[[1]])
  
  # Strip individual legends from plot list
  plot_list <- map(plot_list, ~ .x + theme(legend.position = "none"))
  
  # Assemble grid
  plots_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol)
  
  # Canvas layout
  cowplot::ggdraw() +
    cowplot::draw_plot(plots_grid, x = 0.045, y = 0.12, width = 0.955, height = 0.86) +
    cowplot::draw_label("Estimate", x = 0.015, y = 0.55, angle = 90, fontface = "bold", size = 14) +
    cowplot::draw_label("Event Time", x = 0.52, y = 0.065, fontface = "bold", size = 14) +
    cowplot::draw_plot(shared_legend, x = 0, y = 0.01, width = 1, height = 0.04)
}



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
  es   <- aggte(atts, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = 5)
  
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



###################### Can we fix PTA?  ######################

#### dCDH Controlled ####
dcdh_results_list_fix <- map(outcome_vars, function(var) {
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
    trends_nonparam = "bundesland",
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
df_dcdh_all <- bind_rows(dcdh_results_list_fix)

# Create Plot
create_overlay_plot(df_dcdh_all)


