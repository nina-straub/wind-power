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

#### Data ####
df_did_ready <- readRDS("_data/df_did_ready.rds")

#### For subgroup analysis ####

# Without and Only Bavaria
df_wo_bav <- df_did_ready %>% filter(!str_starts(as.character(ags), "9"))
df_bav <- df_did_ready %>% filter(str_starts(as.character(ags), "9"))
df_bav_10h <- df_wo_bav %>% filter(election_year > 2014)

# East vs. West Germany (wo Bavaria)
df_east <- df_did_ready %>% filter(east_ger == 1)
df_west <- df_wo_bav %>% filter(east_ger == 0)

# Early vs. Late Germany (polarisation hypothesis)
df_early_p <- df_did_ready %>% filter(election_year < 2012)
df_late_p <- df_did_ready %>% filter(election_year > 2012)

# Early vs. Late Germany (tax hypothesis)
df_early_t <- df_did_ready %>% filter(election_year <= 2009)
df_late_t <- df_did_ready %>% filter(election_year > 2009)

# Interaction: Late x West Germany
df_west_late <- df_wo_bav %>% filter(east_ger == 0, election_year > 2015)

# High vs. Low Density
df_hd <- df_did_ready %>% filter(pop_density > quantile(pop_density, 0.5, na.rm = TRUE),
                                 pop_density <= quantile(pop_density, 0.9, na.rm = TRUE))
df_ld <- df_did_ready %>% filter(pop_density < quantile(pop_density, 0.5, na.rm = TRUE))

# High vs. Low AfD
# Filter for municipalities with above-average AfD vote share
df_afd_h <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_afd = mean(afd, na.rm = TRUE)) %>%
                     filter(mean_afd > median(mean_afd, na.rm = TRUE)) %>%
                     pull(ags)))
# Filter for municipalities with below-average AfD vote share
df_afd_l <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_afd = mean(afd, na.rm = TRUE)) %>%
                     filter(mean_afd < median(mean_afd, na.rm = TRUE)) %>%
                     pull(ags)))

# High vs. Low Greens
# Filter for municipalities with above-average Greens vote share
df_gruene_h <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_gruene = mean(gruene, na.rm = TRUE)) %>%
                     filter(mean_gruene > median(mean_gruene, na.rm = TRUE)) %>%
                     pull(ags)))

# Filter for municipalities with below-average Greens vote share
df_gruene_l <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_gruene = mean(gruene, na.rm = TRUE)) %>%
                     filter(mean_gruene < median(mean_gruene, na.rm = TRUE)) %>%
                     pull(ags)))

# Distance Variation
# Create new seq_group over 4 distances
election_years <- sort(unique(df_did_ready$election_year))
distances <- c("0km", "1km", "3km", "5km", "10km")
df_distance <- df_did_ready %>%
  group_by(ags) %>%
  mutate(across(
    all_of(paste0("wind_count_", distances)),
    ~ suppressWarnings(ifelse(is.infinite(min(election_year[. > 0], na.rm = TRUE)),
                              0, min(election_year[. > 0], na.rm = TRUE))),
    .names = "tacs{str_remove(.col, 'wind_count_')}"
  )) %>%
  ungroup() %>%
  mutate(across(
    all_of(paste0("tacs", distances)),
    ~ ifelse(is.na(match(., election_years)), 0, match(., election_years)),
    .names = "seq_group_{str_remove(.col, 'tacs')}"
  ))

# Intensity Variation
# Create new treatment intensity variable over 5 intensities
thresholds <- c(1, 3, 5, 10, 15)
df_intensity <- reduce(thresholds, function(df, thresh) {
  df %>% 
    group_by(ags) %>% 
    mutate(
      !!paste0("treat_", thresh) := case_when(
        seq_group == 0 ~ 0,
        any(seq_time == seq_group & wind_count_3km >= thresh, na.rm = TRUE) ~ seq_group,
        TRUE ~ NA_real_
      )
    ) %>% 
    ungroup()
}, .init = df_did_ready)

# Units treated once
# Discard all units treated more than once
df_treat_once <- df_did_ready %>%
  group_by(ags) %>%
  filter(sum(treat_nonabsorbing, na.rm = TRUE) <= 1) %>%
  ungroup()


#### Outcomes & CS Options ####

outcome_vars <- coutcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "afd", "current_incumbent")

# Drop AfD in Outcomes
outcome_vars_wo_afd <- outcome_vars <- c("turnout", "cdu_csu", "spd", "gruene", "far_right", "current_incumbent")


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


###################### Without and Only Bavaria ######################
# Without Bavaria
did_wo_bav <- run_csdid_pipeline(df_wo_bav, outcome_vars, label = "Without Bavaria", formula = ~ east_ger + pop_density)

# Bavaria only
did_bav <- run_csdid_pipeline(df_bav, outcome_vars, label = "Bavaria Only", formula = ~ pop_density)

# Bavaria post 10H rule
did_bav_10h <- run_csdid_pipeline(df_bav_10h, outcome_vars, label = "Bavaria Post-10H", formula = ~ pop_density)


###################### East vs. West Germany (wo Bavaria) ######################
did_east <- run_csdid_pipeline(df_east, outcome_vars, label = "East Germany", formula = ~ pop_density)
did_west <- run_csdid_pipeline(df_west, outcome_vars, label = "West Germany", formula = ~ pop_density)


###################### Early vs. Late Germany ######################
# Polarisation Hypothesis
did_early_pol <- run_csdid_pipeline(df_early_p, outcome_vars_wo_afd, label = "Early (Polarisation)")
did_late_pol  <- run_csdid_pipeline(df_late_p, outcome_vars_wo_afd, label = "Late (Polarisation)")

# Tax Hypothesis
did_early_tax <- run_csdid_pipeline(df_early_t, outcome_vars_wo_afd, label = "Early - Tax Split")
did_late_tax  <- run_csdid_pipeline(df_late_t, outcome_vars_wo_afd, label = "Late - Tax Split")

# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_honest_smooth.png", plot = plot_honest_smooth, width = 12, height = 8, dpi = 300)
# ggsave(filename = "_results/hyp1_main_analyses/_figures/plot_honest_rm.png", plot = plot_honest_rm, width = 12, height = 8, dpi = 300)


###################### Interaction: Late x West Germany ######################
did_west_late <- run_csdid_pipeline(df_west_late, outcome_vars, label = "Interaction: Late x West Germany", formula = ~ pop_density)


############## High vs. Low Density ##############
# Hypothesis: Only in municipalities with high population density do new wind turbines lead to conflict
did_hd <- run_csdid_pipeline(df_hd, outcome_vars, label = "High Density")
did_ld <- run_csdid_pipeline(df_ld, outcome_vars, label = "Low Density")


############## High vs. Low AfD ##############
# Hypothesis: In municipalities with many afd supporters, effect of wind turbines leads to more backlash
did_afd_h <- run_csdid_pipeline(df_afd_h, outcome_vars, label = "High AfD")
did_afd_l <- run_csdid_pipeline(df_afd_l, outcome_vars, label = "Low AfD")


############## High vs. Low Gruene ##############
# Hypothesis: In municipalities with many green supporters, effect of wind turbines leads to less backlash
did_gruene_h <- run_csdid_pipeline(df_gruene_h, outcome_vars, label = "High Greens")
did_gruene_l <- run_csdid_pipeline(df_gruene_l, outcome_vars, label = "Low Greens")


############## Distance Variation ##############

# Estimation loop
did_dist <- map(distances, function(dist) {
  run_csdid_pipeline(
    data = df_distance, 
    outcomes = outcome_vars, 
    label = paste("Distance:", dist), 
    options = att_options_base, 
    gname = paste0("seq_group_", dist)
  )
}) %>% set_names(distances)


# Plot: one grid per distance
walk(distances, function(dist) {
  plot_list <- map(outcome_vars, function(var) {
    ggdid(did_dist[[dist]][[var]]$es) +
      ggtitle(paste0(var, " (", dist, ")")) +
      theme_minimal()
  })
  grid.arrange(grobs = plot_list, ncol = 3, top = paste("Distance:", dist))
})


############## Intensity Variation ##############

# Summary table
intensity_summary <- map_dfr(c(thresholds), function(thresh) {
  gname <- paste0("treat_", thresh)
  df_thresh <- df_intensity %>%
    filter(!is.na(.data[[gname]])) %>%
    group_by(ags) %>%
    summarise(g = first(.data[[gname]]), .groups = "drop") %>%
    summarise(
      Threshold = paste0(">= ", thresh),
      Units     = n(),
      Treated   = sum(g > 0),
      Control   = sum(g == 0)
    )
})
intensity_summary


# Estimation loop
did_intens <- map(thresholds, function(thresh) {
  gname <- paste0("treat_", thresh)
  df_thresh <- df_intensity %>% filter(!is.na(.data[[gname]]))
  
  run_csdid_pipeline(
    data = df_thresh, 
    outcomes = outcome_vars, 
    label = paste("Intensity Threshold >=", thresh), 
    options = att_options_base,
    gname = gname
  )
}) %>% set_names(paste0("treat_", thresholds))


# Plot: one grid per intensity
walk(paste0("treat_", thresholds), function(thresh) {
  plot_list <- map(outcome_vars, function(var) {
    ggdid(did_intens[[thresh]][[var]]$es) +
      ggtitle(paste0(var, " (", thresh, ")")) +
      theme_minimal()
  })
  grid.arrange(grobs = plot_list, ncol = 3, top = paste("Intensity:", thresh))
})


############## Units treated once ##############
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
did_treat_once <- run_csdid_pipeline(df_treat_once, outcome_vars, label = "One-time Treatment")


###################### Identifying Assumptions and Alterntive Outcomes ######################
# Does treatment affect net migration or agricultural land size? If yes --> Signs for sorting/negative effects of treatment
# Does it affect education or share of female/foreign population? If yes --> Signs for sorting
outcomes_alt <- c("net_migration", "agri_land", "educ", "share_fem", "share_foreign")
did_alt <- run_csdid_pipeline(df_did_ready, outcomes_alt, label = "Test Assumptions", formula = ~ pop_density + east_ger)

# How are socioeconomic indicators affected?
# Could give idea on mechanisms --> If positive effect on tax/employment etc., positive perception of WP
outcomes_se <- c("inc_tax", "busi_tax", "tax_rev", "commute_balance", "dens_work", "purch_pow", "tourism", "hinc", "share_unemp", "share_emp")
did_se <- run_csdid_pipeline(df_did_ready, outcomes_se, label = "Socioeconomic Indicators", formula = ~ pop_density + east_ger)

# Are there east - west differences?
did_see <- run_csdid_pipeline(df_east, outcomes_se, label = "Socioeconomic Indicators- East", formula = ~ pop_density)
did_sew <- run_csdid_pipeline(df_west, outcomes_se, label = "Socioeconomic Indicators- West", formula = ~ pop_density)
# Yes, there are!
# No effects on hinc AND voting in East Germnany but negative effects on both in West Germany

# Control for lagged vote share
# To reduce variance in outcome --> Good control according to Cinelli et al.
df_did_ready <- df_did_ready %>% arrange(ags, election_year) %>% group_by(ags) %>% 
  mutate(lag_far_right = lag(far_right, 1), lag_cdu_csu = lag(cdu_csu, 1), lag_spd = lag(spd, 1),
         lag_gruene = lag(gruene, 1), lag_fdp = lag(fdp, 1), lag_linke_pds = lag(linke_pds, 1), lag_turnout = lag(turnout, 1))
did_lag <- run_csdid_pipeline(df_did_ready, outcome_vars, label = "Controlling for Lag", 
                              formula = ~ pop_density + east_ger + lag_far_right + lag_cdu_csu + lag_spd + lag_gruene
                              + lag_fdp + lag_linke_pds + lag_turnout)


############## Counterfactual Estimator instead of CS-DiD ##############

# Control variables
control_vars <- " ~ treat_absorbing + pop_density"

#### fect estimation function ####
run_fect_abs <- function(outcome, controls, data) {
  message(paste("Running fect (absorbing) for:", outcome))
  formula <- as.formula(paste0(outcome, controls))
  fit <- fect(
    formula,
    data    = data,
    index   = c("ags", "seq_time"),
    method  = "ife",
    force   = "two-way",
    se      = TRUE,
    nboots  = 1000,
    min.T0  = 1
  )
  return(fit)
}


#### Run for all outcomes ####
fect_abs_results <- map(outcome_vars, ~run_fect_abs(.x, control_vars, df_did_ready)) %>% set_names(outcome_vars)

#### Plot results ####
plot_list <- map(outcome_vars, function(var) {plot(fect_abs_results[[var]], main = var)})
grid.arrange(grobs = plot_list, ncol = 3, top = textGrob("FECT Absorbing DiD - All Outcomes"))



############## De-meaning: Unit-Specific Linear Pre-Treatment Trends ##############

detrend_outcome <- function(df, outcome_var, id_var = "ags", 
                            time_var = "seq_time", group_var = "seq_group") {
  
  df %>%
    group_by(across(all_of(id_var))) %>%
    mutate(
      # Pre-treatment = periods before unit's treatment cohort (never-treated: all periods)
      is_pre = (.data[[time_var]] < .data[[group_var]]) | (.data[[group_var]] == 0),
      
      # Estimate unit-specific linear trend on pre-treatment periods only
      trend_coef = {
        pre_data <- cur_data()[cur_data()$is_pre & !is.na(cur_data()[[outcome_var]]), ]
        if (nrow(pre_data) >= 3) {  # need minimum obs for reliable trend
          coef(lm(reformulate(time_var, response = outcome_var), data = pre_data))[[time_var]]
        } else {
          NA_real_
        }
      },
      
      # Center time within unit (avoids large intercept; trend removal is origin-invariant)
      time_centered = .data[[time_var]] - min(.data[[time_var]][is_pre], na.rm = TRUE),
      
      # Subtract extrapolated trend from ALL periods (pre + post)
      "{outcome_var}_detrended" := .data[[outcome_var]] - trend_coef * time_centered
      
    ) %>%
    ungroup() %>%
    select(-is_pre, -trend_coef, -time_centered)
}


#### Apply detrending to all outcomes, add detrended cols to df ####

detrended_vars <- paste0(outcome_vars, "_detrended")

df_did_detrended <- reduce(outcome_vars, function(df, var) {
  detrend_outcome(df, outcome_var = var)
}, .init = df_did_ready)

# Quick sanity check: how many units had too few pre-periods for trend estimation?
map(detrended_vars, function(var) {
  n_na <- sum(is.na(df_did_detrended[[var]]))
  n_total <- nrow(df_did_detrended)
  message(glue::glue("{var}: {n_na} NAs ({round(n_na/n_total*100,1)}%)"))
})


#### Estimate CS-DiD on detrended outcomes ####

detrended_vars <- c("turnout_detrended", "cdu_detrended", "csu_detrended", "spd_detrended", "fdp_detrended",
                    "linke_pds_detrended", "gruene_detrended", "afd_detrended", "current_incumbent_detrended")

did_results_detrended <- map(detrended_vars, function(outcome_var) {
  message(paste("Running detrended CS-DiD for:", outcome_var))
  
  atts <- do.call(att_gt, c(
    list(yname = outcome_var),
    att_options_base           # reuse your existing options but change data input
  ))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(detrended_vars)




#### Create event study plots for all outcomes ####
plot_list <- map(detrended_vars, function(var) {
  p <- ggdid(did_results_detrended[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### Compare pre-trends: original vs detrended (key diagnostic) ####

# For a single outcome, overlay both event studies
compare_pretrends <- function(var) {
  
  orig     <- did_results[[var]]$es
  detrend  <- did_results_detrended[[paste0(var, "_detrended")]]$es
  
  bind_rows(
    tibble(
      e       = orig$egt,
      att     = orig$att.egt,
      se      = orig$se.egt,
      spec    = "Original"
    ),
    tibble(
      e       = detrend$egt,
      att     = detrend$att.egt,
      se      = detrend$se.egt,
      spec    = "Detrended"
    )
  ) %>%
    mutate(
      ci_lo = att - 1.96 * se,
      ci_hi = att + 1.96 * se
    ) %>%
    ggplot(aes(x = e, y = att, color = spec, fill = spec)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c("Original" = "#2C3E50", "Detrended" = "#E74C3C")) +
    scale_fill_manual(values  = c("Original" = "#2C3E50", "Detrended" = "#E74C3C")) +
    labs(
      title    = glue::glue("Event Study: {var}"),
      subtitle = "Original vs. unit-trend detrended outcome",
      x        = "Periods relative to treatment",
      y        = "ATT estimate",
      color    = NULL, fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Plot for all outcomes
comparison_plots <- map(outcome_vars, compare_pretrends)
grid.arrange(grobs = comparison_plots, ncol = 3)


