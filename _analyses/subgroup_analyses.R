###################### Subgroup Analysis of Panel Data ######################

###################### Set the scene  ######################

#### Data ####
# df_did_ready from main analysis script

#### For subgroup analysis ####

# Without and Only Bavaria
df_wo_bav <- df_did_ready %>% filter(!str_starts(as.character(ags), "9"))
df_bav <- df_did_ready %>% filter(str_starts(as.character(ags), "9"))
df_bav_10h <- df_bav %>% filter(election_year > 2014)

# East vs. West Germany (wo Bavaria)
df_east <- df_did_ready %>% filter(east_ger == 1)
df_west <- df_wo_bav %>% filter(east_ger == 0)

# Early vs. Late Germany (polarisation hypothesis)
df_early <- df_did_ready %>% filter(election_year < 2015)
df_late <- df_did_ready %>% filter(election_year > 2015)

# Early vs. Late Germany (tax hypothesis)
df_early <- df_did_ready %>% filter(election_year < 2009)
df_late <- df_did_ready %>% filter(election_year > 2009)

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
distances <- c("0km", "3km", "5km", "10km")
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


#### Outcomes & CS Options ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")

# Drop AfD in Outcomes
outcome_vars_wo_afd <- c("turnout", "cdu_csu", "spd", "fdp", 
                         "linke_pds", "gruene", "current_incumbent", "far_right")


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


###################### Main CS Pipline Estimation Function ######################

run_csdid_pipeline <- function(data, outcomes, label = "", formula = ~ pop_density + east_ger, options = att_options_base, gname = "seq_group") {
  
  # Override default group name in options if a custom one is supplied
  options$gname <- gname
  
  # Run CS-DiD models across all outcomes
  results <- map(outcomes, function(outcome) {
    message(paste("Running CS-DiD for:", outcome, ifelse(label != "", paste("|", label), "")))
    
    # Merge core arguments with base list options
    args_list <- c(list(yname = outcome, data = data, xformla = formula), options)
    atts <- do.call(att_gt, args_list)
    es   <- aggte(atts, type = "dynamic", na.rm = TRUE)
    
    return(list(atts = atts, es = es))
  }) %>% set_names(outcomes)
  
  # Generate plots
  plot_list <- map(outcomes, function(var) {
    ggdid(results[[var]]$es) +
      ggtitle(paste0(var, " (", label, ")")) +
      theme_minimal()
  })
  
  grid.arrange(grobs = plot_list, ncol = 3, top = if(label != "") label else NULL)
  
  # Return results list
  return(results)
}



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
did_early_pol <- run_csdid_pipeline(df_early, outcome_vars_wo_afd, label = "Early (Polarisation)")
did_late_pol  <- run_csdid_pipeline(df_late, outcome_vars_wo_afd, label = "Late (Polarisation)")

# Tax Hypothesis
did_early_tax <- run_csdid_pipeline(df_early, outcome_vars_wo_afd, label = "Early (Tax)")
did_late_tax  <- run_csdid_pipeline(df_late, outcome_vars, label = "Late (Tax)")


###################### Interaction: Late x West Germany ######################
did_west_late <- run_csdid_pipeline(df_west_late, outcome_vars, label = "Interaction: Late x West Germany")


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


