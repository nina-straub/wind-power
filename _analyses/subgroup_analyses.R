###################### Subgroup Analysis of Panel Data ######################

###################### Set the scene  ######################

#### Data ####
df_panel_cov <- readRDS("_data/df_panel_cov.rds")


#### Re-construct treatment variables for different packages & add additional subgroup dummies ####

# Add variables
df_did_ready <- df_panel_cov %>%
  mutate(
    # Create sequential time index (1 to 10) for did package
    seq_time = match(election_year, sort(unique(df_panel_cov$election_year))),
    # Create sequential group index (never treated == 0) for did package
    seq_group = match(treat_absorbing_cs, sort(unique(df_panel_cov$election_year))),
    seq_group = ifelse(is.na(seq_group), 0, seq_group),
    # Create time-to-treatment variable for fixest package
    time_to_treatment = ifelse(seq_group > 0, seq_time - seq_group, -1000),
    # Clean lagged wind count variable of NAs
    cat_cum_lag_wind_count_3km = ifelse(election_year == 1990 & is.na(cat_cum_lag_wind_count_3km),
                                        0, cat_cum_lag_wind_count_3km),
    # Create current incumbent election results
    current_incumbent = case_when(
      election_year %in% c(1990, 1994, 1998, 2013) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(fdp, 0),
      election_year %in% c(2002, 2005) ~ coalesce(spd, 0) + coalesce(gruene, 0),
      election_year %in% c(2009, 2017, 2021) ~ coalesce(cdu, 0) + coalesce(csu, 0) + coalesce(spd, 0),
      election_year == 2025 ~ coalesce(gruene, 0) + coalesce(spd, 0) + coalesce(fdp, 0),
      TRUE ~ NA_real_),
    # Create "other parties" election results
    others = 1 - (coalesce(cdu, 0) + coalesce(afd, 0) + coalesce(csu, 0)
                  + coalesce(spd, 0) + coalesce(fdp, 0) + coalesce(linke_pds, 0)
                  + coalesce(gruene, 0)),
    # Add east Germany dummy
    east_ger = ifelse(substr(ags, 1, 2) %in% c("12", "13", "14", "15", "16"), 1, 0)) %>%
  # Filter for certain states
  # filter(str_starts(ags, "12")) %>%
  # Filter for certain time period
  # filter(!election_year < 1990)
  # Create dose variable for contdid package
  group_by(ags) %>%
  mutate(treat_dose = if (any(seq_group > 0)) {
    max(wind_count_3km[seq_time == seq_group], na.rm = TRUE)
  } else {
    wind_count_3km
  }) %>%
  ungroup() %>%
  # Ensure ID variable is numeric for the did package
  mutate(ags = as.numeric(ags)) %>%
  as.data.frame()


###################### Outcomes ######################

outcome_vars <- c("turnout", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")


#### Set CS options ####
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

df_wo_bav <- df_did_ready %>% filter(!str_starts(as.character(ags), "9"))
df_bav <- df_did_ready %>% filter(str_starts(as.character(ags), "9"))


#### CS-DiD without Bavaria ####
did_wo_bav <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_wo_bav, 
                                 xformla = ~ east_ger + pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_wo_bav[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD Bavaria only ####
did_bav <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_bav, 
                                 xformla = ~ pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_bav[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


# Bavaria post 10H rule
df_bav_10h <- df_bav %>% filter(election_year > 2014)

did_bav_10h <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_bav_10h, 
                                 xformla = ~ pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_bav_10h[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



###################### East vs. West Germany (wo Bavaria) ######################

df_east <- df_did_ready %>% filter(east_ger == 1)
df_west <- df_wo_bav %>% filter(east_ger == 0)


#### CS-DiD East Germany ####
did_east <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_east, 
                                 xformla = ~ pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_east[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD West Germany ####
did_west <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_west, 
                                 xformla = ~ pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_west[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



###################### Early vs. Late Germany (Polarisation Definition) ######################

df_early <- df_did_ready %>% filter(election_year < 2015)
df_late <- df_did_ready %>% filter(election_year > 2015)


#### CS-DiD  Early ####

# Drop AfD in Outcomes
outcome_vars_wo_afd <- c("turnout", "cdu", "csu", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "current_incumbent", "far_right")

did_early <- map(outcome_vars_wo_afd, function(outcome_vars_wo_afd) {
  message(paste("Running CS-DiD for:", outcome_vars_wo_afd))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars_wo_afd, data = df_early, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars_wo_afd)

# Event study plot
plot_list <- map(outcome_vars_wo_afd, function(var) {
  p <- ggdid(did_early[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD Late ####
did_late <- map(outcome_vars_wo_afd, function(outcome_vars_wo_afd) {
  message(paste("Running CS-DiD for:", outcome_vars_wo_afd))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars_wo_afd, data = df_late, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars_wo_afd)

# Event study plot
plot_list <- map(outcome_vars_wo_afd, function(var) {
  p <- ggdid(did_late[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



###################### Early vs. Late Germany (Tax Definition) ######################

df_early <- df_did_ready %>% filter(election_year < 2009)
df_late <- df_did_ready %>% filter(election_year > 2009)


#### CS-DiD  Early ####

# Drop AfD in Outcomes
outcome_vars_wo_afd <- c("turnout", "cdu", "csu", "cdu_csu", "spd", "fdp", 
                         "linke_pds", "gruene", "current_incumbent")

did_early <- map(outcome_vars_wo_afd, function(outcome_vars_wo_afd) {
  message(paste("Running CS-DiD for:", outcome_vars_wo_afd))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars_wo_afd, data = df_early, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars_wo_afd)

# Event study plot
plot_list <- map(outcome_vars_wo_afd, function(var) {
  p <- ggdid(did_early[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD Late ####
did_late <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_late, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_late[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



###################### Interaction: Late x West Germany ######################

df_west_late <- df_wo_bav %>% filter(east_ger == 0, election_year > 2015)

#### CS-DiD West Germany ####
did_west_late <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_west_late, 
                                 xformla = ~ pop_density), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_west_late[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



############## High vs. Low Density ##############
# Hypothesis: Only in municipalities with high population density do new wind turbines lead to conflict

df_hd <- df_did_ready %>% filter(pop_density > quantile(pop_density, 0.5, na.rm = TRUE),
                                 pop_density <= quantile(pop_density, 0.9, na.rm = TRUE))
df_ld <- df_did_ready %>% filter(pop_density < quantile(pop_density, 0.5, na.rm = TRUE))

#### CS-DiD High Density ####
did_hd <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_hd, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_hd[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)

#### CS-DiD Low Density ####
did_ld <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_ld, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_ld[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



############## High vs. Low AfD ##############
# Hypothesis: In municipalities with many afd supporters, effect of wind turbines leads to more backlash

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


#### CS-DiD High AfD ####
did_afd_h <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_afd_h, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_afd_h[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD Low AfD ####
did_afd_l <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_afd_l, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_afd_l[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



############## High vs. Low Gruene ##############
# Hypothesis: In municipalities with many green supporters, effect of wind turbines leads to less backlash

# Filter for municipalities with above-average AfD vote share
df_gruene_h <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_gruene = mean(gruene, na.rm = TRUE)) %>%
                     filter(mean_gruene > median(mean_gruene, na.rm = TRUE)) %>%
                     pull(ags)))

# Filter for municipalities with below-average AfD vote share
df_gruene_l <- df_did_ready %>%
  filter(ags %in% (df_did_ready %>%
                     group_by(ags) %>%
                     summarise(mean_gruene = mean(gruene, na.rm = TRUE)) %>%
                     filter(mean_gruene < median(mean_gruene, na.rm = TRUE)) %>%
                     pull(ags)))


#### CS-DiD High Greens ####
did_gruene_h <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_gruene_h, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_gruene_h[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)


#### CS-DiD Low Greens ####
did_gruene_l <- map(outcome_vars, function(outcome_vars) {
  message(paste("Running CS-DiD for:", outcome_vars))
  
  atts <- do.call(att_gt, c(list(yname = outcome_vars, data = df_gruene_l, 
                                 xformla = ~ pop_density + east_ger), att_options_base))
  
  es <- aggte(atts, type = "dynamic", na.rm = TRUE)
  return(list(atts = atts, es = es))
}) %>% set_names(outcome_vars)

# Event study plot
plot_list <- map(outcome_vars, function(var) {
  p <- ggdid(did_gruene_l[[var]]$es) +
    ggtitle(var) +
    theme_minimal()
  return(p)
})
grid.arrange(grobs = plot_list, ncol = 3)



############## Distance Variation ##############

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

  
#### Set CS options ####
# Set vars for conditional parallel trends with xformla
att_opt_2 <- list(
  tname = "seq_time",
  idname = "ags",
  panel = TRUE, 
  allow_unbalanced_panel = TRUE,
  clustervars = "ags",
  control_group = "nevertreated",
  xformla = ~ pop_density + east_ger,
  anticipation = 0,
  bstrap = TRUE,
  biters = 1000
)

#### CS-DiD Different Distances ####
did_dist <- map(distances, function(dist) {
  gname <- paste0("seq_group_", dist)
  
  map(outcome_vars, function(outcome) {
    message(paste("Running CS-DiD for:", outcome, "| Distance:", dist))
    
    atts <- do.call(att_gt, c(list(yname = outcome, gname = gname, data = df_distance), att_opt_2))
    es   <- aggte(atts, type = "dynamic", na.rm = TRUE)
    return(list(atts = atts, es = es))
  }) %>% set_names(outcome_vars)
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


#### CS-DiD Different Intensities ####
# Use same cs_options as in distance variation analysis
did_intens <- map(thresholds, function(thresh) {
  gname <- paste0("treat_", thresh)
  
  df_thresh <- df_intensity %>% filter(!is.na(.data[[gname]]))
  
  map(outcome_vars, function(outcome) {
    message(paste("Running CS-DiD for:", outcome, "| Intensity:", thresh))
    atts <- do.call(att_gt, c(list(yname = outcome, gname = gname, data = df_thresh), att_opt_2))
    es   <- aggte(atts, type = "dynamic", na.rm = TRUE)
    return(list(atts = atts, es = es))
  }) %>% set_names(outcome_vars)
}) %>% set_names(paste0("treat_", thresholds))


# Plot
walk(paste0("treat_", thresholds), function(thresh) {
  plot_list <- map(outcome_vars, function(var) {
    ggdid(did_intens[[thresh]][[var]]$es) +
      ggtitle(paste0(var, " (", thresh, ")")) +
      theme_minimal()
  })
  grid.arrange(grobs = plot_list, ncol = 3, top = paste("Intensity:", thresh))
})











