# Load packages
library(fixest)
library(did)
library(ggplot2)
library(gridExtra)

###################### Function wrapper for TWFE event study estimation ######################
run_twfe_event_study <- function(outcome, data, fixed_effects = "ags + election_year", covariates = "") {
  # Formula construction: outcome ~ treatment + covariates | fixed effects
  formula_str <- paste0(outcome, " ~ i(time_to_treatment, ref = c(-1, -1000)) ", 
                        covariates, " | ", fixed_effects)
  # Run the TWFE model
  model <- feols(as.formula(formula_str), data = data, cluster = ~ags)
  return(list(model = model))
}

###################### Function wrapper for goodman bacon decomposition ######################
run_bacon_decomposition <- function(outcome, data, covariates) {
  # Drop NAs
  df_clean <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(treat_absorbing), !is.na(ags), !is.na(seq_time), !is.na(pop_density)
    )
  # Find no. of time periods
  expected_periods <- n_distinct(df_clean$seq_time)
  # Balance panel
  df_balanced <- df_clean %>%
    group_by(ags) %>%
    filter(n() == expected_periods) %>%  
    ungroup()
  # Check if any data left after balancing
  if (nrow(df_balanced) == 0) {
    warning(paste("No balanced data left for outcome:", outcome))
    return(NULL)
  }
  # Run decomposition
  fml <- as.formula(paste(outcome, covariates))
  bgd <- bacon(fml, 
               data = df_balanced, 
               id_var = "ags", 
               time_var = "seq_time")
  return(bgd)
}


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

