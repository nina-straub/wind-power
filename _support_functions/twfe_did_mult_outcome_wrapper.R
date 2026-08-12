# Load packages
library(purrr)
library(fixest)
library(bacondecomp)
library(did)
library(DIDmultiplegtDYN)
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


###################### Main dCDH Pipline Estimation Function ######################

run_dcdh_pipeline <- function(data, 
                              outcomes, 
                              label = "", 
                              treatment = "cum_wind_count_3km",
                              group = "ags",
                              time = "seq_time",
                              effects = 3,
                              placebo = 3,
                              controls = NULL,
                              trends_nonparam = NULL,
                              cluster = "ags",
                              normalized = TRUE,
                              same_switchers = FALSE,
                              same_switchers_pl = FALSE,
                              by = NULL,
                              options = list()) {
  
  # Helper to resolve scalar vs outcome-specific arguments
  get_val <- function(param, outcome) {
    if (!is.null(names(param)) && outcome %in% names(param)) {
      return(param[[outcome]])
    }
    if (is.numeric(param) && length(param) == 1 && is.null(names(param))) {
      return(param)
    }
    if ("default" %in% names(param)) return(param[["default"]])
    return(param[[1]])
  }
  
  # 1. Run did_multiplegt_dyn across outcomes
  results <- map(outcomes, function(outcome) {
    
    eff_i <- get_val(effects, outcome)
    plc_i <- get_val(placebo, outcome)
    
    message(paste0("Running dCDH DiD for: ", outcome, 
                   " (effects = ", eff_i, ", placebo = ", plc_i, ")",
                   ifelse(label != "", paste("|", label), "")))
    
    # Define argument list
    base_args <- list(
      df                = data,
      outcome           = outcome,
      group             = group,
      time              = time,
      treatment         = treatment,
      effects           = eff_i,
      placebo           = plc_i,
      controls          = controls,
      trends_nonparam   = trends_nonparam,
      cluster           = cluster,
      normalized        = normalized,
      same_switchers    = same_switchers,
      same_switchers_pl = same_switchers_pl,
      by = by,
      graph_off         = TRUE
    )
    
    args_list <- modifyList(base_args, options)
    
    # Execute call
    res <- do.call(did_multiplegt_dyn, args_list)
    return(res)
  }) %>% set_names(outcomes)
  
  # 2. Extract and format plots
  plot_list <- map(outcomes, function(var) {
    results[[var]]$plot +
      ggtitle(var) +
      theme_minimal()
  })
  
  # 3. Render grid layout
  grid.arrange(
    grobs = plot_list, 
    ncol  = min(3, length(outcomes)), 
    top   = if (label != "") label else NULL
  )
  
  return(results)
}
