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

run_dcdh_pipeline <- function(df, 
                              outcome, 
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
  results <- map(outcome, function(outcome) {
    
    eff_i <- get_val(effects, outcome)
    plc_i <- get_val(placebo, outcome)
    
    message(paste0("Running dCDH DiD for: ", outcome, 
                   " (effects = ", eff_i, ", placebo = ", plc_i, ")",
                   ifelse(label != "", paste("|", label), "")))
    
    # Define argument list
    base_args <- list(
      df                = df,
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
  }) %>% set_names(outcome)
  
  # 2. Extract and format plots
  plot_list <- map(outcome, function(var) {
    results[[var]]$plot +
      ggtitle(var) +
      theme_minimal()
  })
  
  # 3. Render grid layout
  grid.arrange(
    grobs = plot_list, 
    ncol  = min(3, length(outcome)), 
    top   = if (label != "") label else NULL
  )
  
  return(results)
}



###################### Extraction function to apply index shift for dCDH (e = x - 1) ######################

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



###################### Overlay Plot ######################

create_overlay_plot <- function(
    df_all, 
    outcome_vars = unique(df_all$outcome),
    outcome_labels = c(
      "turnout"           = "Voter Turnout",
      "cdu_csu"           = "CDU/CSU",
      "spd"               = "SPD",
      "gruene"            = "Greens",
      "afd"               = "AfD",
      "current_incumbent" = "Current Incumbent",
      "far_right"         = "Far Right"
    ),
    est_cols   = c("CS (did)" = "#2b5c8f", "dCDH (did_multiplegt_dyn)" = "#fd7107", "Winner" = "#2b5c8f", "Loser" = "#fd7107"),
    est_shapes = c("CS (did)" = 16,        "dCDH (did_multiplegt_dyn)" = 17, "Winner" = 16, "Loser" = 17),
    est_lines  = c("CS (did)" = "solid",   "dCDH (did_multiplegt_dyn)" = "solid", "Winner" = "solid", "Loser" = "solid"),
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

