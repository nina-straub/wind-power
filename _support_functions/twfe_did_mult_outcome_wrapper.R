# Load packages
library(purrr)
library(dplyr)
library(stringr)
library(fixest)
library(bacondecomp)
library(did)
library(polars)
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
    es   <- aggte(atts, type = "dynamic", na.rm = TRUE, min_e = -4, max_e = 5)
    
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


###################### Extraction of CS Results for Plotting ######################

extract_cs_df <- function(did_result, label) {
  map_dfr(names(did_result), function(var) {
    es <- did_result[[var]]$es
    crit_val <- es$crit.val.egt
    
    data.frame(
      outcome   = var,
      e         = es$egt,
      estimate  = es$att.egt,
      se        = es$se.egt,
      estimator = label
    ) %>%
      mutate(
        conf.low  = estimate - (crit_val * se),
        conf.high = estimate + (crit_val * se)
      )
  })
}


###################### Main dCDH Pipline Estimation Function ######################

run_dcdh_pipeline <- function(
    data, 
    outcomes, 
    label = "", 
    options = dcdh_options_base,
    effects_default = 6,
    placebo_default = 3,
    effects_afd = 3,
    placebo_afd = 1
) {
  map(outcomes, function(var) {
    message(paste("Running dCDH for:", var, ifelse(label != "", paste("|", label), "")))
    
    # Dynamic window override
    plc_val <- if (var == "afd") placebo_afd else placebo_default
    eff_val <- if (var == "afd") effects_afd else effects_default
    
    # Build argument list
    args_list <- c(
      list(
        df = data,
        outcome = var,
        effects = eff_val,
        placebo = plc_val
      ),
      options
    )
    
    res_dcdh <- do.call(did_multiplegt_dyn, args_list)
    
    extract_dcdh_results(res_dcdh, var, label)
  }) %>% set_names(outcomes)
}



###################### Extraction function to apply index shift for dCDH (e = x - 1) ######################

extract_dcdh_results <- function(dcdh_res, var_name, label) {
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
      estimator = label,
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
    est_cols   = c("CS (did)" = "#2b5c8f", "dCDH (did_multiplegt_dyn)" = "#fd7107",
                   "Winner" = "#2b5c8f", "Loser" = "#fd7107",
                   "East" = "#2b5c8f", "West" = "#fd7107",
                   "Early" = "#2b5c8f", "Late" = "#fd7107",
                   "Treated Once" = "#2b5c8f"),
    est_shapes = c("CS (did)" = 16,        "dCDH (did_multiplegt_dyn)" = 17,
                   "Winner" = 16, "Loser" = 17,
                   "East" = 16, "West" = 17,
                   "Early" = 16, "Late" = 17,
                   "Treated Once" = 16),
    est_lines  = c("CS (did)" = "solid",   "dCDH (did_multiplegt_dyn)" = "solid",
                   "Winner" = "solid", "Loser" = "solid",
                   "East" = "solid", "West" = "solid",
                   "Early" = "solid", "Late" = "solid",
                   "Treated Once" = "solid"),
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

