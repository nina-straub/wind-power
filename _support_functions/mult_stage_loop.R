# Load packages
library(tidyverse)
library(patchwork)
library(purrr)
library(ggplot2)


# Function for creating arbitrary number of sequential treatment "stages"
# In first stage, it compares units treated for the first time vs. never treated
# Treated units in each stage are the sample for the subsequent stage, comparing treated vs. non-treated
run_sequential_stages <- function(outcome, att_opts, stages) {
  
  get_seq_time <- function(year) {
    df_did_ready %>% filter(election_year == year) %>% pull(seq_time) %>% unique()
  }
  
  att_opts$yname <- outcome
  
  # Stage 1: derive year and seq_group from stages[[1]]
  s1_year     <- stages[[1]]$year
  s1_seq_time <- get_seq_time(s1_year)
  
  # Create stage 1 panel data
  current_panel <- df_did_ready %>%
    filter(seq_group %in% c(s1_seq_time, 0)) %>%
    group_by(ags) %>%
    mutate(
      s1_treat = as.integer(any(election_year == s1_year & treat_nonabsorbing == 1)),
      s1_group = s1_treat * s1_seq_time
    ) %>%
    ungroup()
  
  stage_plots_data <- list()
  
  # Loop over stages and estimate effect with CS DiD
  for (i in seq_along(stages)) {
    cfg    <- stages[[i]]
    t_var  <- paste0(cfg$suffix, "_treat")
    g_var  <- paste0(cfg$suffix, "_group")
    
    # Update panel data for stages > 1
    if (i > 1) {
      prev_t_var    <- paste0(stages[[i-1]]$suffix, "_treat")
      time_current  <- get_seq_time(cfg$year)
      
      current_panel <- current_panel %>%
        filter(.data[[prev_t_var]] == 1) %>%
        group_by(ags) %>%
        mutate(
          !!t_var := as.integer(any(election_year == cfg$year & treat_nonabsorbing == 1)),
          !!g_var := .data[[t_var]] * time_current
        ) %>%
        ungroup()
    }
    
    # CS estimator
    cs_res  <- do.call(att_gt, c(list(gname = g_var, data = current_panel), att_opts))
    agg_res <- aggte(cs_res, type = "dynamic", na.rm = TRUE)
    
    # Create df summary for each stage
    stage_plots_data[[i]] <- tibble(
      event_time = agg_res$egt,
      att        = agg_res$att,
      se         = agg_res$se,
      ci_lo      = agg_res$att - agg_res$crit.val.egt * agg_res$se,
      ci_hi      = agg_res$att + agg_res$crit.val.egt * agg_res$se,
      stage      = paste0("Stage ", i, " (", cfg$year, ")")
    )
  }
  
  plot_df <- bind_rows(stage_plots_data)
  
  return(plot_df)
  
}


# Generate summary table to check if pre-set stages will yield a reasonable sample size
generate_pipeline_summary <- function(df, stages, initial_filter_groups = c(4, 0)) {
  # Initialize panel with Stage 1 base filtering
  current_panel <- df %>% 
    filter(seq_group %in% initial_filter_groups)
  
  # List to hold summary row for each stage
  summary_rows <- list()
  
  for (i in seq_along(stages)) {
    cfg <- stages[[i]]
    t_var <- paste0(cfg$suffix, "_treat")
    
    # Apply sequential filtering and treatment assignment
    if (i == 1) {
      current_panel <- current_panel %>%
        group_by(ags) %>%
        mutate(!!t_var := as.integer(any(election_year == cfg$year & treat_nonabsorbing == 1))) %>%
        ungroup()
    } else {
      prev_cfg <- stages[[i-1]]
      prev_t_var <- paste0(prev_cfg$suffix, "_treat")
      
      # Filter to only include those treated in the previous stage
      current_panel <- current_panel %>%
        filter(.data[[prev_t_var]] == 1) %>%
        group_by(ags) %>%
        mutate(!!t_var := as.integer(any(election_year == cfg$year & treat_nonabsorbing == 1))) %>%
        ungroup()
    }
    
    # Calculate the summary stats for the current stage
    stage_summary <- current_panel %>%
      group_by(ags) %>%
      summarise(is_treated = max(.data[[t_var]], na.rm = TRUE), .groups = "drop") %>%
      summarise(
        Stage = paste0("Stage ", i, " (", cfg$year, ")"),
        Units = n(),
        Treated = sum(is_treated == 1),
        Control = sum(is_treated == 0)
      )
    
    summary_rows[[i]] <- stage_summary
  }
  
  # Combine all stages into a single data frame
  pipeline_summary_table <- bind_rows(summary_rows)
  return(pipeline_summary_table)
}

# Function to generate plots for each outcome and stores them in a list
generate_sequential_plots <- function(results_list) {
  
  stage_colors <- c("#2166AC", "#D6604D", "#4DAC26", "#E69F00", "#999999")
  stage_shapes <- c(16, 17, 15, 18, 8)
  
  outcome_labels <- c(
    "turnout"           = "Voter Turnout",
    "cdu"               = "CDU",
    "csu"               = "CSU",
    "cdu_csu"           = "CDU/CSU",
    "spd"               = "SPD",
    "fdp"               = "FDP",
    "linke_pds"         = "Die Linke/PDS",
    "gruene"            = "Greens",
    "afd"               = "AfD",
    "far_right"         = "Far Right",
    "current_incumbent" = "Current Incumbent"
  )
  
  n_out <- length(results_list)
  ncol  <- 2
  
  # Iterate over data and outcomes
  plots <- purrr::imap(results_list, function(df, outcome_name) {
    
    idx <- which(names(results_list) == outcome_name)
    
    # Axis titles only on the left column / bottom row, mimicking global labels
    y_lab <- if ((idx - 1) %% ncol == 0) "Estimate"   else NULL
    x_lab <- if (idx > n_out - ncol)     "Event Time" else NULL
    
    # Keep stage order as produced by run_sequential_stages()
    df <- dplyr::mutate(df, stage = factor(stage, levels = unique(df$stage)))
    n_stages <- nlevels(df$stage)
    
    ggplot(df, aes(x = event_time, y = att,
                   color = stage, fill = stage, shape = stage, group = stage)) +
      geom_hline(yintercept = 0,    linetype = "dashed", color = "gray70", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray70", linewidth = 0.4) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 2.0) +
      scale_x_continuous(breaks = seq(-6, 6, 1)) +
      scale_y_continuous(labels = function(x) format(x, scientific = FALSE, trim = TRUE)) +
      scale_color_manual(values = stage_colors[seq_len(n_stages)]) +
      scale_fill_manual(values  = stage_colors[seq_len(n_stages)]) +
      scale_shape_manual(values = stage_shapes[seq_len(n_stages)]) +
      labs(
        title    = unname(outcome_labels[outcome_name]) %||% outcome_name,
        x        = x_lab,
        y        = y_lab,
        color    = NULL,
        fill     = NULL,
        shape    = NULL
      ) +
      guides(color = guide_legend(nrow = 1),
             fill  = guide_legend(nrow = 1),
             shape = guide_legend(nrow = 1)) +
      theme_classic() +
      theme(
        plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.text        = element_text(size = 13, color = "black"),
        axis.title       = element_text(face = "bold", size = 14),
        legend.title     = element_blank(),
        legend.text      = element_text(size = 13, face = "bold"),
        legend.key.width = unit(1.0, "cm")
      )
  })
  
  # Return named list of ggplot objects
  return(plots)
}
