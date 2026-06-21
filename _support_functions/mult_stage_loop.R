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


# Function to generate plots for each outcome and stores them in a list
generate_sequential_plots <- function(results_list) {
  
  stage_colors <- c("#2166AC", "#D6604D", "#4DAC26", "#E69F00", "#999999")
  
  # Iterate over data and outcomes
  plots <- purrr::imap(results_list, function(df, outcome_name) {
    
    ggplot(df, aes(x = event_time, y = att, color = stage, fill = stage, group = stage)) +
      geom_hline(yintercept = 0,    linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.0) +
      scale_color_manual(values = stage_colors) +
      scale_fill_manual(values  = stage_colors) +
      labs(
        title = outcome_name, # Dynamically uses the sublist/outcome name
        subtitle = "Conditioned on identical prior history",
        x = "Periods relative to treatment", 
        y = "ATT", 
        color = NULL, 
        fill = NULL
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5, color = "gray30")
      )
  })
  
  # Return named list of ggplot objects
  return(plots)
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

