# Load packages
library(HonestDiD)
library(gridExtra)
library(cowplot)
library(purrr)
library(gridExtra)


# Function for HonestDiD Smoothness Sensitivity Test for Multiple Outcomes
run_honest_smoothness <- function(var_name, results_list, ...) {
  # Extract event study object from CS-DiD results
  es_obj <- results_list[[var_name]]$es
  # Run smoothness sensitivity analysis for e = 0
  sens_smooth <- honest_did(
    es = es_obj,
    e = 1,
    type = "smoothness",
    ...
  )
  return(sens_smooth)
}


# Function for HonestDiD Relative Magnitudes Sensitivity for Multiple Outcomes
run_honest_rm <- function(var_name, results_list, mbar_seq = seq(0, 0.5, by = 0.05)) {
  # Extract event study object from CS-DiD results
  es_obj <- results_list[[var_name]]$es
  # Run relative magnitude sensitivity analysis for e = 0
  sens_rm <- honest_did(
    es = es_obj,
    e = 1,
    type = "relative_magnitude",
    Mbarvec = mbar_seq
  )
  return(sens_rm)
}


# Function to generate a grid of sensitivity plots for multiple outcomes
generate_sensitivity_grid_plot <- function(
    outcome_vars, 
    results_list, 
    type = c("smooth", "rm"), 
    ncol = 2,
    outcome_labels = c(
      "turnout"           = "Voter Turnout",
      "cdu"               = "CDU",
      "csu"               = "CSU",
      "cdu_csu"           = "CDU/CSU",
      "spd"               = "SPD",
      "gruene"            = "Greens",
      "afd"               = "AfD",
      "fdp"               = "FDP",
      "linke_pds"         = "Die Linke/PDS",
      "current_incumbent" = "Current Incumbent"
    ),
    method_cols = c("Original" = "#2b5c8f", "C-LF" = "#fd7107", "FLCI" = "#fd7107", "C-BF" = "#e7298a")
) {
  type <- match.arg(type)
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Safeguard: filter out variables not present in results_list
  missing_vars <- setdiff(outcome_vars, names(results_list))
  if (length(missing_vars) > 0) {
    warning("Skipping outcomes not found in results_list: ", paste(missing_vars, collapse = ", "))
    outcome_vars <- intersect(outcome_vars, names(results_list))
  }
  
  # Build individual subplots
  plots <- map(outcome_vars, function(var_name) {
    sens_data <- results_list[[var_name]]
    
    if (type == "smooth") {
      p <- HonestDiD::createSensitivityPlot(sens_data$robust_ci, sens_data$orig_ci)
    } else {
      p <- HonestDiD::createSensitivityPlot_relativeMagnitudes(sens_data$robust_ci, sens_data$orig_ci)
    }
    
    pretty_title <- outcome_labels[[var_name]] %||% var_name
    
    # Suppress redundant color scale warnings while applying custom styling
    p <- suppressWarnings(
      p +
        scale_color_manual(values = method_cols, na.value = "#fd7107") +
        scale_y_continuous(labels = function(x) format(x, scientific = FALSE, trim = TRUE)) +
        labs(title = pretty_title) +
        theme_classic() +
        theme(
          plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
          axis.text        = element_text(size = 13, color = "black"),
          axis.title       = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.text      = element_text(size = 14, face = "bold"),
          legend.key.width = unit(1.4, "cm")
        )
    )
    
    return(p)
  })
  
  # Shared legend extraction
  shared_legend <- cowplot::get_legend(plots[[1]])
  plots_no_legend <- map(plots, ~ .x + theme(legend.position = "none"))
  
  # Grid assembly
  plots_grid <- cowplot::plot_grid(plotlist = plots_no_legend, ncol = ncol)
  grid_title <- ""
  
  cowplot::ggdraw() +
    cowplot::draw_label(grid_title, x = 0.5, y = 0.98, fontface = "bold", size = 18) +
    cowplot::draw_plot(plots_grid, x = 0, y = 0.07, width = 1, height = 0.89) +
    cowplot::draw_plot(shared_legend, x = 0, y = 0.01, width = 1, height = 0.05)
}


