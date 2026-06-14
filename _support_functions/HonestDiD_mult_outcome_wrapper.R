# Load packages
library(HonestDiD)
library(gridExtra)
library(cowplot)
library(purrr)
library(gridExtra)


# Function for HonestDiD Smoothness Sensitivity Test for Multiple Outcomes
run_honest_smoothness <- function(var_name, results_list) {
  # Extract event study object from CS-DiD results
  es_obj <- results_list[[var_name]]$es
  # Run smoothness sensitivity analysis for e = 0
  sens_smooth <- honest_did(
    es = es_obj,
    e = 0,
    type = "smoothness"
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
    e = 0,
    type = "relative_magnitude",
    Mbarvec = mbar_seq
  )
  return(sens_rm)
}


# Function to generate a grid of sensitivity plots for multiple outcomes
generate_sensitivity_grid_plot <- function(outcome_vars, results_list, type = c("smooth", "rm"), ncol = 3) {
  type <- match.arg(type)
  # Generate individual plots
  plots <- map(outcome_vars, function(var_name) {
    sens_data <- results_list[[var_name]]
    if (type == "smooth") {
      p <- createSensitivityPlot(sens_data$robust_ci, sens_data$orig_ci) +
        labs(title = paste(var_name))
    } else if (type == "rm") {
      p <- createSensitivityPlot_relativeMagnitudes(sens_data$robust_ci, sens_data$orig_ci) +
        labs(title = paste(var_name))
    }
    # Apply minimal theme and hide individual legends
    p <- p + theme_minimal() + theme(legend.position = "none")
    
    return(p)
  })
  # Re-create one plot temporarily to extract legend
  sample_data <- results_list[[outcome_vars[1]]]
  if (type == "smooth") {
    legend_plot <- createSensitivityPlot(sample_data$robust_ci, sample_data$orig_ci)
  } else {
    legend_plot <- createSensitivityPlot_relativeMagnitudes(sample_data$robust_ci, sample_data$orig_ci)
  }
  # Extract legend
  shared_legend <- cowplot::get_legend(
    legend_plot + 
      theme_minimal() + 
      theme(legend.position = "bottom") +
      labs(color = "Method")
  )
  # Main title
  grid_title <- if (type == "smooth") "Smoothness Sensitivity Analysis" else "Relative Magnitude Sensitivity Analysis"
  # Arrange into grid with the shared legend placed at the bottom
  grid.arrange(
    grobs = plots, 
    ncol = ncol, 
    top = grid_title,
    bottom = shared_legend
  )
}


