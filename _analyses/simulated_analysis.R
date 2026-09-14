###################### Implementation Dip Analysis with Simulated Data ######################

#### Packages ####
library(did)
library(DIDmultiplegtDYN)
library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(cowplot)

source("_analyses/simulated_data.R")
source("_support_functions/twfe_did_mult_outcome_wrapper.R")


###################### Reference series ######################

DIP_W <- c(1.0, 1.25, 0.9, 0.4)   # Matches KERNEL in DGP

ref_series <- function(df, min_e = -4, max_e = 5) {
  d <- df %>% filter(G > 0) %>% mutate(e = time - G) %>%
    filter(e >= min_e, e <= max_e)
  
  # The estimand CS and dCDH target.
  agg <- d %>% group_by(e) %>%
    summarise(estimate = mean(tau_true), .groups = "drop")
  
  # Single-event dip: at e = 0 no treatment yet, so the mean of tau_true there is the first event's magnitude. Scale it by the kernel.
  dip <- data.frame(e = min_e:max_e) %>%
    left_join(data.frame(e = 0:(length(DIP_W) - 1), w = DIP_W), by = "e") %>%
    mutate(estimate = agg$estimate[agg$e == 0] * coalesce(w, 0)) %>%
    select(e, estimate)
  
  bind_rows(
    dip %>% mutate(estimator = "True dip"),
    agg %>% mutate(estimator = "True estimand")
  ) %>%
    mutate(outcome = "Y", se = 0, conf.low = estimate, conf.high = estimate) %>%
    select(outcome, e, estimate, se, conf.low, conf.high, estimator)
}


df_ref_A <- ref_series(simA)
df_ref_B <- ref_series(simB)


###################### Set Up ######################

# CS options
att_options_sim <- list(
  tname                  = "time",
  idname                 = "unit",
  gname                  = "G",
  panel                  = TRUE,
  allow_unbalanced_panel = TRUE,
  clustervars            = "unit",
  control_group          = "nevertreated",
  anticipation           = 0,
  base_period            = "universal",
  bstrap                 = TRUE,
  biters                 = 1000
)

# dCDH options
dcdh_options_sim <- list(
  treatment            = "D_cum",
  group                = "unit",
  time                 = "time",
  cluster              = "unit",
  normalized           = FALSE,
  only_never_switchers = TRUE,
  same_switchers       = FALSE,
  same_switchers_pl    = FALSE,
  graph_off            = TRUE
)

# Normalised dCDH
dcdh_options_simN  <- modifyList(dcdh_options_sim,  list(normalized = TRUE))

# Two reference series: the dip (what H1 is about) and the estimand CS and
# dCDH actually target. The gap between them in Sim B is the result.
est_cols_sim   <- c("True dip"          = "#111111",
                    "True estimand"     = "#1b9e77",
                    "CS"                = "#2b5c8f",
                    "dCDH"              = "#fd7107",
                    "dCDH (norm.)"      = "#7a3e9d")
est_shapes_sim <- c("True dip"          = 15,
                    "True estimand"     = 18,
                    "CS"                = 16,
                    "dCDH"              = 17,
                    "dCDH (norm.)"      = 4)
est_lines_sim  <- c("True dip"          = "dashed",
                    "True estimand"     = "dotted",
                    "CS"                = "solid",
                    "dCDH"              = "solid",
                    "dCDH (norm.)"      = "longdash")


###################### Sim A: single treatment event ######################

# CS
did_cs_A <- run_csdid_pipeline(simA, "Y", label = "Sim A", formula = NULL,
                               options = att_options_sim, gname = "G")
df_cs_A <- extract_cs_df(did_cs_A, "CS")

# dCDH
dcdh_A <- do.call(did_multiplegt_dyn,
                  c(list(df = simA, outcome = "Y", effects = 6, placebo = 3),
                    dcdh_options_sim))
df_dcdh_A <- extract_dcdh_results(dcdh_A, "Y", "dCDH")

# Normalised dCDH
dcdhN_A <- do.call(did_multiplegt_dyn,
                   c(list(df = simA, outcome = "Y", effects = 6, placebo = 3),
                     dcdh_options_simN))
df_dcdhN_A <- extract_dcdh_results(dcdhN_A, "Y", "dCDH (norm.)")


###################### Sim B: multiple treatment events ######################

# CS
did_cs_B <- run_csdid_pipeline(simB, "Y", label = "Sim B", formula = NULL,
                               options = att_options_sim, gname = "G")
df_cs_B <- extract_cs_df(did_cs_B, "CS")

# dCDH
dcdh_B <- do.call(did_multiplegt_dyn,
                  c(list(df = simB, outcome = "Y", effects = 6, placebo = 3),
                    dcdh_options_sim))
df_dcdh_B <- extract_dcdh_results(dcdh_B, "Y", "dCDH")

# normalised dCDH
dcdhN_B <- do.call(did_multiplegt_dyn,
                   c(list(df = simB, outcome = "Y", effects = 6, placebo = 3),
                     dcdh_options_simN))
df_dcdhN_B <- extract_dcdh_results(dcdhN_B, "Y", "dCDH (norm.)")


###################### Overlay plot of both simulations ######################

df_AB <- bind_rows(
  bind_rows(df_ref_A, df_cs_A, df_dcdh_A, df_dcdhN_A) %>% mutate(outcome = "A"),
  bind_rows(df_ref_B, df_cs_B, df_dcdh_B, df_dcdhN_B) %>% mutate(outcome = "B")
)

plot_AB <- create_overlay_plot(
  df_AB,
  outcome_labels = c("A" = "Simulation A",
                     "B" = "Simulation B"),
  est_cols = est_cols_sim, est_shapes = est_shapes_sim, est_lines = est_lines_sim,
  ncol = 2
)
plot_AB

# Single nAVSQ plot
df_norm_B <- bind_rows(df_ref_B, df_cs_B, df_dcdhN_B)

plot_B_norm_only <- create_overlay_plot(
  df_dcdhN_B,
  outcome_labels = c("Y" = "Sim B: normalised dCDH (nAVSQ)"),
  est_cols = est_cols_sim, est_shapes = est_shapes_sim, est_lines = est_lines_sim,
  ncol = 1
)
plot_B_norm_only



