###################### Hyp1 and Hyp2 using dCDH ######################

###################### Set the scene  ######################

#### Packages ####
# install.packages(c("rlang", "S7"))
# install.packages("polars", repos = "https://rpolars.r-universe.dev")
library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(DIDmultiplegtDYN)
library(polars)

#### Helper functions ####
source("_support_functions/twfe_did_mult_outcome_wrapper.R")

#### Data ####
df_did_ready <- readRDS("_data/df_did_ready.rds")

# Add categorical variable "state" (bundesland)
df_did_ready <- df_did_ready %>%
  mutate(
    ags_char   = str_pad(as.character(ags), width = 8, side = "left", pad = "0"),
    bundesland = substr(ags_char, 1, 2)
  )

#### Outcomes ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")


###################### Hyp 1 ######################

# Option 1: Most simple setup, not controls
res_dcdh_1 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 1 (baseline)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = 3,
  placebo         = c(default = 3, afd = 1),
  cluster         = "ags",
  normalized      = FALSE
)

# Individual outcome estimates, e.g. gruene
summary(dcdh_results[["gruene"]])


# Option 2: Controls as in CS DiD
res_dcdh_2 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 2 (with controls)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = 3,
  placebo         = c(default = 3, afd = 1),
  controls        = "pop_density",
  trends_nonparam = "east_ger",
  cluster         = "ags",
  normalized      = TRUE
)


# Option 3: Only units with same trajectory
res_dcdh_3 <- run_dcdh_pipeline(
  data            = df_did_ready,
  outcomes        = outcome_vars,
  label           = "Model 3 (same trajectory)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = 3,
  placebo         = c(default = 3, afd = 1),
  controls        = "pop_density",
  trends_nonparam = "east_ger",
  cluster         = "ags",
  normalized      = TRUE,
  same_switchers = TRUE,
  same_switchers_pl = TRUE
)


# Open Question:
# Normalisation yes/no
# Which controls and how (trends_nonparam vs. controls)?
# Same switchers/same_witchers_pl to ensure fully balanced panel across pre- and post periods? 


