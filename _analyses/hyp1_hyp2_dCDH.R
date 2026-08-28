###################### Hyp1 and Hyp2 using dCDH ######################

###################### Set the scene  ######################

#### Packages ####
# install.packages(c("rlang", "S7"))
# install.packages("polars", repos = "https://rpolars.r-universe.dev")
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(did)
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

# Create static baseline strata for population density
df_did_ready <- df_did_ready %>%
  group_by(ags) %>%
  mutate(pop_density_1990 = pop_density[seq_time == min(seq_time)]) %>%
  ungroup() %>%
  mutate(
    pop_quartile_1990 = ntile(pop_density_1990, 4),
    state_pop_strata  = paste(bundesland, pop_quartile_1990, sep = "_")
  )


#### Outcomes and Controls ####

outcome_vars <- c("turnout", "cdu_csu", "spd", "fdp", 
                  "linke_pds", "gruene", "afd", "current_incumbent")



###################### Hyp 1 ######################

# Option 1: Most simple setup, no controls
res_dcdh_1 <- run_dcdh_pipeline(
  df            = df_did_ready,
  outcome        = outcome_vars,
  label           = "Model 1 (Baseline)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = c(default = 5, afd = 3),
  placebo         = c(default = 4, afd = 1),
  cluster         = "ags",
  normalized      = TRUE
)

# Individual outcome estimates, e.g. gruene
summary(res_dcdh_1[["gruene"]])


# Option 2: Controls as in CS DiD
res_dcdh_2 <- run_dcdh_pipeline(
  df            = df_did_ready,
  outcome        = outcome_vars,
  label           = "Model 2 (With Controls)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = c(default = 5, afd = 3),
  placebo         = c(default = 4, afd = 1),
  controls        = "pop_density",
  trends_nonparam = "bundesland",
  by              = "east_ger",
  cluster         = "ags",
  normalized      = TRUE
)


# Option 3: Only units with same trajectory
res_dcdh_3 <- run_dcdh_pipeline(
  df            = df_did_ready,
  outcome        = outcome_vars,
  label           = "Model 3 (Only Complete Observations)",
  treatment       = "cum_wind_count_3km",
  group           = "ags",
  time            = "seq_time",
  effects         = 3,
  placebo         = c(default = 3, afd = 0),
  controls        = "pop_density",
  trends_nonparam = "east_ger",
  cluster         = "ags",
  normalized      = TRUE,
  same_switchers = TRUE,
  same_switchers_pl = TRUE
)


# Investigate further:
# Normalised weights
# Treatment paths


