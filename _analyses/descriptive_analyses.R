###################### Descriptive Analysis of Panel Data ######################

###################### Load packages & data ######################

# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(panelView)

#### Data ####

df_panel_cov <- readRDS('_data/df_panel_cov.rds')


###################### Treatment Variables ######################

# Table treated vs. non-treated units, and comparing if treat_absorbing/treat_absorbing_cs align
comparison_table <- df_panel_cov %>%
  mutate(
    treat_absorbing_cs = ifelse(treat_absorbing_cs <= 1990 & treat_absorbing == 0, NA, treat_absorbing_cs)
  ) %>%
  group_by(election_year) %>%
  summarize(
    total_observations = n(),
    units_treated_via_CS = sum(treat_absorbing_cs <= election_year, na.rm = TRUE),
    units_treated_via_cumsum = sum(treat_absorbing == 1, na.rm = TRUE)
  )

print(comparison_table)

# Panel view
panelview(data = df_panel_cov,
          index = c("ags", "election_year"),
          D = "treat_absorbing",
          main = "Panel View: Staggered & Absorbing Treatment",
          xlab = "Year", 
          ylab = "Municipalities",
          display.all = TRUE)






