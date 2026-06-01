###################### Create Covariate Data ######################
#' Covariates Municipality Level:
#' - average age, 2002-2023
#' - population female, 1995-2023
#' 
#' - income tax, 1995-2023
#' - business tax, 1995-2023
#' - tax revenue, 1995-2023
#' 
#' - employees subject to social security contributions in their place of residence, 1998-2023
#' - workplace density, 2000-2023
#' - number of unemployed, 1998-2023
#' - commuter balance, 1998-2023
#' 
#' - number of overnight stays in tourist accommodations, 2009-2023
#' - net migration, 1995-2023
#' - purchasing power, 2013-2023
#' - agricultural land, 2017-2023
#' 
#' Covariates County Level:
#' - average monthly household income, 2000-2022
#' - share of foreign population 1995-2023
#' - share of working population with vocational qualifications, 2012-2023
#' 
#' Covariates Constructed:
#' - share of female population, 1995-2023
#' - unemployment rate, 1998-2023
#' - share of population subject to social security contributions through employment, 1998-2023

###################### Prep ######################

# Packages
library(tidyverse)


# Define inputs for loop
target_years <- c(1994, 1998, 2002, 2005, 2009, 2013, 2017, 2021, 2025)

year_recode <- function(year) {
  case_when(
    year > 2021        ~ 2025L,
    year %in% c(1999, 2000) ~ 1998L,
    year %in% c(1995, 1996) ~ 1994L,
    year %in% c(2004) ~ 2005L,
    TRUE               ~ as.integer(year)
  )
}


###################### Municipality Level ######################

file_list_muni <- list(
  average_age = "_data/_covariate_data/average_age.csv",
  pop_fem = "_data/_covariate_data/pop_fem.csv",
  unemp = "_data/_covariate_data/unemp.csv",
  
  inc_tax = "_data/_covariate_data/inc_tax.csv",
  busi_tax = "_data/_covariate_data/busi_tax.csv",
  tax_rev = "_data/_covariate_data/tax_rev.csv",
  
  agri_land = "_data/_covariate_data/agri_land.csv",
  commute_balance = "_data/_covariate_data/commute_balance.csv",
  dens_work = "_data/_covariate_data/dens_work.csv",
  employees = "_data/_covariate_data/employees.csv",
  net_migration = "_data/_covariate_data/net_migration.csv",
  purch_pow = "_data/_covariate_data/purch_pow.csv",
  tourism = "_data/_covariate_data/tourism.csv"
)

cov_list_muni <- imap(file_list_muni, function(path, cov_name) {
  
  df_raw <- read.csv2(path, skip = 1, colClasses = "character")
  df_raw$ags <- str_pad(as.character(df_raw$X), width = 8, side = "left", pad = "0")
  
  df_long <- df_raw %>%
    rename(ags_name = X.1) %>%
    select(-X, -X.2) %>%
    pivot_longer(
      cols      = starts_with("X"),
      names_to  = "election_year",
      values_to = cov_name
    ) %>%
    mutate(
      election_year = as.integer(str_remove(election_year, "X")),
      election_year = year_recode(election_year),
      !!cov_name   := trimws(as.character(.data[[cov_name]])),
      !!cov_name   := na_if(as.character(.data[[cov_name]]), "-"),
      !!cov_name   := as.numeric(gsub(",", ".", gsub("\\.", "", gsub("\\s+", "", .data[[cov_name]]))))
    ) %>%
    filter(election_year %in% target_years) %>%
    select(ags, ags_name, election_year, all_of(cov_name))
  
  all_combos <- expand_grid(
    ags           = unique(df_long$ags),
    election_year = target_years
  ) %>%
    left_join(df_long %>% distinct(ags, ags_name), by = "ags")
  
  all_combos %>%
    left_join(df_long, by = c("ags", "ags_name", "election_year"))
})

df_covariates_muni <- cov_list_muni %>%
  reduce(full_join, by = c("ags", "ags_name", "election_year")) %>%
  mutate(county_ags = substr(ags, 1, 5),
         average_age = average_age / 100)


###################### County Level ######################

file_list_county <- list(
  hinc = "_data/_covariate_data/hinc.csv",
  share_foreign = "_data/_covariate_data/share_foreign.csv",
  educ = "_data/_covariate_data/educ.csv"
)

cov_list_county <- imap(file_list_county, function(path, cov_name) {
  
  df_raw <- read.csv2(path, skip = 1, colClasses = "character")
  df_raw$county_ags <- str_pad(as.character(df_raw$X), width = 5, side = "left", pad = "0")
  
  df_long <- df_raw %>%
    rename(county_name = X.1) %>%
    select(-X, -X.2) %>%
    pivot_longer(
      cols      = starts_with("X"),
      names_to  = "election_year",
      values_to = cov_name
    ) %>%
    mutate(
      election_year = as.integer(str_remove(election_year, "X")),
      election_year = year_recode(election_year),
      !!cov_name   := trimws(as.character(.data[[cov_name]])),
      !!cov_name   := na_if(as.character(.data[[cov_name]]), "-"),
      !!cov_name   := as.numeric(gsub(",", ".", gsub("\\.", "", gsub("\\s+", "", .data[[cov_name]]))))
    ) %>%
    filter(election_year %in% target_years) %>%
    select(county_ags, county_name, election_year, all_of(cov_name))
  
  all_combos <- expand_grid(
    county_ags           = unique(df_long$county_ags),
    election_year = target_years
  ) %>%
    left_join(df_long %>% distinct(county_ags, county_name), by = "county_ags")
  
  all_combos %>%
    left_join(df_long, by = c("county_ags", "county_name", "election_year"))
})

df_covariates_county <- cov_list_county %>%
  reduce(full_join, by = c("county_ags", "county_name", "election_year")) %>%
  mutate(share_foreign = share_foreign/10000,
         educ = educ/10000)


###################### Merge municipality and county level data ######################

df_covariates <- df_covariates_muni %>%
  left_join(df_covariates_county, by = c("county_ags", "election_year")) %>%
  select(-county_name, -county_ags)

#saveRDS(df_covariates, '_data/df_covariates.rds')

