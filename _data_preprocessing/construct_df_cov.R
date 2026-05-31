###################### Create Covariate Data ######################
#' Covariates Municipality Level:
#' - average age, 2002-2023
#' - population female, 1995-2023
#' - number of unemployed, 1998-2023
#' 
#' 
#' Covariates County Level:
#' - average monthly household income, 2000-2022
#' - share of foreign populationm 1995-2023
#' - share of working population with vocational qualifications, 2012-2023
#' 
#' 
#' Covariates Constructed:
#' - share of female population, 1995-2023
#' - unemployment rate, 1998-2023

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
    TRUE               ~ as.integer(year)
  )
}


###################### Municipality Level ######################

file_list_muni <- list(
  average_age = "_data/_covariate_data/average_age.csv",
  pop_fem = "_data/_covariate_data/pop_fem.csv",
  unemp = "_data/_covariate_data/unemp.csv"
)

cov_list_muni <- imap(file_list_muni, function(path, cov_name) {
  
  df_raw <- read.csv2(path, skip = 1)
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
  
  df_raw <- read.csv2(path, skip = 1)
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


###################### Merge covariate data to df_panel_wo_cov ######################

# Check ags overlap in both directions
# In df_panel_wo_cov but not in df_covariates
cat("AGS in df_panel_wo_cov but NOT in df_covariates:\n")
anti_join(df_panel_wo_cov, df_covariates, by = "ags") %>% distinct(ags) %>% print()

# In df_covariates but not in panel:
in_cov_not_panel <- unique(df_covariates$ags)[!unique(df_covariates$ags) %in% unique(df_panel_wo_cov$ags)]
cat("AGS in df_covariates but NOT in df_panel_wo_cov:", length(in_cov_not_panel), "\n")
ags_in_cov_not_panel <- anti_join(df_covariates, df_panel_wo_cov, by = "ags") %>% 
  distinct(ags, ags_name)

# Fix wrong ags of "Obergeckler", "Uder" and "Berga/Elster, Stadt" in df_av_age
df_covariates <- df_covariates %>%
  mutate(ags = case_when(
    ags_name == "Obergeckler" ~ "07232503",
    ags_name == "Uder" ~ "16061119",
    ags_name == "Berga/Elster, Stadt" ~ "16076094",
    TRUE ~ ags
  ))

# Left join
df_panel_cov <- df_panel_wo_cov %>%
  left_join(df_covariates %>% select(-ags_name),
            by = c("ags", "election_year"))

# Verify NA patterns look right
df_panel_cov %>%
  group_by(election_year) %>%
  summarise(n_missing = sum(is.na(share_foreign)),
            n_total = n()) %>%
  print()


###################### Calculate Shares ######################

df_panel_cov <- df_panel_cov %>%
  mutate(
    share_fem   = pop_fem / (population*1000),
    share_unemp = unemp / (population*1000)
  )


