###################### Combine GERDA Election Data with Wind Power Count Data ######################

###################### Load packages & data ######################

# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(gerda)

#### Data ####
# Election results
gerda_data_list(print_table = TRUE)
df_elec_pre <- load_gerda_web("federal_muni_harm_25")

# Wind power counts
df_wind_pre <- readRDS('_data/df_wp_count.rds')

# Covariates
df_covariates <- readRDS('_data/df_covariates.rds')


###################### Filter df_elec and df_wind for columns of interest ######################

df_elec <- df_elec_pre %>%
  mutate(pop_density = (population/area)*1000) %>%
  select(ags, ags_name, election_year, county, turnout, cdu, csu, spd, fdp, linke_pds, gruene, afd,
         cdu_csu, far_right, far_left, far_left_w_linke, freie_wahler, population, pop_density)

df_wind <- df_wind_pre %>%
  rename(ags = AGS, election_year = year) %>%
  select(ags, GEN, election_year, wind_count_0km, wind_count_1km, wind_count_3km,
         wind_count_5km, wind_count_10km, cum_wind_count_3km, lag_wind_count_3km,
         cum_lag_wind_count_3km, wind_count_muni, cum_wind_count_muni, lag_wind_count_muni,
         cum_lag_wind_count_muni)


###################### Check ags overlap between df_elec and df_wind ######################

cat("Unique ags in df_elec:", length(unique(df_elec$ags)), "\n")
cat("Unique ags in df_wind:", length(unique(df_wind$ags)), "\n")
cat("Matching ags:", sum(unique(df_elec$ags) %in% unique(df_wind$ags)), "\n")

# Which df_elec ags are unmatched by df_wind?
unmatched_elec <- unique(df_elec$ags)[!unique(df_elec$ags) %in% unique(df_wind$ags)]
cat("Unmatched ags in df_elec:", length(unmatched_elec), "\n")
# What are the unmatched ags in df_elec?
df_elec %>%
  filter(ags %in% unmatched_elec) %>%
  distinct(ags, ags_name) %>%
  print(n = Inf)
# Obergeckler = 07232503

# Which df_wind ags not in df_elec?
unmatched_wind <- unique(df_wind$ags)[!unique(df_wind$ags) %in% unique(df_elec$ags)]
cat("df_wind ags not in df_elec:", length(unmatched_wind), "\n")
unmatched_df_wind <- df_wind %>%
  filter(ags %in% unmatched_wind) %>%
  distinct(ags, GEN)
# 205 units, all of which are uninhabited or have very low number of inhabitants (e.g. Wiedenborstel)

# Can savely proceed in merging without loosing units
# But fix ags of Obergeckler in df_wind and name in df_elec
df_wind <- df_wind %>% mutate(ags = ifelse(trimws(toupper(GEN)) == "OBERGECKLER", "07232503", ags))
df_elec <- df_elec %>% 
  mutate(ags_name = case_when(
    ags == "07232503" ~ "Obergeckler",
    ags == "16076094" ~ "Berga-Wünschendorf, Stadt",
    ags == "16061119" ~ "Uder",
    TRUE ~ ags_name
  ))


###################### Join df_elec and df_wind ######################
# This leads to 205 municipalities in df_wind being dropped
# 204 of them are uninhabited, one (Wiedenborstel) has only 10 inhabitants

df_panel_pre <- df_elec %>% left_join(df_wind %>% select(-GEN), by = c("ags", "election_year"))


###################### Construct treatment and other variables ######################

# Absorbing treatment indicator
df_panel_wo_cov <- df_panel_pre %>%
  mutate(
    dummy_lag_wind_count_3km = ifelse(cum_lag_wind_count_3km > 0, 1, 0),
    cat_cum_lag_wind_count_3km = cut(cum_lag_wind_count_3km,
                                     breaks = c(-1, 0, seq(10, 100, by = 10), Inf),
                                     labels = c("0", "1–10", "11–20", "21–30", "31–40", "41–50",
                                                "51–60", "61–70", "71–80", "81–90", "91–100", ">100"),
                                     right = TRUE),
    treat_nonabsorbing = as.integer(wind_count_3km > 0)) %>%
  arrange(ags, election_year) %>%
  group_by(ags) %>%
  mutate(treat_absorbing = as.integer(cumsum(wind_count_3km) > 0)) %>%
  ungroup()


###################### Save data ######################
#saveRDS(df_panel_wo_cov, '_data/df_panel_wo_cov.rds')


###################### Add covariates ######################

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
# This leads to 229 municipalities in df_covariates being dropped
# 204 of them are uninhabited, one (Wiedenborstel) has only 10 inhabitants
# 19 are municipalities merged or split between 2023 and 2025 harmonisation standards, and crosswalks are not yet available
# Most importantly, no municipalities for which election results and wind power counts exist are dropped
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
    share_fem   = pop_fem/(population*1000),
    share_unemp = unemp/(population*1000),
    share_emp = employees/(population*1000)
  )

###################### Save data ######################
#saveRDS(df_panel_cov, '_data/df_panel_cov.rds')








