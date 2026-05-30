###################### Data pre-processing ######################

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
df_wind_pre <- readRDS('/Users/Nina/Documents/_Daten/_Studium/Master SDS/4. Sommer 2026/P5 Master Thesis/_data/df_wp_count.rds')

# Covariates



###################### Filter df_elec and df_wind for columns of interest ######################

df_elec <- df_elec_pre %>%
  mutate(pop_density = (population/area)*1000) %>%
  select(ags, ags_name, election_year, county, turnout, cdu, csu, spd, fdp, linke_pds, gruene, afd,
         cdu_csu, far_right, far_left, far_left_w_linke, freie_wahler, pop_density)

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
df_elec <- df_elec %>% mutate(ags_name = ifelse(ags == "07232503", "Obergeckler", ags_name))


###################### Join df_elec and df_wind ######################

df_panel <- df_elec %>% left_join(df_wind %>% select(-GEN), by = c("ags", "election_year"))

# Save data
#saveRDS(df_panel, '/Users/Nina/Documents/_Daten/_Studium/Master SDS/4. Sommer 2026/P5 Master Thesis/_data/df_panel_wo_cov.rds')


###################### Add covariates ######################






