############ Count Variable Wind Power Plants - All Periods ############

########### Set Up ###########

#### Load packages ####
library(sf)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

#### Load municipality spatial data files ####
shapefile_path <- '/Users/Nina/Documents/_Daten/_Studium/Master SDS/4. Sommer 2026/P5 Master Thesis/_data/_shapefiles/VG250_GEM.shp'
municipalities <- st_read(shapefile_path)
municipalities <- municipalities %>% mutate(AGS = as.character(AGS))
filtered_municipalities <- municipalities

#### Load wind power plant unit data ####
wind_units <- readRDS('/Users/Nina/Documents/_Daten/_Studium/Master SDS/4. Sommer 2026/P5 Master Thesis/_data/df_wind_units_all.rds')

########### Construct Count Variable ###########

wind_units <- wind_units %>%
  mutate(Inbetriebnahmedatum.der.Einheit = as.Date(Inbetriebnahmedatum.der.Einheit, format = "%d.%m.%Y"))

# Define time periods (extended back to 1990)
time_periods <- list(
  "1990" = c(NA,                       as.Date("1990-12-02")),
  "1994" = c(as.Date("1990-12-02"),    as.Date("1994-10-16")),
  "1998" = c(as.Date("1994-10-16"),    as.Date("1998-09-27")),
  "2002" = c(as.Date("1998-09-27"),    as.Date("2002-09-22")),
  "2005" = c(as.Date("2002-09-22"),    as.Date("2005-09-18")),
  "2009" = c(as.Date("2005-09-18"),    as.Date("2009-09-27")),
  "2013" = c(as.Date("2009-09-27"),    as.Date("2013-09-22")),
  "2017" = c(as.Date("2013-09-22"),    as.Date("2017-09-24")),
  "2021" = c(as.Date("2017-09-24"),    as.Date("2021-09-26")),
  "2025" = c(as.Date("2021-09-26"),    NA))

# Buffer distances
buffer_distances <- c(0, 1000, 3000, 5000, 10000)
buffer_labels    <- c("0km", "1km", "3km", "5km", "10km")

results_list <- list()

for (year_label in names(time_periods)) {
  
  period     <- time_periods[[year_label]]
  start_date <- period[1]
  end_date   <- period[2]
  
  wind_data <- wind_units %>%
    filter(is.na(start_date) | Inbetriebnahmedatum.der.Einheit >= start_date) %>%
    filter(is.na(end_date)   | Inbetriebnahmedatum.der.Einheit <  end_date)
  
  wind_sf <- st_as_sf(wind_data, coords = c("longitude", "latitude"), crs = 4326)
  wind_sf <- st_transform(wind_sf, crs = 3857)
  
  filtered_municipalities <- st_transform(filtered_municipalities, crs = st_crs(wind_sf))
  
  # Compute intersection counts for each buffer distance
  for (i in seq_along(buffer_distances)) {
    
    dist  <- buffer_distances[i]
    label <- buffer_labels[i]
    
    if (dist == 0) {
      wind_geom <- wind_sf
    } else {
      wind_geom <- st_buffer(wind_sf, dist = dist, endCapStyle = "ROUND")
    }
    
    intersections <- st_intersects(filtered_municipalities, wind_geom, sparse = FALSE)
    filtered_municipalities[[paste0("wind_count_", label)]] <- rowSums(intersections)
  }
  
  filtered_municipalities_df <- filtered_municipalities %>%
    mutate(longitude = st_coordinates(st_centroid(geometry))[,1],
           latitude  = st_coordinates(st_centroid(geometry))[,2]) %>%
    st_drop_geometry() %>%
    mutate(year = as.numeric(year_label)) %>%
    group_by(AGS) %>%
    arrange(desc(year)) %>%
    slice(1) %>%
    ungroup()
  
  results_list[[year_label]] <- filtered_municipalities_df
}

#### Combine all years and create wind count variables ####

all_periods <- c("1990", "1994", "1998", "2002", "2005", "2009", "2013", "2017", "2021", "2025")

wind_long_df <- bind_rows(results_list) %>%
  arrange(AGS, year) %>%
  group_by(AGS) %>%
  mutate(
    # Cumulative and lag only for 3km buffer
    cum_wind_count_3km     = cumsum(wind_count_3km),
    lag_wind_count_3km     = lag(wind_count_3km,          order_by = year),
    cum_lag_wind_count_3km = lag(cumsum(wind_count_3km),  order_by = year)) %>%
  ungroup()


########### Robustness Check: WP Count at Municipality Level (No Buffer) ###########

wind_units_muni <- wind_units %>%
  rename(construction_date = Inbetriebnahmedatum.der.Einheit,
         AGS = Gemeindeschlüssel) %>%
  mutate(
    construction_date = as.Date(construction_date, format = "%Y-%m-%d"),
    AGS = str_pad(as.character(AGS), width = 8, side = "left", pad = "0"))

wind_muni_df <- wind_units_muni %>%
  mutate(year = case_when(
    construction_date < as.Date("1990-12-02") ~ "1990",
    construction_date >= as.Date("1990-12-02") & construction_date < as.Date("1994-10-16") ~ "1994",
    construction_date >= as.Date("1994-10-16") & construction_date < as.Date("1998-09-27") ~ "1998",
    construction_date >= as.Date("1998-09-27") & construction_date < as.Date("2002-09-22") ~ "2002",
    construction_date >= as.Date("2002-09-22") & construction_date < as.Date("2005-09-18") ~ "2005",
    construction_date >= as.Date("2005-09-18") & construction_date < as.Date("2009-09-27") ~ "2009",
    construction_date >= as.Date("2009-09-27") & construction_date < as.Date("2013-09-22") ~ "2013",
    construction_date >= as.Date("2013-09-22") & construction_date < as.Date("2017-09-24") ~ "2017",
    construction_date >= as.Date("2017-09-24") & construction_date < as.Date("2021-09-26") ~ "2021",
    construction_date >= as.Date("2021-09-26") ~ "2025")) %>%
  group_by(AGS, year) %>%
  summarise(wind_count_muni = n(),
            Gemeinde = first(Gemeinde),
            .groups = "drop") %>%
  complete(AGS, year = all_periods, fill = list(wind_count_muni = 0)) %>%
  arrange(AGS, year) %>%
  mutate(year = as.integer(year)) %>%
  group_by(AGS) %>%
  mutate(Gemeinde                = first(na.omit(Gemeinde)),
         cum_wind_count_muni     = cumsum(wind_count_muni),
         lag_wind_count_muni     = lag(wind_count_muni,         default = 0),
         cum_lag_wind_count_muni = lag(cum_wind_count_muni,     default = 0)) %>%
  ungroup()

# Two-stage join
unmatched <- unique(wind_muni_df$AGS)[!unique(wind_muni_df$AGS) %in% unique(wind_long_df$AGS)]

joined_stage1 <- wind_long_df %>%
  left_join(
    wind_muni_df %>%
      select(AGS, year, wind_count_muni, cum_wind_count_muni,
             lag_wind_count_muni, cum_lag_wind_count_muni),
    by = c("AGS", "year"))

name_lookup <- wind_muni_df %>%
  filter(AGS %in% unmatched) %>%
  mutate(Gemeinde = trimws(toupper(Gemeinde))) %>%
  select(Gemeinde, year, wind_count_muni, cum_wind_count_muni,
         lag_wind_count_muni, cum_lag_wind_count_muni)

# Inspect name matches
joined_stage1 %>%
  mutate(GEN = trimws(toupper(GEN))) %>%
  inner_join(name_lookup %>% distinct(Gemeinde), by = c("GEN" = "Gemeinde")) %>%
  select(AGS, GEN, year, wind_count_muni) %>%
  arrange(GEN, year)

# Document unresolvable municipalities
unresolved <- wind_muni_df %>%
  filter(AGS %in% unmatched) %>%
  mutate(Gemeinde_upper = trimws(toupper(Gemeinde))) %>%
  filter(!Gemeinde_upper %in% trimws(toupper(wind_long_df$GEN))) %>%
  distinct(AGS, Gemeinde) %>%
  left_join(
    wind_units_muni %>%
      group_by(AGS) %>%
      summarise(n_turbines         = n(),
                first_construction = min(construction_date),
                last_construction  = max(construction_date),
                .groups = "drop"),
    by = "AGS")

cat("Municipalities unresolvable by AGS or name:", nrow(unresolved), "\n")
cat("Total turbines lost:", sum(unresolved$n_turbines), "\n")

# AGS with counts from wrong municipality (matched by name, non-NA from stage 1)
wrong_ags <- joined_stage1 %>%
  mutate(GEN = trimws(toupper(GEN))) %>%
  inner_join(name_lookup %>% distinct(Gemeinde), by = c("GEN" = "Gemeinde")) %>%
  filter(!is.na(wind_count_muni)) %>%
  distinct(AGS) %>%
  pull(AGS)

# Apply stage 2 fill-in
wind_long_df <- joined_stage1 %>%
  mutate(GEN = trimws(toupper(GEN))) %>%
  left_join(name_lookup, by = c("GEN" = "Gemeinde", "year")) %>%
  mutate(
    wind_count_muni = case_when(
      AGS %in% wrong_ags ~ rowSums(cbind(wind_count_muni.x, wind_count_muni.y), na.rm = TRUE),
      TRUE ~ coalesce(wind_count_muni.x, wind_count_muni.y)),
    cum_wind_count_muni = case_when(
      AGS %in% wrong_ags ~ rowSums(cbind(cum_wind_count_muni.x, cum_wind_count_muni.y), na.rm = TRUE),
      TRUE ~ coalesce(cum_wind_count_muni.x, cum_wind_count_muni.y)),
    lag_wind_count_muni = case_when(
      AGS %in% wrong_ags ~ rowSums(cbind(lag_wind_count_muni.x, lag_wind_count_muni.y), na.rm = TRUE),
      TRUE ~ coalesce(lag_wind_count_muni.x, lag_wind_count_muni.y)),
    cum_lag_wind_count_muni = case_when(
      AGS %in% wrong_ags ~ rowSums(cbind(cum_lag_wind_count_muni.x, cum_lag_wind_count_muni.y), na.rm = TRUE),
      TRUE ~ coalesce(cum_lag_wind_count_muni.x, cum_lag_wind_count_muni.y))) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  mutate(
    wind_count_muni         = ifelse(is.na(wind_count_muni), 0, wind_count_muni),
    cum_wind_count_muni     = ifelse(is.na(cum_wind_count_muni), 0, cum_wind_count_muni),
    lag_wind_count_muni     = ifelse(year == 1990, NA, ifelse(is.na(lag_wind_count_muni), 0, lag_wind_count_muni)),
    cum_lag_wind_count_muni = ifelse(year == 1990, NA, ifelse(is.na(cum_lag_wind_count_muni), 0, cum_lag_wind_count_muni)))


# Save data
#saveRDS(wind_long_df, '/Users/Nina/Documents/_Daten/_Studium/Master SDS/4. Sommer 2026/P5 Master Thesis/_data/df_wp_count.rds')

