###################### Descriptive Analysis of Panel Data ######################

###################### Load packages & data ######################

# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(sf)
library(panelView)
library(ggplot2)
library(gridExtra)

#### Data ####

df_panel_cov <- readRDS('_data/df_panel_cov.rds')

###################### Technical Check: Treatment Variables ######################

# Table treated vs. non-treated units, and comparing if treat_absorbing/treat_absorbing_cs align
comparison_table <- df_panel_cov %>%
  mutate(
    treat_absorbing_cs = ifelse(treat_absorbing_cs <= 1990 & treat_absorbing == 0, NA, treat_absorbing_cs)
  ) %>%
  group_by(election_year) %>%
  summarize(
    total_observations = n(),
    units_treated = sum(treat_nonabsorbing == 1, na.rm = TRUE),
    new_units_treated = sum(treat_absorbing_cs == election_year, na.rm = TRUE),
    cum_treated_via_CS = sum(treat_absorbing_cs <= election_year, na.rm = TRUE),
    cum_treated_via_cumsum = sum(treat_absorbing == 1, na.rm = TRUE)
  )

print(comparison_table)

# Panel view absorbing treatment
panelview(data = df_panel_cov,
          index = c("ags", "election_year"),
          D = "treat_absorbing",
          main = "Panel View: Staggered & Absorbing Treatment",
          xlab = "Year", 
          ylab = "Municipalities",
          display.all = TRUE)

# Panel view non-absorbing treatment
panelview(data = df_panel_cov,
          index = c("ags", "election_year"),
          D = "treat_nonabsorbing",
          main = "Panel View: Staggered & Non-Absorbing Treatment",
          xlab = "Year", 
          ylab = "Municipalities",
          display.all = TRUE)


###################### Wind Power over the Years ######################
# Plot coded with help of claude.ai

df_wp <- df_panel_cov %>%
  mutate(east_ger = ifelse(substr(ags, 1, 2) %in% c("12", "13", "14", "15", "16"), 1, 0)) %>%
  group_by(election_year) %>%
  summarise(
    total = sum(wind_count_0km, na.rm = TRUE),
    east  = sum(wind_count_0km[east_ger == 1], na.rm = TRUE),
    west  = sum(wind_count_0km[east_ger == 0], na.rm = TRUE)
  ) %>%
  arrange(election_year) %>%
  mutate(
    total_cum = cumsum(total),
    east_cum  = cumsum(east),
    west_cum  = cumsum(west)
  )

colors <- c(total = "#333333", west = "#2166AC", east = "#D6604D")
alphas <- c(total = 1, west = 0.45, east = 0.45)

p1 <- df_wp %>%
  pivot_longer(c(total, east, west), names_to = "region", values_to = "count") %>%
  mutate(region = factor(region, levels = c("total", "west", "east"))) %>%
  ggplot(aes(x = election_year, y = count, colour = region, alpha = region)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = colors, labels = c("Total", "West", "East")) +
  scale_alpha_manual(values = alphas, guide = "none") +
  labs(title = "New Wind Turbines per Year", x = NULL, y = "Number of WTs", colour = NULL) +
  theme_minimal(base_size = 12)

p2 <- df_wp %>%
  pivot_longer(c(total_cum, east_cum, west_cum), names_to = "region", values_to = "count") %>%
  mutate(region = factor(region, levels = c("total_cum", "west_cum", "east_cum"))) %>%
  ggplot(aes(x = election_year, y = count, colour = region, alpha = region)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = setNames(colors, paste0(names(colors), "_cum")),
                      labels = c("Total", "West", "East")) +
  scale_alpha_manual(values = setNames(alphas, paste0(names(alphas), "_cum")), guide = "none") +
  labs(title = "Cumulative Wind Turbines", x = NULL, y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "right")

grid.arrange(p1, p2, ncol = 2)
dev.off()


###################### Election Outcomes over the Years ######################
# Plot coded with help of claude.ai

party_colors <- c(
  cdu       = "#2C2C2C",
  spd       = "#E3000F",
  fdp       = "#FFCC00",
  linke_pds = "#BE3075",
  gruene    = "#64A12D",
  afd       = "#009EE0"
)

df_panel_cov %>%
  select(election_year, cdu, spd, fdp, linke_pds, gruene, afd) %>%
  pivot_longer(-election_year, names_to = "party", values_to = "share") %>%
  group_by(election_year, party) %>%
  summarise(share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = election_year, y = share, colour = party)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = party_colors,
    labels = c(cdu = "CDU", spd = "SPD", fdp = "FDP",
               linke_pds = "Linke/PDS", gruene = "Grüne", afd = "AfD")
  ) +
  labs(title = "Mean vote share by party", x = NULL, y = "Vote share (%)", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")


###################### Spatial distribution Wind Power in Germany ######################
# Plot coded with help of claude.ai

shapefile_path <- '_data/_shapefiles/VG250_GEM.shp'
municipalities <- st_read(shapefile_path) %>%
  mutate(AGS = as.character(AGS))

# Germany outline by dissolving all municipalities
germany_outline <- municipalities %>% st_union()

df_cum <- df_panel_cov %>%
  arrange(ags, election_year) %>%
  group_by(ags) %>%
  mutate(wind_cum = cumsum(wind_count_0km)) %>%
  ungroup() %>%
  mutate(ags = as.character(ags))

cap <- quantile(df_cum$wind_cum, 0.99, na.rm = TRUE)

make_map <- function(year) {
  df_year <- df_cum %>%
    filter(election_year == year) %>%
    select(ags, wind_cum) %>%
    mutate(wind_cum_capped = pmin(wind_cum, cap))
  
  map_data <- municipalities %>%
    left_join(df_year, by = c("AGS" = "ags"))
  
  if (all(is.na(map_data$wind_cum_capped))) {
    df_year <- df_year %>% mutate(ags = str_pad(ags, 8, "left", "0"))
    map_data <- municipalities %>%
      mutate(AGS = str_pad(AGS, 8, "left", "0")) %>%
      left_join(df_year, by = c("AGS" = "ags"))
  }
  
  ggplot() +
    geom_sf(data = map_data, aes(fill = wind_cum_capped), colour = NA) +
    geom_sf(data = germany_outline, fill = NA, colour = "black", linewidth = 0.3) +
    scale_fill_gradient(
      low      = "#f7fbff",
      high     = "#084594",
      na.value = "grey92",
      limits   = c(0, cap),
      name     = "Turbines",
      breaks   = c(0, cap/3, cap/1.5, cap),
      labels   = c("0", round(cap/3), round(cap/1.5), paste0(">", round(cap))),
      guide    = guide_colorbar(barwidth = 0.5, barheight = 4, title.position = "top")
    ) +
    labs(title = year) +
    theme_void(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8)
    )
}

p1 <- make_map(1998)
p2 <- make_map(2009)
p3 <- make_map(2017)
p4 <- make_map(2025)

grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()


###################### Availability of Covariates across periods ######################
# Plot coded with help of claude.ai

years <- c(1990, 1994, 1998, 2002, 2005, 2009, 2013, 2017, 2021, 2025)

covariates <- tribble(
  ~variable,                        ~group,          ~start, ~end,  ~imp_1994, ~imp_2025,
  # Outcomes
  "Turnout / CDU / SPD / FDP / etc","Outcomes",       1990,   2025,  FALSE,     FALSE,
  "AfD",                            "Outcomes",       2013,   2025,  FALSE,     FALSE,
  # Municipality (alphabetical)
  "Agricultural land",              "Municipality",   2017,   2025,  FALSE,     TRUE,
  "Average age",                    "Municipality",   2002,   2025,  FALSE,     TRUE,
  "Business tax",                   "Municipality",   1994,   2025,  TRUE,      TRUE,
  "Commuter balance",               "Municipality",   1998,   2025,  FALSE,     TRUE,
  "Income tax",                     "Municipality",   1994,   2025,  TRUE,      TRUE,
  "Net migration",                  "Municipality",   1994,   2025,  TRUE,      TRUE,
  "Overnight stays",                "Municipality",   2009,   2025,  FALSE,     TRUE,
  "Population density",             "Municipality",   1990,   2025,  FALSE,     FALSE,
  "Purchasing power",               "Municipality",   2013,   2025,  FALSE,     TRUE,
  "Share employed (social sec.)",   "Municipality",   1998,   2025,  FALSE,     TRUE,
  "Share female population",        "Municipality",   1994,   2025,  TRUE,      TRUE,
  "Tax revenue",                    "Municipality",   1994,   2025,  TRUE,      TRUE,
  "Unemployment rate",              "Municipality",   1998,   2025,  FALSE,     TRUE,
  "Workplace density",              "Municipality",   2000,   2025,  FALSE,     TRUE,
  # County (alphabetical)
  "Household income (monthly)",     "County",         2000,   2025,  FALSE,     TRUE,
  "Share foreign population",       "County",         1994,   2025,  TRUE,      TRUE,
  "Vocational qualifications",      "County",         2013,   2025,  FALSE,     TRUE
)

df_plot <- covariates %>%
  mutate(variable = factor(variable, levels = rev(variable))) %>%
  crossing(year = years) %>%
  mutate(
    available = year >= start & year <= end,
    imputed   = available & ((year == 1994 & imp_1994) | (year == 2025 & imp_2025))
  )

group_colors <- c(
  "Outcomes"     = "#4a6fa5",
  "Municipality" = "#4a7c59",
  "County"       = "#8b5e3c"
)

ggplot(df_plot, aes(x = factor(year), y = variable)) +
  geom_tile(aes(fill = group, alpha = available), colour = "white", linewidth = 0.4) +
  geom_text(data = filter(df_plot, available & !imputed),
            label = "o", size = 3, colour = "white") +
  geom_text(data = filter(df_plot, imputed),
            label = "~", size = 3.4, colour = "white") +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.07), guide = "none") +
  labs(title = "Covariate Availability by Election Year",
       subtitle = "o = available    ~ imputed",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_text(size = 8.5, colour = "grey30"),
    axis.text.x     = element_text(size = 8.5, colour = "grey30"),
    legend.position = "top",
    legend.text     = element_text(size = 9, colour = "grey30"),
    plot.title      = element_text(size = 11, face = "bold", colour = "grey20", margin = margin(b = 4)),
    plot.subtitle   = element_text(size = 8.5, colour = "grey50", margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", colour = NA)
  )

