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
library(patchwork)
library(cowplot)


#### Data ####

df_panel_cov <- readRDS('_data/df_panel_cov.rds')
df_did_ready <- readRDS('_data/df_did_ready.rds')


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

# Panel View: Staggered and Absorbing Treatment
panel <- panelview(data = df_panel_cov,
          index = c("ags", "election_year"),
          D = "treat_absorbing",
          main = "",
          xlab = "Year", 
          ylab = "Municipalities",
          display.all = TRUE,
          cex.main = 10,
          cex.lab = 18,
          cex.axis = 18,
          cex.legend = 18)

# Panel view non-absorbing treatment
panelview(data = df_panel_cov,
          index = c("ags", "election_year"),
          D = "treat_nonabsorbing",
          main = "Panel View: Staggered & Non-Absorbing Treatment",
          xlab = "Year", 
          ylab = "Municipalities",
          display.all = TRUE)

# ggsave(filename = "_results//descriptive_results/plot_panel_treat.png", plot = panel, width = 12, height = 8, dpi = 300)


###################### Overlap of Propensity Scores for CS-DiD by Cohort ######################

#### Formula for conditional parallel trends ####
covariates_formula <- ~ pop_density + east_ger

# Initialize list
plot_list <- list()

# Loop through cohorts 2 to 10
for (g in 2:10) {
  # Set period to the one before treatment
  current_seq_time <- g - 1
  # Filter data for specific cohort and never-treated group at correct pre-treatment time
  df_filtered <- df_did_ready %>% 
    filter(seq_time == current_seq_time & (seq_group == g | seq_group == 0))
  # Estimate propensity score
  ps_model_g <- glm(ifelse(seq_group == g, 1, 0) ~ pop_density + east_ger,
                    data = df_filtered, family = binomial())
  # Predict propensity scores and assign to df
  df_filtered$pscore <- predict(ps_model_g, type = "response")
  # Create plot
  p <- ggplot(df_filtered, aes(x = pscore, fill = factor(ifelse(seq_group > 0, 1, 0)))) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c("#2a4d7c", "#e2b13c"),
                      labels = c("Control", "Treated"),
                      name = "") +
    labs(title = paste("Cohort", g, "(Period", current_seq_time, ")"),
         x = "Propensity Score", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") # Remove individual legends
  # Store plot in list
  plot_list[[g - 1]] <- p
}

# Extract the legend from one of the plots
shared_legend_plot <- ggplot(df_filtered, aes(x = pscore, fill = factor(ifelse(seq_group > 0, 1, 0)))) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("#2a4d7c", "#e2b13c"),
                    labels = c("Never Treated", "Treated"),
                    name = "") +
  theme_minimal() +
  theme(legend.position = "bottom")
shared_legend <- get_legend(shared_legend_plot)

# Arrange 3x3 grid
main_grid <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 3)

# Combine main grid and shared legend
final_ps_plot <- plot_grid(main_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.05))
# Display the final plot
final_ps_plot

# ggsave(filename = "_results//descriptive_results/plot_ps_overlap.png", plot = final_ps_plot, width = 12, height = 8, dpi = 300)


###################### Wind Power over the Years ######################

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

# Define shared styling vectors
colors <- c(total = "#333333", west = "#084594", east = "#e2b13c")
alphas <- c(total = 1, west = 0.45, east = 0.45)

# --- Plot 1 ---
p1 <- df_wp %>%
  pivot_longer(c(total, east, west), names_to = "region", values_to = "count") %>%
  mutate(region = factor(region, levels = c("total", "west", "east"))) %>%
  ggplot(aes(x = election_year, y = count, colour = region, alpha = region)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  scale_colour_manual(values = colors, labels = c("Total", "West", "East")) +
  scale_alpha_manual(values = alphas, guide = "none") +
  labs(title = "A.)", x = NULL, y = "Number of WTs", colour = NULL) +
  theme_minimal(base_size = 20)

# --- Plot 2 ---
p2 <- df_wp %>%
  pivot_longer(c(total_cum, east_cum, west_cum), names_to = "region", values_to = "count") %>%
  mutate(region = case_when(
    region == "total_cum" ~ "total",
    region == "west_cum"  ~ "west",
    region == "east_cum"  ~ "east"
  )) %>%
  mutate(region = factor(region, levels = c("total", "west", "east"))) %>%
  ggplot(aes(x = election_year, y = count, colour = region, alpha = region)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  # FIX: Now using the exact same color vector and scale as Plot 1
  scale_colour_manual(values = colors, labels = c("Total", "West", "East")) +
  scale_alpha_manual(values = alphas, guide = "none") +
  labs(title = "B.)", x = NULL, y = NULL, colour = NULL) +
  theme_minimal(base_size = 20)

# --- Combine using patchwork ---
plot_wt_dev <- (p1 + p2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# ggsave(filename = "_results//descriptive_results/plot_wt_dev.png", plot = plot_wt_dev, width = 15, height = 8, dpi = 300)


###################### Election Outcomes over the Years ######################

party_colors <- c(
  cdu_csu   = "#2C2C2C",
  spd       = "#E3000F",
  fdp       = "#FFCC00",
  linke_pds = "#BE3075",
  gruene    = "#64A12D",
  afd       = "#009EE0"
)

df_panel_cov %>%
  select(election_year, cdu_csu, spd, fdp, linke_pds, gruene, afd) %>%
  pivot_longer(-election_year, names_to = "party", values_to = "share") %>%
  group_by(election_year, party) %>%
  summarise(share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = election_year, y = share, colour = party)) +
  geom_hline(yintercept = 0.05, colour = "grey", linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = party_colors,
    labels = c(cdu_csu = "CDU/CSU", spd = "SPD", fdp = "FDP", linke_pds = "Linke/PDS", gruene = "Grüne", afd = "AfD")
  ) +
  labs(x = NULL, y = "Vote share (%)", colour = NULL) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom")

# ggsave(filename = "_results//descriptive_results/plot_voteshare_panel.png", plot = plot_voteshare_panel, width = 12, height = 8, dpi = 300)


###################### Spatial distribution Wind Power in Germany ######################

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
      guide    = guide_colorbar(barwidth = 0.8, barheight = 6, title.position = "top")
    ) +
    labs(title = year) +
    theme_void(base_size = 18) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position = "right",
      legend.text     = element_text(size = 11),
      legend.title    = element_text(size = 12)
    )
}

p1 <- make_map(1998)
p2 <- make_map(2009)
p3 <- make_map(2017)
p4 <- make_map(2025)

plot_wt_dist <- grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
# ggsave(filename = "_results//descriptive_results/plot_wt_dist.png", plot = plot_wt_dist, width = 12, height = 8, dpi = 300)



###################### Availability of Covariates across periods ######################

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
    imputed   = available & ((year == 1994 & imp_1994) | (year == 2025 & imp_2025)),
    status    = case_when(
      imputed ~ "imputed",
      available & !imputed ~ "available",
      TRUE ~ NA_character_
    )
  )

group_colors <- c(
  "Outcomes"     = "#2a4d7c",
  "Municipality" = "#629460",
  "County"       = "#e2b13c"
)


# Plotting
plot_data_avail <- ggplot(df_plot, aes(x = factor(year), y = variable)) +
  geom_tile(aes(fill = group, alpha = available), colour = "white", linewidth = 0.4) +
  
  # Replace geom_text layers with a single geom_point layer mapping shapes
  geom_point(data = filter(df_plot, !is.na(status)),
             aes(shape = status), colour = "white", size = 5.5) +
  
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.07), guide = "none") +
  
  # Add override.aes to make the legend text/shapes visible
  scale_shape_manual(
    values = c("available" = "o", "imputed" = "~"), 
    name = NULL,
    guide = guide_legend(override.aes = list(colour = "grey30", size = 4))
  ) +
  
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_text(size = 11, colour = "grey10"),
    axis.text.x     = element_text(size = 11, colour = "grey10"),
    legend.position = "bottom",
    legend.box      = "horizontal", 
    legend.text     = element_text(size = 11, colour = "grey10"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# ggsave(filename = "_results//descriptive_results/plot_data_avail.png", plot = plot_data_avail, width = 12, height = 8, dpi = 300)


###################### Theory Plot ######################

# --- Generate data ---
x <- seq(0, 12, length.out = 800)

# --- OPTIMIZED PARAMETERS ---
main_amp <- 0.15
main_center <- 4.5
left_width <- 1.4   # steeper drop before dip
right_width <- 2.2  # slower recovery after dip

y_main <- ifelse(
  x < main_center,
  0.5 - main_amp * exp(-((x - main_center)^2) / left_width),
  0.5 - main_amp * exp(-((x - main_center)^2) / right_width)
)

# Green curve (Benefits > Costs)
green_amp <- 0.07
green_width <- 2
green_center <- 4.5
green_exp_amp <- 0.15
green_exp_rate <- 0.3
y_green_shape <- 0.5 - green_amp * exp(-((x - green_center)^2) / green_width)
y_green <- y_green_shape + ifelse(x >= green_center, green_exp_amp * (1 - exp(-green_exp_rate * (x - green_center))), 0)
y_green[x < 2.2] <- NA

# Red curve (Costs > Benefits)
red_amp <- 0.17
red_exp_rate <- 0.25
dip_index <- which.min(y_main)
dip_x <- x[dip_index]
y_red <- ifelse(x >= dip_x,
                y_main - red_amp * (1 - exp(-red_exp_rate * (x - dip_x))),
                NA)

# --- Combine into data frame ---
df <- data.frame(x, y_main, y_green, y_red)

# --- Plot Parameters ---
y_offset <- 0.03
x_line1 <- 3.8
x_line2 <- 6.6

# --- Plot ---
ggplot(df, aes(x = x)) +
  geom_line(aes(y = y_main + y_offset), color = "black", linewidth = 1.3) +
  geom_line(aes(y = y_green + y_offset, color = "Benefits > Costs", linetype = "Benefits > Costs"), linewidth = 1.1) +
  geom_line(aes(y = y_red + y_offset, color = "Costs > Benefits", linetype = "Costs > Benefits"), linewidth = 1.1) +
  
  # Dotted line 1: Construction finishes
  geom_vline(xintercept = x_line1, linetype = "dotted", color = "gray30", linewidth = 1) +
  annotate("text", x = x_line1, y = 0.72 + y_offset, label = "Construction", angle = 90, vjust = -0.5, size = 7) +
  
  # Dotted line 2: Habituation
  geom_vline(xintercept = x_line2, linetype = "dotted", color = "gray30", linewidth = 1) +
  annotate("text", x = x_line2, y = 0.72 + y_offset, label = "Habituation", angle = 90, vjust = -0.5, size = 7) +
  
  scale_color_manual(
    name = NULL,
    values = c("Benefits > Costs" = "darkgreen", "Costs > Benefits" = "red")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Benefits > Costs" = "dashed", "Costs > Benefits" = "dashed")
  ) +
  labs(x = "Time", y = "Level of Acceptance") +
  coord_cartesian(ylim = c(0.1, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", size = 23),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.0, "cm"),
    legend.text = element_text(size = 18)
  )

