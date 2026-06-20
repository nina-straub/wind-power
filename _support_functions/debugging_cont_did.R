################ Debug script contdid application ################
# When trying to implement analysis with contdid, encountered multiple errors
# This script documents these errors and their root causes
# Insight: Some parts of the package (level + dose and slope + dose) cannot be used bc. of these bugs

#### Clean df ####
# Make sure the df is clean to prevent potential bugs due to NA etc.
df_clean <- df_did_ready %>%
  group_by(ags) %>%
  mutate(
    # Take the (first observed) post-treatment dose for each unit and carry it
    # across all periods, including pre-treatment ones.
    treat_dose = if (any(seq_group > 0)) {
      max(wind_count_3km[seq_time == seq_group], na.rm = TRUE)
    } else {
      wind_count_3km
    }
  ) %>%
  ungroup() %>%
  filter(
    !is.na(df_did_ready$turnout),
    !is.na(treat_dose),
    !is.na(ags),
    !is.na(seq_time),
    !is.na(pop_density),
    election_year > 1994,                    # drop periods without within-period-dose variation
    !seq_group == 2                          # drop group without within-period-dose variation
  ) %>%
  as.data.frame()                           # turn into data.frame

#### Balance panel ####
expected_periods <- n_distinct(df_clean$seq_time)

df_balanced <- df_clean %>%
  group_by(ags) %>%
  filter(n() == expected_periods) %>%
  ungroup()

message(paste0(
  "  Balanced panel: ", n_distinct(df_balanced$ags), " units × ",
  expected_periods, " periods (", nrow(df_balanced), " obs)"
))


# Check: Are there any groups without variation in within-period-dosage? --> Nope
df_balanced %>%
  filter(seq_group > 0) %>%
  group_by(seq_group) %>%
  summarise(
    unique_doses = n_distinct(treat_dose),
    min_dose = min(treat_dose),
    max_dose = max(treat_dose)
  )


#### Run different estimation specifications ####
# Check which specifications work in the clean, balance setup

# 1. slope - eventstudy
res1 <- cont_did(
  yname           = "turnout",
  tname           = "seq_time",
  idname          = "ags",
  dname           = "treat_dose",
  gname           = "seq_group",
  data            = df_balanced,
  target_parameter = "slope",
  aggregation     = "eventstudy",
  treatment_type  = "continuous",
  control_group   = "notyettreated",
  biters          = 1000,
  cband           = TRUE,
  num_knots       = 2,
  degree          = 5
)

summary(res1)
ggcont_did(res1, type = "att")

# Runs without issues (yeay) and (sometimes) without confidence bands warnings when using non-linear spline


#### 2. level - dose ####
res2 <- cont_did(
  yname           = "turnout",
  tname           = "seq_time",
  idname          = "ags",
  dname           = "treat_dose",
  gname           = "seq_group",
  data            = df_balanced,
  target_parameter = "level",
  aggregation     = "dose",
  treatment_type  = "continuous",
  control_group   = "notyettreated",
  biters          = 100,
  cband           = TRUE,
  num_knots       = 0,
  degree          = 1
)

# "Error in overall_weights(att_gt, ...): something's going wrong calculating overall weights"

# Simple sanity checks:
# Are never-treated units properly coded as 0? --> Yes
# Any cohort with only one unit? --> No
# Treat_dose truly time-invariant per unit? --> Yes
# Also: "seq(min(dose[dose > 0]), max(dose[dose > 0]), length.out = 50)" from npiv handles skewed distribution

# Look at source code:
if (aggregation == "dose" && target_parameter == "level") {
  attgt_fun <- cont_did_acrt  # computes both att and acrt
  gt_type <- "dose"           # <-- dose gt_type
}
if (aggregation == "eventstudy" && target_parameter == "slope") {
  attgt_fun <- cont_did_acrt
  gt_type <- "att"            # <-- att gt_type
}

# Since 1.) slope + eventstudy works, but 2.) dose + level fails, sth. withing gt_type = dose causes issues
# Debug report: https://github.com/bcallaway11/contdid/issues/6
# --> Floating point issue when small values > 1 cause if (sum(out_weight) != 1) to break

# No chance to fix this --> Can't use level - dose setup


#### 3. level - eventstudy ####
res3 <- cont_did(
  yname           = "turnout",
  tname           = "seq_time",
  idname          = "ags",
  dname           = "treat_dose_squish",
  gname           = "seq_group",
  data            = df_balanced,
  target_parameter = "level",
  aggregation     = "eventstudy",
  treatment_type  = "continuous",
  control_group   = "nevertreated",
  biters          = 1000,
  cband           = TRUE,
  num_knots       = 0,
  degree          = 1
)

# "Error in eval(family$initialize) : y values must be 0 <= y <= 1"

# Sanity checks on outcome turnout
# Something NA? --> Nope
any(is.infinite(df_balanced$turnout))
# Problems with turnout above 1? Switching it to percent-scaling? --> doesn't help
df_balanced <-  df_balanced %>% mutate(turnout = pmin(0.99999, pmax(0.00001, turnout)))
df_balanced <-  df_balanced %>% mutate(turnout = 100*pmin(0.99999, pmax(0.00001, turnout)))
# tbl_df factor? --> restructured to data.frame
# Some id_name multiple times? --> Nope
table(table(df_balanced$ags, df_balanced$seq_time))

# Error also emerges when trying different control_group or discarding bootstrapping

# Diagnostic attempts:
# What is y? --> Error stems from glm() function, where is this used?
# Error thrown by ptetools::did_attgt using DR and calls glm() for propensity scores
# I.e. y is the treatment indicator

# Check what did_attgt sees as its "treatment"
df_balanced %>%
  mutate(treated = as.numeric(seq_group > 0)) %>%
  summarise(
    min = min(treated),
    max = max(treated),
    any_na = any(is.na(treated)),
    class = class(treated)
  )
# Check seq_group coding
table(df_balanced$seq_group)
# --> Looks fine

# Check what did_attgt does
getFromNamespace("did_attgt", "ptetools")
# --> Problem. Internal DRDID::drdid_panel() uses D as treatment indicator, which is continuous dose variable
# BUT: It expects binary variable

# Proof of diagnosis: Squish dose to (0,1) fixes issue
df_balanced <- df_balanced %>% mutate (treat_dose_squish = treat_dose / (1 + treat_dose))

summary(res3)
ggcont_did(res3, type = "att")

# --> Error wasn't detected bc. "simulate_contdid_data" only simulates doses (0, 1)
# Bug was reported: https://github.com/bcallaway11/contdid/issues/11


# 4. slope - dose
res4 <- cont_did(
  yname           = "turnout",
  tname           = "seq_time",
  idname          = "ags",
  dname           = "treat_dose",
  gname           = "seq_group",
  data            = df_balanced,
  target_parameter = "slope",
  aggregation     = "dose",
  treatment_type  = "continuous",
  control_group   = "notyettreated",
  biters          = 1000,
  cband           = TRUE,
  num_knots       = 2,
  degree          = 5
)

# Same floating point issue as case 2.)



#### Miscellaneous ####
# Typo in line 348/350 of contdid/R/cont_did.R (https://github.com/bcallaway11/contdid/blob/main/R/cont_did.R):
# dose_est_method = dose_est_method is set twice

# Typo in line 162 of pte/R/process_att_gt.R (https://github.com/bcallaway11/pte/blob/master/R/process_att_gt.R):
# "using pointwise confidence interal" --> "using pointwise confidence interval"


#### Minimal Reproducible Example for "level + eventstudy" ####
# Change D <- runif(n, 0, 1) to D <- runif(n, 0, 1.1)
sim_contdid_dat_changed <- function(
    n = 5000,
    num_time_periods = 4,
    num_groups = num_time_periods,
    pg = rep(1 / num_groups, num_groups - 1),
    pu = 1 / (num_groups),
    dose_linear_effect = 0,
    dose_quadratic_effect = 0) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for this function but is not installed.
         Please install it with install.packages('tidyr').", call. = FALSE)
  }
  time_periods <- 1:num_time_periods
  groups <- c(0, time_periods[-1])
  p <- c(pu, pg)
  G <- sample(groups, n, replace = TRUE, prob = p)
  D <- runif(n, 0, 1.1)
  eta <- rnorm(n, mean = G)
  time_effects <- 1:num_time_periods
  Y0t <- sapply(1:num_time_periods, function(tp) {
    time_effects[tp] + eta + rnorm(n)
  })
  Y1t <- sapply(1:num_time_periods, function(tp) {
    dose_linear_effect * D + dose_quadratic_effect * D^2 + time_effects[tp] + eta + rnorm(n)
  })
  post_mat <- sapply(1:num_time_periods, function(tp) {
    1 * ((G <= tp) & G != 0)
  })
  Y <- post_mat * Y1t + (1 - post_mat) * Y0t
  df <- as.data.frame(Y)
  colnames(df) <- paste0("Y_", 1:num_time_periods)
  df$id <- 1:n
  df$G <- G
  df$D <- D
  df2 <- tidyr::pivot_longer(df,
                             cols = tidyr::starts_with("Y"),
                             names_to = "time_period",
                             names_prefix = "Y_",
                             names_transform = list(time_period = as.numeric),
                             values_to = "Y"
  ) |> as.data.frame()
  df2$D[df2$G == 0] <- 0
  df2
}


set.seed(1234)
# Simulate data (same setup as in example)
df <- sim_contdid_dat_changed(
  n = 5000,
  num_time_periods = 4,
  num_groups = 4,
  dose_linear_effect = 0,
  dose_quadratic_effect = 0
)

# Try to estimate level - eventstudy setting (same as in example)
cd_res <- cont_did(
  yname = "Y",
  tname = "time_period",
  idname = "id",
  dname = "D",
  data = df,
  gname = "G",
  target_parameter = "level",
  aggregation = "eventstudy",
  treatment_type = "continuous",
  control_group = "notyettreated",
  biters = 100,
  cband = TRUE,
  num_knots = 1,
  degree = 3,
)

# Gives "Error in eval(family$initialize) : y values must be 0 <= y <= 1"
# Setup level + eventstudy only constellation calling "ptetools::did_attgt"
# ptetools::did_attgt seems to require D = [0,1] for propensity score estimation

# Test hypothesis: squish D ~ [0, 1.1] to D ~ [0, 1]
df <- df %>% mutate (D_bin = D / (1 + D))

# Re-run cont_did with level + eventstudy
cd_res <- cont_did(
  yname = "Y",
  tname = "time_period",
  idname = "id",
  dname = "D_bin",
  data = df,
  gname = "G",
  target_parameter = "level",
  aggregation = "dose",
  treatment_type = "continuous",
  control_group = "notyettreated",
  biters = 100,
  cband = TRUE,
  num_knots = 1,
  degree = 3,
)

# --> No error, runs smoothly


