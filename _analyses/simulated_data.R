###################### Simulation of Data Following the Implementation Dip ######################

# Two DGPs to check whether CS and dCDH can recover an implementation dip.
# Sim A: Staggered and non-binary, but absorbing
# Sim B: Additionally non-absorbing, i.e. units can be treated multiple times


###################### Parameters ######################

#### Panel size ####
N  <- 10000    # units
TT <- 10       # periods


#### Treatment assignment ####
COHORTS     <- 3:9                                   # possible first-treatment periods
COHORT_W    <- c(.16, .18, .18, .16, .14, .10, .08)  # tilted early, so more cohorts
COHORT_W    <- COHORT_W / sum(COHORT_W)              # observe the full recovery
SELECTION   <- 4.0    # selection on unit FE level, 0 = none. Does not break PT.
MAX_CUM     <- 150    # cap on cumulative treatment
MAX_PER_PER <- 50     # cap on treatment added in a single period
DELTA_P     <- 1/12   # geometric parameter for treatment intensity (mean ~12)
P_EVENT     <- 0.45   # Sim B only, per-period prob. of receiving another treatment


#### Treatment effect ####
# Dose-response is concave in the cumulative count. An event's magnitude is proportional to the change in log(1 + D), i.e. the first turbine matters most
KAPPA <- 0.009

# Event-time path of a single treatment event, i.e. the implementation dip.
KERNEL <- c(1.0, 1.25, 0.9, 0.4)  # e = 0, 1, 2, 3
KMAX   <- length(KERNEL) - 1

# Effect is negative
DIP_SIGN <- -1

# Calendar gradient, an event in period s is scaled by 1 + RHO*(s-1)/(TT-1) --> later events have larger effects
RHO      <- 0.6

# Unit-level effect heterogeneity (lognormal, mean 1)
SD_ALPHA <- 0.4


#### Outcome ####
MU       <- 0.30    # baseline
SD_UNIT  <- 0.05    # unit fixed effects
SD_NOISE <- 0.03    # noise
DELTA_T  <- c(0, .01, .015, .005, -.005, .01, .02, .015, .025, .03)  # time FE

stopifnot(length(DELTA_T) == TT)


###################### Simulation Function ######################

simulate_panel <- function(sim = c("A", "B")) {
  
  sim <- match.arg(sim)
  
  #### Unit heterogeneity ####
  u     <- rnorm(N, 0, SD_UNIT)
  alpha <- rlnorm(N, meanlog = -0.5 * SD_ALPHA^2, sdlog = SD_ALPHA)  # E[alpha] = 1
  
  #### Who is treated, and when ####
  # Treatment depends on the level of the unit FE, i.e. selection on levels, keeping PTA intact, since u is time-invariant and time effects common
  p_treat <- plogis(SELECTION * u)
  ever    <- runif(N) < p_treat
  
  G <- integer(N)   # 0 = never treated
  G[ever] <- sample(COHORTS, sum(ever), replace = TRUE, prob = COHORT_W)
  
  #### Treatment path ####
  # Cumulative, weakly increasing and capped. Event indicator matrix first.
  E <- matrix(FALSE, N, TT)
  for (t in seq_len(TT)) {
    first <- ever & (G == t)
    E[first, t] <- TRUE
    if (sim == "B") {
      later <- ever & (G < t)
      E[later, t] <- runif(sum(later)) < P_EVENT
    }
  }
  
  # Increment sizes, capped per period
  Draw <- matrix(pmin(1L + rgeom(N * TT, DELTA_P), MAX_PER_PER), N, TT)
  
  # Accumulate, respecting the overall cap
  D   <- matrix(0, N, TT)
  cum <- numeric(N)
  for (t in seq_len(TT)) {
    add <- ifelse(E[, t], pmin(Draw[, t], MAX_CUM - cum), 0)
    add <- pmax(add, 0)
    cum <- cum + add
    D[, t] <- cum
  }
  
  #### Treatment effects ####
  # Change in the concave dose transform, per period
  Fd <- log1p(D)
  dF <- Fd - cbind(0, Fd[, -TT, drop = FALSE])
  
  # Every event launches its own path and the paths add up.
  # tau = all events layered --> what CS and dCDH target
  # tau_dip = first event only --> the implementation dip itself
  tau     <- matrix(0, N, TT)
  tau_dip <- matrix(0, N, TT)
  for (s in seq_len(TT)) {
    gam   <- 1 + RHO * (s - 1) / (TT - 1)
    mag   <- DIP_SIGN * alpha * KAPPA * gam * dF[, s]
    first <- (G == s)
    for (e in 0:KMAX) {
      tt <- s + e
      if (tt <= TT) {
        tau[, tt]     <- tau[, tt]     + mag * KERNEL[e + 1]
        tau_dip[, tt] <- tau_dip[, tt] + mag * KERNEL[e + 1] * first
      }
    }
  }
  
  #### Outcomes ####
  Y0 <- MU +
    matrix(u, N, TT) +
    matrix(DELTA_T, N, TT, byrow = TRUE) +
    matrix(rnorm(N * TT, 0, SD_NOISE), N, TT)
  Y  <- pmin(pmax(Y0 + tau, 0), 1)
  
  #### Long format ####
  df <- data.frame(
    unit     = rep(seq_len(N), times = TT),
    time     = rep(seq_len(TT), each  = N),
    G        = rep(G, times = TT),               # CS group, 0 = never treated
    D_cum    = as.vector(D),                     # cumulative treatment count
    treated  = as.integer(as.vector(D) > 0),     # absorbing 0/1 indicator
    Y        = as.vector(Y),                     # observed outcome
    Y0       = as.vector(Y0),                    # untreated potential outcome
    tau_true = as.vector(tau),                   # effect of all events
    tau_dip  = as.vector(tau_dip)                # effect of the first event only
  )
  df <- df[order(df$unit, df$time), ]
  rownames(df) <- NULL
  df
}


###################### Ground Truth ######################

# Dynamic effect relative to first treatment.
# var = "tau_true" --> the aggregate estimand CS and dCDH target
# var = "tau_dip"  --> the single-event implementation dip
true_event_study <- function(df, var = "tau_true", min_e = -4, max_e = 5) {
  d <- df[df$G > 0, ]
  d$e <- d$time - d$G
  d <- d[d$e >= min_e & d$e <= max_e, ]
  out <- aggregate(d[[var]], by = list(e = d$e), FUN = mean)
  names(out)[2] <- "att_true"
  out
}

# True ATT(g,t) for post-treatment cells, i.e. the CS estimand cell by cell
true_att_gt <- function(df) {
  d <- df[df$G > 0 & df$time >= df$G, ]
  out <- aggregate(tau_true ~ G + time, data = d, FUN = mean)
  names(out)[3] <- "att_true"
  out[order(out$G, out$time), ]
}


###################### Simulate and Describe ######################

#### Descriptives ####
describe <- function(df, label) {
  last <- df[df$time == TT, ]
  ntr  <- sum(last$G > 0)
  inc  <- ave(df$D_cum, df$unit, FUN = function(x) c(x[1], diff(x)))
  nev  <- tapply(inc > 0, df$unit, sum)
  cat("\n===== SIM", label, "=====\n")
  cat("Ever treated       :", ntr, sprintf("(%.1f%%)\n", 100 * ntr / N))
  cat("Cohort sizes       :"); print(table(last$G[last$G > 0]))
  cat("Final dose (treated): mean", round(mean(last$D_cum[last$G > 0]), 1),
      "| max", max(last$D_cum), "\n")
  cat("Max single-period increment:", max(inc), "\n")
  cat("Events per treated unit    : mean",
      round(mean(nev[last$G > 0]), 2), "| max", max(nev), "\n")
  cat("Outcome range      :", round(range(df$Y), 3), "\n")
}


#### Build both panels ####

set.seed(0309)

simA <- simulate_panel("A")
simB <- simulate_panel("B")

describe(simA, "A")
describe(simB, "B")


#### Check imposed truth ####
esA <- true_event_study(simA); esB <- true_event_study(simB)

# TRUE dynamic ATT relative to first treatment
merge(setNames(esA, c("e", "simA")), setNames(esB, c("e", "simB")), by = "e")

