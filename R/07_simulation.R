# =============================================================================
# 07_simulation.R
# Project : Population PK Portfolio — Theophylline
# Phase 7 : Clinical simulation
# Purpose : Use the qualified model to simulate (1) a typical patient,
#           (2) weight subgroups, and (3) a population with IIV, via rxode2.
# Usage   : source(here::here("R", "07_simulation.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
stopifnot(requireNamespace("rxode2", quietly = TRUE))
library(rxode2)

SIM_DOSE <- 300   # fixed oral dose (mg) for the simulations
SIM_TMAX <- 24    # hours
REF_WT   <- 70    # allometric reference weight (kg)

# -----------------------------------------------------------------------------
# 1. get_typical_params(): pull typical estimates + omega from the final fit
# -----------------------------------------------------------------------------
get_typical_params <- function() {
  for (f in c("cov_fit.rds", "iiv_fit.rds", "base_fit.rds")) {
    p <- fs::path(PATH_MODELS, f)
    if (fs::file_exists(p)) { fit <- readRDS(p); src <- f; break }
  }
  if (!exists("fit")) stop("No fitted model found; run Phases 3-5 first.")

  th <- tryCatch(fit$theta, error = function(e) NULL)
  if (is.null(th)) th <- nlmixr2::fixef(fit)
  getp <- function(nm, default) if (nm %in% names(th)) as.numeric(th[[nm]]) else default

  om <- tryCatch(fit$omega, error = function(e) NULL)

  list(
    source = src,
    ka  = exp(getp("tka", log(1.5))),
    cl  = exp(getp("tcl", log(2.8))),
    v   = exp(getp("tv",  log(32))),
    dWTCL = getp("dWTCL", 0.75),
    dWTV  = getp("dWTV",  1.0),
    omega = om
  )
}

# -----------------------------------------------------------------------------
# 2. rxode2 models
# -----------------------------------------------------------------------------
# Typical / subgroup model: parameters supplied directly (no random effects).
.mod_typ <- rxode2({
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
})

# Population model: random effects multiply the typical parameters.
.mod_iiv <- rxode2({
  ka <- ka_typ * exp(eta.ka)
  cl <- cl_typ * (WT / 70)^dWTCL * exp(eta.cl)
  v  <- v_typ  * (WT / 70)^dWTV  * exp(eta.v)
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
})

.ev <- function(dose = SIM_DOSE, tmax = SIM_TMAX) {
  rxode2::et(amt = dose, cmt = "depot") |>
    rxode2::et(seq(0, tmax, by = 0.1))
}

# -----------------------------------------------------------------------------
# 3. Simulations
# -----------------------------------------------------------------------------
#' Typical patient at a single reference weight (no IIV).
sim_typical <- function(pars, wt = REF_WT, dose = SIM_DOSE) {
  cl <- pars$cl * (wt / 70)^pars$dWTCL
  v  <- pars$v  * (wt / 70)^pars$dWTV
  s  <- rxode2::rxSolve(.mod_typ, params = c(ka = pars$ka, cl = cl, v = v),
                        events = .ev(dose))
  tibble::tibble(time = s$time, cp = s$cp, WT = wt)
}

#' Typical curve across weight subgroups (same dose).
sim_weight_subgroups <- function(pars, weights = c(55, 70, 85), dose = SIM_DOSE) {
  purrr::map_dfr(weights, function(w) sim_typical(pars, wt = w, dose = dose))
}

#' Population simulation with IIV: median + 90% prediction interval.
sim_population <- function(pars, wt = REF_WT, dose = SIM_DOSE, nsub = 500) {
  if (is.null(pars$omega)) {
    message("• No omega available — population simulation skipped.")
    return(NULL)
  }
  th <- c(ka_typ = pars$ka, cl_typ = pars$cl, v_typ = pars$v,
          dWTCL = pars$dWTCL, dWTV = pars$dWTV, WT = wt)
  set.seed(20260627)
  s <- rxode2::rxSolve(.mod_iiv, params = th, events = .ev(dose),
                       omega = pars$omega, nSub = nsub)
  tibble::as_tibble(s) |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      median = stats::median(cp),
      lo     = stats::quantile(cp, 0.05),
      hi     = stats::quantile(cp, 0.95),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# 4. Plot builders
# -----------------------------------------------------------------------------
sim_plot_typical <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(time, cp)) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Typical-patient simulation",
      subtitle = glue::glue("70 kg subject, single {SIM_DOSE} mg oral dose")
    )
}

sim_plot_weight <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(time, cp, colour = factor(WT))) +
    ggplot2::geom_line(linewidth = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo",
                                  name = "Weight (kg)") +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Effect of body weight on exposure",
      subtitle = glue::glue("Typical patients, same {SIM_DOSE} mg dose")
    )
}

sim_plot_population <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(time, median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "#3498DB", alpha = 0.25) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Population simulation with inter-individual variability",
      subtitle = "Median and 90% prediction interval (500 virtual subjects, 70 kg)"
    )
}

# -----------------------------------------------------------------------------
# 5. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  pars <- get_typical_params()
  message("Simulating from: ", pars$source)

  typ <- sim_typical(pars)
  wts <- sim_weight_subgroups(pars)
  pop <- sim_population(pars)

  readr::write_csv(wts, fs::path(PATH_DATA_PROC, "sim_weight_subgroups.csv"))
  save_fig(sim_plot_typical(typ), "07_sim_typical")
  save_fig(sim_plot_weight(wts),  "07_sim_weight")
  if (!is.null(pop)) {
    readr::write_csv(pop, fs::path(PATH_DATA_PROC, "sim_population_pi.csv"))
    save_fig(sim_plot_population(pop), "07_sim_population")
  }

  cmax_typ <- round(max(typ$cp), 2)
  message("\n✓ Phase 7 simulations complete. Typical Cmax = ", cmax_typ,
          " mg/L | figures in outputs/figures/.")
}
