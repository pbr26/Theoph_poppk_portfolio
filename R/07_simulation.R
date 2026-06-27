# =============================================================================
# 07_simulation.R
# Project : Population PK Portfolio — Theophylline
# Phase 7 : Clinical simulation
# Purpose : Use the SELECTED model (the IIV model — Phase 5 found weight was not
#           statistically supported) to simulate (1) a typical patient,
#           (2) weight subgroups via standard allometric scaling (illustrative),
#           and (3) a population with IIV, all via rxode2.
# Usage   : source(here::here("R", "07_simulation.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
stopifnot(requireNamespace("rxode2", quietly = TRUE))
library(rxode2)

SIM_TMAX <- 24    # hours
REF_WT   <- 70    # allometric reference weight (kg)
ALLO_CL  <- 0.75  # theoretical allometric exponent for clearance
ALLO_V   <- 1.0   # theoretical allometric exponent for volume

# -----------------------------------------------------------------------------
# 1. rep_dose(): a representative (median) actual dose from the dataset
# -----------------------------------------------------------------------------
rep_dose <- function() {
  p <- fs::path(PATH_DATA_PROC, "theoph_clean.csv")
  if (!fs::file_exists(p)) return(320)
  df <- readr::read_csv(p, show_col_types = FALSE)
  round(stats::median(dplyr::distinct(df, id, amt)$amt))
}

# -----------------------------------------------------------------------------
# 2. get_typical_params(): typical estimates + omega from the SELECTED (IIV) fit
# -----------------------------------------------------------------------------
get_typical_params <- function() {
  fit <- NULL; src <- NA_character_
  for (f in c("iiv_fit.rds", "cov_fit.rds", "base_fit.rds")) {
    p <- fs::path(PATH_MODELS, f)
    if (fs::file_exists(p)) { fit <- readRDS(p); src <- f; break }
  }
  if (is.null(fit)) stop("No fitted model found; run Phases 3-4 first.")

  th <- tryCatch(fit$theta, error = function(e) nlmixr2::fixef(fit))
  gp <- function(nm, default) if (nm %in% names(th)) as.numeric(th[[nm]]) else default

  list(
    source = src,
    ka  = exp(gp("tka", log(1.5))),
    cl  = exp(gp("tcl", log(2.8))),
    v   = exp(gp("tv",  log(32))),
    omega = tryCatch(fit$omega, error = function(e) NULL)
  )
}

# -----------------------------------------------------------------------------
# 3. rxode2 models
# -----------------------------------------------------------------------------
.mod_typ <- rxode2({
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
})

.mod_iiv <- rxode2({
  ka <- ka_typ * exp(eta.ka)
  cl <- cl_typ * exp(eta.cl)
  v  <- v_typ  * exp(eta.v)
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
})

.ev <- function(dose, tmax = SIM_TMAX) {
  rxode2::et(amt = dose, cmt = "depot") |>
    rxode2::et(seq(0, tmax, by = 0.1))
}

# -----------------------------------------------------------------------------
# 4. Simulations
# -----------------------------------------------------------------------------
#' Typical patient. Weight enters only through THEORETICAL allometric scaling
#' (illustrative) — the selected model itself has no weight covariate.
sim_typical <- function(pars, wt = REF_WT, dose = rep_dose()) {
  cl <- pars$cl * (wt / REF_WT)^ALLO_CL
  v  <- pars$v  * (wt / REF_WT)^ALLO_V
  s  <- rxode2::rxSolve(.mod_typ, params = c(ka = pars$ka, cl = cl, v = v),
                        events = .ev(dose))
  tibble::tibble(time = s$time, cp = s$cp, WT = wt)
}

#' Illustrative weight subgroups via standard allometric exponents.
sim_weight_subgroups <- function(pars, weights = c(55, 70, 85), dose = rep_dose()) {
  purrr::map_dfr(weights, function(w) sim_typical(pars, wt = w, dose = dose))
}

#' Population simulation with IIV (at the reference weight): median + 90% PI.
sim_population <- function(pars, wt = REF_WT, dose = rep_dose(), nsub = 500) {
  if (is.null(pars$omega)) {
    message("• No omega available — population simulation skipped.")
    return(NULL)
  }
  th <- c(ka_typ = pars$ka, cl_typ = pars$cl, v_typ = pars$v)
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
# 5. Plot builders
# -----------------------------------------------------------------------------
sim_plot_typical <- function(df, dose = rep_dose()) {
  ggplot2::ggplot(df, ggplot2::aes(time, cp)) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Typical-patient simulation",
      subtitle = glue::glue("70 kg subject, single {dose} mg oral dose")
    )
}

sim_plot_weight <- function(df, dose = rep_dose()) {
  ggplot2::ggplot(df, ggplot2::aes(time, cp, colour = factor(WT))) +
    ggplot2::geom_line(linewidth = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo",
                                  name = "Weight (kg)") +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Illustrative effect of body weight (allometric scaling)",
      subtitle = glue::glue("Same {dose} mg dose; theoretical exponents 0.75 (CL), 1.0 (V)")
    )
}

sim_plot_population <- function(df, dose = rep_dose()) {
  ggplot2::ggplot(df, ggplot2::aes(time, median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "#3498DB", alpha = 0.25) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = "Population simulation with inter-individual variability",
      subtitle = glue::glue("Median and 90% PI (500 virtual subjects, 70 kg, {dose} mg)")
    )
}

# -----------------------------------------------------------------------------
# 6. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  pars <- get_typical_params()
  dose <- rep_dose()
  message("Simulating from: ", pars$source, " | representative dose = ", dose, " mg")

  typ <- sim_typical(pars, dose = dose)
  wts <- sim_weight_subgroups(pars, dose = dose)
  pop <- sim_population(pars, dose = dose)

  readr::write_csv(wts, fs::path(PATH_DATA_PROC, "sim_weight_subgroups.csv"))
  save_fig(sim_plot_typical(typ, dose), "07_sim_typical")
  save_fig(sim_plot_weight(wts, dose),  "07_sim_weight")
  if (!is.null(pop)) {
    readr::write_csv(pop, fs::path(PATH_DATA_PROC, "sim_population_pi.csv"))
    save_fig(sim_plot_population(pop, dose), "07_sim_population")
  }

  cmax_typ <- round(max(typ$cp), 2)
  message("\n✓ Phase 7 simulations complete. Typical Cmax = ", cmax_typ,
          " mg/L | figures in outputs/figures/.")
}
