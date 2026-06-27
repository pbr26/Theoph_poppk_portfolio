# =============================================================================
# 03_base_model.R
# Project : Population PK Portfolio — Theophylline
# Phase 3 : Base population PK model
# Purpose : Build the NONMEM-convention analysis dataset, define the base
#           one-compartment oral model (first-order absorption), fit it with
#           nlmixr2 (FOCEi), and persist the fit + parameter table.
# Usage   : source(here::here("R", "03_base_model.R"))   # defines functions
#           Sourcing directly runs the pipeline (demo block at the bottom).
# =============================================================================

source(here::here("R", "00_setup.R"))
stopifnot(requireNamespace("nlmixr2", quietly = TRUE))
library(nlmixr2)

# -----------------------------------------------------------------------------
# 1. load_pk_data(): the cleaned dataset from Phase 1 (rebuild if missing)
# -----------------------------------------------------------------------------
load_pk_data <- function() {
  clean_path <- fs::path(PATH_DATA_PROC, "theoph_clean.csv")
  if (fs::file_exists(clean_path)) {
    readr::read_csv(clean_path, show_col_types = FALSE)
  } else {
    source(here::here("R", "01_eda.R"))
    prepare_theoph()
  }
}

# -----------------------------------------------------------------------------
# 2. build_nonmem_data(): event-level (NONMEM-convention) dataset
# -----------------------------------------------------------------------------
#' Convert the tidy concentration table into an event dataset.
#'
#' Columns (NONMEM convention):
#'   ID    subject
#'   TIME  hours since dose
#'   DV    observed concentration (NA on dosing rows)
#'   AMT   dose amount in mg (only on the dosing row)
#'   EVID  1 = dose, 0 = observation
#'   CMT   1 = depot (dose), 2 = central (observation)
#'   MDV   1 = missing DV (dose rows), 0 = observation rows
#'   WT    body weight (kg), carried for the covariate model (Phase 5)
#'
#' One dosing row per subject at TIME 0, then all observation rows.
#' @return A tibble ordered by ID, TIME, with dose before obs at t = 0.
build_nonmem_data <- function(df = load_pk_data()) {
  obs <- df |>
    dplyr::transmute(
      ID = id, TIME = time, DV = conc,
      AMT = 0, EVID = 0L, CMT = 2L, MDV = 0L, WT = wt
    )

  dose <- df |>
    dplyr::distinct(id, amt, wt) |>
    dplyr::transmute(
      ID = id, TIME = 0, DV = NA_real_,
      AMT = amt, EVID = 1L, CMT = 1L, MDV = 1L, WT = wt
    )

  dplyr::bind_rows(dose, obs) |>
    # dose (EVID 1) must precede the observation at the same time
    dplyr::arrange(ID, TIME, dplyr::desc(EVID))
}

# -----------------------------------------------------------------------------
# 3. base_model(): structural 1-compartment oral model (IIV on CL)
# -----------------------------------------------------------------------------
#' nlmixr2 model function. Fixed effects are estimated on the log scale so the
#' back-transformed CL/V/Ka are guaranteed positive. Initial estimates are
#' seeded from the Phase 2 NCA (CL ~ 2.74 L/h, V ~ 31.6 L, Ka ~ 1.5 / h).
base_model <- function() {
  ini({
    tka    <- log(1.5)    # log typical absorption rate constant (1/h)
    tcl    <- log(2.8)    # log typical clearance (L/h)
    tv     <- log(32)     # log typical central volume (L)
    eta.cl ~ 0.1          # between-subject variability on CL (variance)
    add.sd  <- 0.3        # additive residual error (mg/L)
    prop.sd <- 0.1        # proportional residual error (fraction)
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# -----------------------------------------------------------------------------
# 4. fit_base_model(): estimate with FOCEi
# -----------------------------------------------------------------------------
fit_base_model <- function(dat = build_nonmem_data()) {
  # nlmixr2 exposes the estimation entry point as nlmixr() (attached by library)
  nlmixr(
    base_model, data = as.data.frame(dat), est = "focei",
    control = foceiControl(print = 0)
  )
}

# -----------------------------------------------------------------------------
# 5. base_param_table(): tidy back-transformed fixed-effects table
# -----------------------------------------------------------------------------
#' Build a readable parameter table from a fitted nlmixr2 object.
base_param_table <- function(fit) {
  pf <- as.data.frame(fit$parFixed)
  pf$Parameter <- rownames(pf)
  tibble::as_tibble(pf) |>
    dplyr::relocate(Parameter)
}

# -----------------------------------------------------------------------------
# 6. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- build_nonmem_data()
  readr::write_csv(dat, fs::path(PATH_DATA_PROC, "theoph_nonmem.csv"))

  fit <- fit_base_model(dat)
  saveRDS(fit, fs::path(PATH_MODELS, "base_fit.rds"))

  ofv <- tryCatch(as.numeric(fit$objDf$OBJF[1]), error = function(e) NA_real_)
  print(fit)
  message("\n✓ Phase 3 base model fitted",
          if (!is.na(ofv)) paste0(" (OFV = ", round(ofv, 2), ")") else "",
          ". Fit saved to outputs/models/base_fit.rds")
}
