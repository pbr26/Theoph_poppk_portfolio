# =============================================================================
# 08_two_compartment.R
# Project : Population PK Portfolio — Theophylline
# Phase 9 : One- vs two-compartment structural comparison
# Purpose : Fit a two-compartment oral model (Ka, CL, Vc, Q, V2) and test
#           whether the extra peripheral compartment is justified over the
#           selected one-compartment IIV model.
# Usage   : source(here::here("R", "08_two_compartment.R"))
# =============================================================================

source(here::here("R", "04_iiv_model.R"))   # build_nonmem_data(), model_metrics(),
                                             # compare_models(), fit_iiv_model(), setup

# -----------------------------------------------------------------------------
# 1. two_cmt_model(): two-compartment oral, first-order absorption
# -----------------------------------------------------------------------------
#' Central + peripheral compartments. IIV is placed on CL, Vc, and Ka (mirroring
#' the 1-cmt IIV model) for a fair comparison; Q and V2 are fixed effects only,
#' which aids convergence on a small dataset.
two_cmt_model <- function() {
  ini({
    tka    <- log(1.5)   # absorption rate (1/h)
    tcl    <- log(2.8)   # clearance (L/h)
    tvc    <- log(32)    # central volume (L)
    tq     <- log(2)     # inter-compartmental clearance (L/h)
    tvp    <- log(10)    # peripheral volume (L)
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.vc ~ 0.1
    add.sd  <- 0.2
    prop.sd <- 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    vc <- exp(tvc + eta.vc)
    q  <- exp(tq)
    vp <- exp(tvp)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / vc) * center -
                     (q / vc) * center + (q / vp) * periph
    d/dt(periph) <-  (q / vc) * center - (q / vp) * periph
    cp <- center / vc
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

fit_two_cmt_model <- function(dat = build_nonmem_data()) {
  nlmixr(
    two_cmt_model, data = as.data.frame(dat), est = "focei",
    control = foceiControl(print = 0)
  )
}

# -----------------------------------------------------------------------------
# 2. load_iiv_fit(): the 1-compartment reference (rebuild if missing)
# -----------------------------------------------------------------------------
load_iiv_fit <- function() {
  p <- fs::path(PATH_MODELS, "iiv_fit.rds")
  if (fs::file_exists(p)) readRDS(p) else fit_iiv_model()
}

# -----------------------------------------------------------------------------
# 3. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- build_nonmem_data()

  iiv_fit  <- load_iiv_fit()
  two_fit  <- fit_two_cmt_model(dat)
  saveRDS(two_fit, fs::path(PATH_MODELS, "two_cmt_fit.rds"))

  cmp <- compare_models(
    model_metrics(iiv_fit, "1-compartment (selected)"),
    model_metrics(two_fit, "2-compartment")
  )
  readr::write_csv(cmp, fs::path(PATH_TABLES, "model_comparison_1cmt_2cmt.csv"))

  print(cmp)
  message("\n✓ Phase 9 two-compartment fit complete. dOFV vs 1-cmt = ",
          cmp$dOFV[2], " | fit saved to outputs/models/two_cmt_fit.rds")
}
