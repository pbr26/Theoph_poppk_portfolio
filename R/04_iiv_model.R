# =============================================================================
# 04_iiv_model.R
# Project : Population PK Portfolio — Theophylline
# Phase 4 : Inter-individual variability (IIV) model
# Purpose : Extend the base model with between-subject random effects on
#           CL, V, and Ka; refit; compare to the base model; inspect the
#           eta distributions.
# Usage   : source(here::here("R", "04_iiv_model.R"))   # defines functions
#           Sourcing directly runs the pipeline (demo block at the bottom).
# =============================================================================

source(here::here("R", "03_base_model.R"))   # brings build_nonmem_data(), setup

# -----------------------------------------------------------------------------
# 1. iiv_model(): structural model + IIV on CL, V, Ka (diagonal omega)
# -----------------------------------------------------------------------------
#' Each parameter gets its own random effect (eta) drawn from N(0, omega^2).
#' Starting omegas are modest variances; fixed-effect inits reuse the base.
iiv_model <- function() {
  ini({
    tka    <- log(1.5)
    tcl    <- log(2.8)
    tv     <- log(32)
    eta.ka ~ 0.3      # IIV on absorption rate
    eta.cl ~ 0.1      # IIV on clearance
    eta.v  ~ 0.1      # IIV on volume
    add.sd  <- 0.2
    prop.sd <- 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

fit_iiv_model <- function(dat = build_nonmem_data()) {
  nlmixr(
    iiv_model, data = as.data.frame(dat), est = "focei",
    control = foceiControl(print = 0)
  )
}

# -----------------------------------------------------------------------------
# 2. model_metrics(): pull OFV / AIC / BIC / #params from a fit
# -----------------------------------------------------------------------------
#' Robustly extract comparison metrics from an nlmixr2 fit.
model_metrics <- function(fit, name) {
  od <- as.data.frame(fit$objDf)
  pick <- function(cols) {
    hit <- intersect(cols, names(od))
    if (length(hit)) suppressWarnings(as.numeric(od[[hit[1]]][1])) else NA_real_
  }
  tibble::tibble(
    Model = name,
    OFV   = pick(c("OBJF", "objf", "OFV")),
    AIC   = pick(c("AIC")),
    BIC   = pick(c("BIC")),
    nPar  = tryCatch(length(fit$ini$est[!is.na(fit$ini$est)]), error = function(e) NA_integer_)
  )
}

#' Build a base-vs-IIV comparison table (adds dOFV vs the first row).
compare_models <- function(...) {
  tbl <- dplyr::bind_rows(...)
  tbl |>
    dplyr::mutate(
      dOFV = round(OFV - OFV[1], 2),
      dAIC = round(AIC - AIC[1], 2)
    )
}

# -----------------------------------------------------------------------------
# 3. iiv_etas() + iiv_plot_eta(): extract and visualise random effects
# -----------------------------------------------------------------------------
iiv_etas <- function(fit) {
  e <- as.data.frame(fit$eta)
  id_col <- intersect(c("ID", "id"), names(e))
  eta_cols <- grep("^eta", names(e), value = TRUE)
  e |>
    dplyr::select(dplyr::all_of(c(id_col, eta_cols))) |>
    tibble::as_tibble()
}

#' Boxplot + jitter of each eta; should be centred on 0 and roughly symmetric.
iiv_plot_eta <- function(fit) {
  e <- iiv_etas(fit)
  id_col <- intersect(c("ID", "id"), names(e))
  long <- e |>
    tidyr::pivot_longer(-dplyr::all_of(id_col), names_to = "eta", values_to = "value")
  ggplot2::ggplot(long, ggplot2::aes(eta, value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#3498DB", alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.75, colour = "#2C3E50") +
    ggplot2::labs(
      x = NULL, y = "Random effect (eta)",
      title = "Empirical-Bayes eta distributions",
      subtitle = "Each eta should be centred on 0 and roughly symmetric"
    )
}

# -----------------------------------------------------------------------------
# 4. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- build_nonmem_data()

  base_fit <- readRDS(fs::path(PATH_MODELS, "base_fit.rds"))
  iiv_fit  <- fit_iiv_model(dat)
  saveRDS(iiv_fit, fs::path(PATH_MODELS, "iiv_fit.rds"))

  cmp <- compare_models(
    model_metrics(base_fit, "Base (IIV on CL)"),
    model_metrics(iiv_fit,  "IIV on CL, V, Ka")
  )
  readr::write_csv(cmp, fs::path(PATH_TABLES, "model_comparison_base_iiv.csv"))

  save_fig(iiv_plot_eta(iiv_fit), "04_iiv_eta_distributions")

  print(cmp)
  message("\n✓ Phase 4 IIV model fitted. dOFV vs base = ",
          cmp$dOFV[2], " | fit saved to outputs/models/iiv_fit.rds")
}
