# =============================================================================
# 05_covariate_model.R
# Project : Population PK Portfolio — Theophylline
# Phase 5 : Covariate model — body weight on CL and V (allometric)
# Purpose : Add an allometric weight effect (normalised to 70 kg) on CL and V,
#           estimate the exponents, compare to the IIV model, and check whether
#           the random IIV on CL/V shrinks.
# Usage   : source(here::here("R", "05_covariate_model.R"))
# =============================================================================

source(here::here("R", "04_iiv_model.R"))   # brings IIV model, metrics, setup

REF_WT <- 70   # reference body weight (kg) for allometric normalisation

# -----------------------------------------------------------------------------
# 1. cov_model(): allometric weight on CL and V (exponents estimated)
# -----------------------------------------------------------------------------
#' CL = exp(tcl + eta.cl) * (WT/REF_WT)^dWTCL
#' V  = exp(tv  + eta.v ) * (WT/REF_WT)^dWTV
#' Exponents start at allometric theory (0.75 for CL, 1.0 for V).
cov_model <- function() {
  ini({
    tka    <- log(1.5)
    tcl    <- log(2.8)
    tv     <- log(32)
    dWTCL  <- 0.75    # weight exponent on clearance
    dWTV   <- 1.0     # weight exponent on volume
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    add.sd  <- 0.2
    prop.sd <- 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl) * (WT / 70)^dWTCL
    v  <- exp(tv  + eta.v ) * (WT / 70)^dWTV
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

fit_cov_model <- function(dat = build_nonmem_data()) {
  nlmixr(
    cov_model, data = as.data.frame(dat), est = "focei",
    control = foceiControl(print = 0)
  )
}

# -----------------------------------------------------------------------------
# 2. load_iiv_fit(): reference (no-covariate) model
# -----------------------------------------------------------------------------
load_iiv_fit <- function() {
  p <- fs::path(PATH_MODELS, "iiv_fit.rds")
  if (fs::file_exists(p)) readRDS(p) else fit_iiv_model()
}

# -----------------------------------------------------------------------------
# 3. iiv_cv(): extract between-subject CV% on CL and V from a fit
# -----------------------------------------------------------------------------
#' Pull the BSV(CV%) for named etas from parFixed, robust to formatting.
iiv_cv <- function(fit) {
  pf <- as.data.frame(fit$parFixed)
  col <- grep("BSV", names(pf), value = TRUE)
  if (!length(col)) return(tibble::tibble(eta = character(), cv = numeric()))
  vals <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pf[[col[1]]])))
  tibble::tibble(Parameter = rownames(pf), `BSV CV%` = vals) |>
    dplyr::filter(!is.na(`BSV CV%`))
}

# -----------------------------------------------------------------------------
# 4. cov_plot_eta_wt(): residual random effects vs weight
# -----------------------------------------------------------------------------
#' After including weight, eta.cl / eta.v should show little residual trend
#' with weight (a flat cloud) — the covariate has captured the size effect.
cov_plot_eta_wt <- function(fit, dat = build_nonmem_data()) {
  wt <- dat |>
    dplyr::distinct(ID, WT) |>
    dplyr::mutate(ID = as.integer(as.character(ID)))
  e  <- iiv_etas(fit)
  id_col <- intersect(c("ID", "id"), names(e))
  e <- dplyr::rename(e, ID = dplyr::all_of(id_col)) |>
    dplyr::mutate(ID = as.integer(as.character(ID)))
  d <- dplyr::left_join(e, wt, by = "ID") |>
    dplyr::select(ID, WT, dplyr::matches("^eta\\.(cl|v)$")) |>
    tidyr::pivot_longer(dplyr::starts_with("eta"),
                        names_to = "eta", values_to = "value")
  ggplot2::ggplot(d, ggplot2::aes(WT, value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "#E74C3C",
                         formula = y ~ x) +
    ggplot2::geom_point(colour = "#2C3E50", size = 2) +
    ggplot2::facet_wrap(~ eta, scales = "free_y") +
    ggplot2::labs(
      x = "Body weight (kg)", y = "Random effect (eta)",
      title = "Residual random effects vs weight (covariate model)",
      subtitle = "A flat red line = weight effect captured by the covariate"
    )
}

# -----------------------------------------------------------------------------
# 5. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- build_nonmem_data()

  iiv_fit <- load_iiv_fit()
  cov_fit <- fit_cov_model(dat)
  saveRDS(cov_fit, fs::path(PATH_MODELS, "cov_fit.rds"))

  cmp <- compare_models(
    model_metrics(iiv_fit, "IIV (no covariate)"),
    model_metrics(cov_fit, "Weight on CL & V")
  )
  readr::write_csv(cmp, fs::path(PATH_TABLES, "model_comparison_iiv_cov.csv"))

  save_fig(cov_plot_eta_wt(cov_fit, dat), "05_cov_eta_vs_weight")

  print(cmp)
  message("\n✓ Phase 5 covariate model fitted. dOFV vs IIV = ",
          cmp$dOFV[2], " | fit saved to outputs/models/cov_fit.rds")
}
