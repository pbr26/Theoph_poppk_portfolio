# =============================================================================
# 06_diagnostics.R
# Project : Population PK Portfolio — Theophylline
# Phase 6 : Model diagnostics
# Purpose : Goodness-of-fit plots (DV~PRED, DV~IPRED, CWRES~TIME, CWRES~PRED)
#           and a Visual Predictive Check (VPC) for the selected model.
# Usage   : source(here::here("R", "06_diagnostics.R"))
# =============================================================================

source(here::here("R", "00_setup.R"))
stopifnot(requireNamespace("nlmixr2", quietly = TRUE))
library(nlmixr2)

# -----------------------------------------------------------------------------
# 1. load_final_fit(): the most complete model available
# -----------------------------------------------------------------------------
#' Prefer the covariate model, then IIV, then base.
load_final_fit <- function() {
  for (f in c("cov_fit.rds", "iiv_fit.rds", "base_fit.rds")) {
    p <- fs::path(PATH_MODELS, f)
    if (fs::file_exists(p)) {
      message("Diagnostics on: ", f)
      return(readRDS(p))
    }
  }
  stop("No fitted model found in outputs/models/. Run Phases 3-5 first.")
}

# -----------------------------------------------------------------------------
# 2. gof_frame(): observation-level predictions / residuals
# -----------------------------------------------------------------------------
gof_frame <- function(fit) {
  d <- tibble::as_tibble(as.data.frame(fit))
  needed <- c("DV", "PRED", "IPRED", "CWRES", "TIME")
  miss <- setdiff(needed, names(d))
  if (length(miss)) {
    warning("Missing columns in fit data: ", paste(miss, collapse = ", "))
  }
  d
}

# -----------------------------------------------------------------------------
# 3. GOF plot builders
# -----------------------------------------------------------------------------
.gof_obs_pred <- function(d, xvar, xlab) {
  rng <- range(c(d$DV, d[[xvar]]), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(.data[[xvar]], DV)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_point(alpha = 0.6, colour = "#2C3E50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#E74C3C",
                         linewidth = 0.7, formula = y ~ x) +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(x = xlab, y = "Observed (mg/L)")
}

.gof_cwres <- function(d, xvar, xlab) {
  ggplot2::ggplot(d, ggplot2::aes(.data[[xvar]], CWRES)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                        colour = "grey70") +
    ggplot2::geom_point(alpha = 0.6, colour = "#2C3E50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#E74C3C",
                         linewidth = 0.7, formula = y ~ x) +
    ggplot2::labs(x = xlab, y = "CWRES")
}

#' Four-panel GOF figure assembled with patchwork.
gof_panel <- function(fit) {
  d <- gof_frame(fit)
  p1 <- .gof_obs_pred(d, "PRED",  "Population prediction (mg/L)")
  p2 <- .gof_obs_pred(d, "IPRED", "Individual prediction (mg/L)")
  p3 <- .gof_cwres(d, "TIME", "Time (h)")
  p4 <- .gof_cwres(d, "PRED", "Population prediction (mg/L)")
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2) +
    patchwork::plot_annotation(
      title = "Goodness-of-fit diagnostics",
      subtitle = "Top: observed vs predictions (hug the dashed line). Bottom: CWRES (flat band within ±2)."
    )
}

# -----------------------------------------------------------------------------
# 4. Visual Predictive Check (guarded — needs the 'vpc' package)
# -----------------------------------------------------------------------------
make_vpc <- function(fit, n = 500) {
  if (!requireNamespace("vpc", quietly = TRUE)) {
    message("• 'vpc' package not installed — skipping VPC. ",
            "Install with: renv::install('vpc')")
    return(NULL)
  }
  # vpcPlot lives in nlmixr2plot (installed alongside nlmixr2)
  vpc_fun <- tryCatch(
    get("vpcPlot", envir = asNamespace("nlmixr2plot")),
    error = function(e) NULL
  )
  if (is.null(vpc_fun)) {
    message("• vpcPlot() not found in nlmixr2plot — skipping VPC.")
    return(NULL)
  }
  tryCatch(
    vpc_fun(fit, n = n, show = list(obs_dv = TRUE)) +
      ggplot2::labs(title = "Visual Predictive Check",
                    x = "Time (h)", y = "Concentration (mg/L)"),
    error = function(e) {
      message("• VPC failed: ", conditionMessage(e))
      NULL
    }
  )
}

# -----------------------------------------------------------------------------
# 5. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  fit <- load_final_fit()

  gof <- gof_panel(fit)
  save_fig(gof, "06_gof_panel", w = 9, h = 8)

  vpc <- make_vpc(fit, n = 500)
  if (!is.null(vpc)) save_fig(vpc, "06_vpc", w = 8, h = 6)

  message("\n✓ Phase 6 diagnostics complete. ",
          "GOF panel saved; VPC ", if (is.null(vpc)) "skipped." else "saved.")
}
