# =============================================================================
# 10_forest_plot.R
# Project : Population PK Portfolio — Theophylline
# Phase 11 : Covariate forest plot
# Purpose : Summarise the weight effect on CL and V as parameter ratios
#           (relative to a 70 kg typical subject) with 95% CIs, for a light
#           and a heavy subject — the standard clinical-relevance figure.
# Usage   : source(here::here("R", "10_forest_plot.R"))
# =============================================================================

source(here::here("R", "05_covariate_model.R"))   # fit_cov_model(), build_nonmem_data(), setup

REF_WT  <- 70                 # reference weight (kg)
BAND_LO <- 0.8                # clinical "no-effect" band (lower)
BAND_HI <- 1.25               # clinical "no-effect" band (upper)

# -----------------------------------------------------------------------------
# 1. load_cov_fit(): covariate model (refit if the saved fit is missing)
# -----------------------------------------------------------------------------
load_cov_fit <- function() {
  p <- fs::path(PATH_MODELS, "cov_fit.rds")
  if (fs::file_exists(p)) readRDS(p) else fit_cov_model()
}

# -----------------------------------------------------------------------------
# 2. forest_data(): parameter ratios vs typical, with 95% CIs
# -----------------------------------------------------------------------------
#' Ratio of a parameter at weight w to its typical (70 kg) value is
#' (w/REF_WT)^theta. The CI is propagated from the exponent's 95% CI
#' (theta ± 1.96·SE), taking the min/max across the exponent bounds.
forest_data <- function(fit = load_cov_fit(),
                        weights = c(light = 55, heavy = 85)) {
  th <- fit$theta
  se <- sqrt(diag(as.matrix(fit$cov)))   # SE of fixed effects

  spec <- tibble::tribble(
    ~param, ~exp_name, ~label,
    "CL",   "dWTCL",   "Clearance (CL)",
    "V",    "dWTV",    "Volume (V)"
  )

  purrr::pmap_dfr(spec, function(param, exp_name, label) {
    e  <- as.numeric(th[[exp_name]])
    s  <- as.numeric(se[[exp_name]])
    lo <- e - 1.96 * s
    hi <- e + 1.96 * s
    purrr::imap_dfr(weights, function(w, wname) {
      rs <- (w / REF_WT)^c(e, lo, hi)
      tibble::tibble(
        Parameter = label,
        Scenario  = paste0(wname, " (", w, " kg)"),
        ratio     = rs[1],
        lo        = min(rs[2], rs[3]),
        hi        = max(rs[2], rs[3])
      )
    })
  })
}

# -----------------------------------------------------------------------------
# 3. forest_plot(): the figure
# -----------------------------------------------------------------------------
forest_plot <- function(df) {
  df <- df |>
    dplyr::mutate(row = paste(Parameter, "—", Scenario),
                  row = forcats::fct_rev(factor(row)))
  ggplot2::ggplot(df, ggplot2::aes(ratio, row)) +
    ggplot2::annotate("rect", xmin = BAND_LO, xmax = BAND_HI,
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                            height = 0.18, colour = "#2C3E50") +
    ggplot2::geom_point(ggplot2::aes(colour = Parameter), size = 3) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo", end = 0.7,
                                  guide = "none") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Parameter ratio vs typical 70 kg subject (log scale)",
      y = NULL,
      title = "Covariate forest plot — body weight effect on CL and V",
      subtitle = "Shaded band = 0.8–1.25 (no clinically relevant effect); dashed line = no change"
    )
}

# -----------------------------------------------------------------------------
# 4. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  fit <- load_cov_fit()
  fd  <- forest_data(fit)

  readr::write_csv(fd, fs::path(PATH_TABLES, "covariate_forest.csv"))
  save_fig(forest_plot(fd), "10_covariate_forest", w = 9, h = 5)

  print(fd)
  message("\n✓ Phase 11 forest plot complete. ",
          "Figure: outputs/figures/10_covariate_forest.png")
}
