# =============================================================================
# 09_bootstrap.R
# Project : Population PK Portfolio — Theophylline
# Phase 10 : Bootstrap confidence intervals
# Purpose : Case (subject) resampling bootstrap of the selected IIV model to
#           obtain empirical 95% CIs, compared with the asymptotic CIs.
# Usage   : source(here::here("R", "09_bootstrap.R"))
# Note    : Refits the model nboot times — expect a few minutes.
# =============================================================================

source(here::here("R", "04_iiv_model.R"))   # fit_iiv_model(), build_nonmem_data(), setup

# Fixed effects estimated on the log scale (back-transform with exp()).
LOG_PARS <- c("tka", "tcl", "tv")

# -----------------------------------------------------------------------------
# 1. load_iiv_fit(): original fit (point estimates + asymptotic CIs)
# -----------------------------------------------------------------------------
load_iiv_fit <- function() {
  p <- fs::path(PATH_MODELS, "iiv_fit.rds")
  if (fs::file_exists(p)) readRDS(p) else fit_iiv_model()
}

# -----------------------------------------------------------------------------
# 2. .resample(): one subject-resampled dataset (new sequential IDs)
# -----------------------------------------------------------------------------
.resample <- function(dat) {
  ids <- unique(dat$ID)
  samp <- sample(ids, length(ids), replace = TRUE)
  purrr::imap_dfr(samp, function(id, i) {
    d <- dat[dat$ID == id, , drop = FALSE]
    d$ID <- i          # unique IDs so resampled subjects stay distinct
    d
  })
}

# -----------------------------------------------------------------------------
# 3. run_bootstrap(): refit on nboot resamples, collect fixed effects
# -----------------------------------------------------------------------------
#' @return list(theta = matrix [nsuccess x npar], n_success, n_attempt)
run_bootstrap <- function(dat = build_nonmem_data(), nboot = 50, seed = 20260627) {
  set.seed(seed)
  rows <- vector("list", nboot)
  ok <- 0L
  for (b in seq_len(nboot)) {
    bd  <- .resample(dat)
    fitb <- tryCatch(
      suppressWarnings(suppressMessages(fit_iiv_model(bd))),
      error = function(e) NULL
    )
    th <- tryCatch(fitb$theta, error = function(e) NULL)
    if (!is.null(th)) { ok <- ok + 1L; rows[[b]] <- th }
  }
  mat <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  list(theta = mat, n_success = ok, n_attempt = nboot)
}

# -----------------------------------------------------------------------------
# 4. boot_summary(): asymptotic vs bootstrap 95% CI (back-transformed)
# -----------------------------------------------------------------------------
boot_summary <- function(boot, fit = load_iiv_fit()) {
  mat <- boot$theta
  th0 <- fit$theta
  params <- intersect(colnames(mat), names(th0))

  purrr::map_dfr(params, function(p) {
    x   <- mat[, p]
    bt  <- p %in% LOG_PARS
    tf  <- if (bt) exp else identity
    est <- tf(as.numeric(th0[[p]]))
    qs  <- stats::quantile(tf(x), c(0.025, 0.5, 0.975), names = FALSE)
    tibble::tibble(
      Parameter      = if (bt) sub("^t", "", p) else p,
      `Estimate`     = est,
      `Boot median`  = qs[2],
      `Boot 2.5%`    = qs[1],
      `Boot 97.5%`   = qs[3]
    )
  })
}

# -----------------------------------------------------------------------------
# 5. boot_plot(): bootstrap distribution of the structural parameters
# -----------------------------------------------------------------------------
boot_plot <- function(boot) {
  mat <- boot$theta
  keep <- intersect(LOG_PARS, colnames(mat))
  long <- tibble::as_tibble(exp(mat[, keep, drop = FALSE])) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "param", values_to = "value") |>
    dplyr::mutate(param = sub("^t", "", param))
  ci <- long |>
    dplyr::group_by(param) |>
    dplyr::summarise(lo = stats::quantile(value, 0.025),
                     hi = stats::quantile(value, 0.975), .groups = "drop")
  ggplot2::ggplot(long, ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 25, fill = "#3498DB", colour = "white", alpha = 0.7) +
    ggplot2::geom_vline(data = ci, ggplot2::aes(xintercept = lo),
                        linetype = "dashed", colour = "#E74C3C") +
    ggplot2::geom_vline(data = ci, ggplot2::aes(xintercept = hi),
                        linetype = "dashed", colour = "#E74C3C") +
    ggplot2::facet_wrap(~ param, scales = "free") +
    ggplot2::labs(
      x = "Bootstrap estimate", y = "Count",
      title = "Bootstrap distributions of structural parameters",
      subtitle = "Dashed red lines = empirical 95% percentile interval"
    )
}

# -----------------------------------------------------------------------------
# 6. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  dat <- build_nonmem_data()
  fit <- load_iiv_fit()

  NBOOT <- 50
  message("Running ", NBOOT, "-replicate bootstrap (this takes a few minutes)…")
  boot <- run_bootstrap(dat, nboot = NBOOT)

  smry <- boot_summary(boot, fit)
  readr::write_csv(smry, fs::path(PATH_TABLES, "bootstrap_summary.csv"))
  save_fig(boot_plot(boot), "09_bootstrap_distributions")

  print(smry)
  message("\n✓ Phase 10 bootstrap complete. ",
          boot$n_success, "/", boot$n_attempt, " replicates converged.")
}
