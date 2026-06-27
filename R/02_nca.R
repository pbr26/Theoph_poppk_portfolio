# =============================================================================
# 02_nca.R
# Project : Population PK Portfolio — Theophylline
# Phase 2 : Noncompartmental Analysis (NCA) with PKNCA
# Purpose : Compute Cmax, Tmax, AUClast, AUCinf, t1/2 (and CL/Vz) per subject,
#           tidy the results, summarise, and build NCA figures.
# Usage   : source(here::here("R", "02_nca.R"))   # defines functions
#           Sourcing directly runs the pipeline (demo block at the bottom).
# =============================================================================

source(here::here("R", "00_setup.R"))
stopifnot(requireNamespace("PKNCA", quietly = TRUE))
library(PKNCA)

# -----------------------------------------------------------------------------
# 1. load_pk_data(): read the cleaned dataset from Phase 1 (or rebuild it)
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
# 2. run_nca(): full PKNCA pipeline -> tidy per-subject parameter table
# -----------------------------------------------------------------------------
#' Run noncompartmental analysis on the Theoph data.
#'
#' Builds PKNCAconc (concentration) and PKNCAdose (dose) objects, defines a
#' single 0–Inf interval requesting the key NCA parameters, and runs pk.nca().
#'
#' @param df Prepared data (id, time, conc, amt). Defaults to load_pk_data().
#' @return list(result = <PKNCAresults>, wide = <tibble one row per subject>,
#'              long = <tidy long results>)
run_nca <- function(df = load_pk_data()) {

  # Concentration object: conc vs time, grouped by subject
  conc_obj <- PKNCA::PKNCAconc(as.data.frame(df), conc ~ time | id)

  # Dose object: one administration per subject at time 0 (absolute mg)
  dose_df <- df |>
    dplyr::distinct(id, amt) |>
    dplyr::mutate(time = 0) |>
    as.data.frame()
  dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | id)

  # Interval 0 -> Inf requesting the parameters of interest
  intervals <- data.frame(
    start      = 0,
    end        = Inf,
    cmax       = TRUE,   # peak concentration
    tmax       = TRUE,   # time of peak
    auclast    = TRUE,   # AUC to last measurable conc
    aucinf.obs = TRUE,   # AUC extrapolated to infinity (obs-based)
    half.life  = TRUE,   # terminal half-life (from lambda.z)
    lambda.z   = TRUE,   # terminal rate constant
    cl.obs     = TRUE,   # apparent clearance = dose / AUCinf
    vz.obs     = TRUE,   # apparent volume = CL / lambda.z
    aucpext.obs = TRUE   # % AUC extrapolated (quality metric)
  )

  data_obj <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
  res      <- PKNCA::pk.nca(data_obj)

  long <- as.data.frame(res) |>
    tibble::as_tibble() |>
    dplyr::select(id, PPTESTCD, PPORRES)

  wide <- long |>
    tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
    dplyr::arrange(id)

  list(result = res, wide = wide, long = long)
}

# -----------------------------------------------------------------------------
# 3. nca_param_summary(): geometric mean + CV% across subjects
# -----------------------------------------------------------------------------
#' Geometric mean and geometric CV% for a numeric vector.
geo_mean <- function(x) exp(mean(log(x[x > 0 & is.finite(x)])))
geo_cv   <- function(x) sqrt(exp(stats::var(log(x[x > 0 & is.finite(x)]))) - 1) * 100

#' Summarise key NCA parameters across subjects (geomean, CV%, range).
nca_param_summary <- function(wide) {
  params <- intersect(
    c("cmax", "tmax", "auclast", "aucinf.obs", "half.life",
      "cl.obs", "vz.obs", "aucpext.obs"),
    names(wide)
  )
  purrr::map_dfr(params, function(p) {
    x <- wide[[p]]
    tibble::tibble(
      Parameter   = p,
      N           = sum(is.finite(x)),
      `Geo. mean` = geo_mean(x),
      `Geo. CV%`  = geo_cv(x),
      Min         = min(x, na.rm = TRUE),
      Median      = stats::median(x, na.rm = TRUE),
      Max         = max(x, na.rm = TRUE)
    )
  })
}

# -----------------------------------------------------------------------------
# 4. NCA plot builders
# -----------------------------------------------------------------------------

#' Per-subject terminal half-life, ordered.
nca_plot_halflife <- function(wide) {
  d <- wide |> dplyr::mutate(id = forcats::fct_reorder(factor(id), half.life))
  ggplot2::ggplot(d, ggplot2::aes(half.life, id)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = half.life, yend = id),
                          colour = "grey75") +
    ggplot2::geom_point(ggplot2::aes(colour = half.life), size = 3) +
    viridis::scale_colour_viridis(option = "turbo", guide = "none") +
    ggplot2::geom_vline(xintercept = geo_mean(wide$half.life),
                        linetype = "dashed", colour = "#E74C3C") +
    ggplot2::labs(
      x = "Terminal half-life (h)", y = "Subject",
      title = "Terminal half-life by subject",
      subtitle = "Dashed line = geometric mean; theophylline t½ is typically ~6–9 h"
    )
}

#' AUClast vs AUCinf — closeness indicates little extrapolation.
nca_plot_auc <- function(wide) {
  rng <- range(c(wide$auclast, wide$aucinf.obs), na.rm = TRUE)
  ggplot2::ggplot(wide, ggplot2::aes(auclast, aucinf.obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = factor(id)), size = 2.6) +
    ggrepel::geom_text_repel(ggplot2::aes(label = id), size = 3, seed = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo", guide = "none") +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(
      x = "AUC(last) (mg·h/L)", y = "AUC(inf, obs) (mg·h/L)",
      title = "AUClast vs AUCinf",
      subtitle = "Points near the line of identity = minimal extrapolation"
    )
}

#' Apparent clearance vs body weight (motivates the weight covariate).
nca_plot_cl_weight <- function(wide, df = load_pk_data()) {
  wt <- df |> dplyr::distinct(id, wt)
  d  <- dplyr::left_join(wide, wt, by = "id")
  ggplot2::ggplot(d, ggplot2::aes(wt, cl.obs)) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#3498DB",
                         fill = "#3498DB", alpha = 0.15, formula = y ~ x) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(id)), size = 2.6) +
    ggrepel::geom_text_repel(ggplot2::aes(label = id), size = 3, seed = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo", guide = "none") +
    ggplot2::labs(
      x = "Body weight (kg)", y = "Apparent clearance CL/F (L/h)",
      title = "NCA clearance vs body weight",
      subtitle = "Exploratory check of weight as a clearance covariate (Phase 5)"
    )
}

#' Boxplots of the key NCA parameters (distribution overview).
nca_plot_distributions <- function(wide) {
  keep <- intersect(c("cmax", "tmax", "auclast", "aucinf.obs",
                      "half.life", "cl.obs", "vz.obs"), names(wide))
  long <- wide |>
    dplyr::select(id, dplyr::all_of(keep)) |>
    tidyr::pivot_longer(-id, names_to = "param", values_to = "value")
  ggplot2::ggplot(long, ggplot2::aes(param, value)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#3498DB", alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.7, colour = "#2C3E50") +
    ggplot2::facet_wrap(~ param, scales = "free", ncol = 4) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Distribution of NCA parameters across subjects"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

# -----------------------------------------------------------------------------
# 5. Demo / pipeline block
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  df  <- load_pk_data()
  nca <- run_nca(df)
  smry <- nca_param_summary(nca$wide)

  readr::write_csv(nca$wide, fs::path(PATH_TABLES, "nca_parameters.csv"))
  readr::write_csv(smry,     fs::path(PATH_TABLES, "nca_summary.csv"))

  save_fig(nca_plot_halflife(nca$wide),        "02_nca_halflife")
  save_fig(nca_plot_auc(nca$wide),             "02_nca_auc")
  save_fig(nca_plot_cl_weight(nca$wide, df),   "02_nca_cl_weight")
  save_fig(nca_plot_distributions(nca$wide),   "02_nca_distributions")

  print(smry)
  message("✓ Phase 2 NCA pipeline complete. ",
          "Parameters in outputs/tables/nca_parameters.csv.")
}
