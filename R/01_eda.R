# =============================================================================
# 01_eda.R
# Project : Population PK Portfolio — Theophylline
# Phase 1 : Exploratory Data Analysis
# Purpose : Prepare the Theoph dataset, run data-quality checks, build a
#           per-subject summary, and produce the core EDA figures.
# Usage   : source(here::here("R", "01_eda.R"))   # defines functions
#           Then call prepare_theoph(), run_dq_checks(), subject_summary(),
#           and the eda_plot_*() builders. Running this file directly (via
#           Rscript / Source) executes the demo block at the bottom.
# =============================================================================

source(here::here("R", "00_setup.R"))

# -----------------------------------------------------------------------------
# 1. prepare_theoph(): tidy the built-in Theoph dataset
# -----------------------------------------------------------------------------
#' Load and tidy the Theoph dataset.
#'
#' The base `datasets::Theoph` object stores `Subject` as an ordered factor
#' whose level order is NOT 1..12. We coerce to an integer ID via the factor
#' label (the actual subject number), and derive the absolute dose in mg.
#'
#' Columns returned:
#'   id     integer subject id (1..12)
#'   wt     body weight (kg)
#'   dose   administered dose (mg/kg)
#'   amt    absolute dose (mg) = dose * wt
#'   time   time since dose (h)
#'   conc   theophylline concentration (mg/L)
#'
#' @return A tibble ordered by id, time.
prepare_theoph <- function() {
  raw <- datasets::Theoph

  df <- tibble::as_tibble(raw) |>
    janitor::clean_names() |>
    dplyr::transmute(
      id   = as.integer(as.character(subject)),
      wt   = as.numeric(wt),
      dose = as.numeric(dose),
      amt  = as.numeric(dose) * as.numeric(wt),   # mg/kg * kg = mg
      time = as.numeric(time),
      conc = as.numeric(conc)
    ) |>
    dplyr::arrange(id, time)

  df
}

# -----------------------------------------------------------------------------
# 2. run_dq_checks(): data-quality / integrity checks
# -----------------------------------------------------------------------------
#' Run a battery of data-quality checks on a prepared Theoph tibble.
#'
#' @param df Output of prepare_theoph().
#' @return A tibble with one row per check: check, status (PASS/FLAG), detail.
run_dq_checks <- function(df) {
  add <- function(check, ok, detail) {
    tibble::tibble(check = check,
                   status = ifelse(ok, "PASS", "FLAG"),
                   detail = detail)
  }

  n_subj   <- dplyr::n_distinct(df$id)
  per_subj <- df |> dplyr::count(id, name = "n_obs")
  dup      <- df |> dplyr::count(id, time) |> dplyr::filter(n > 1)
  miss     <- colSums(is.na(df))
  neg_conc <- df |> dplyr::filter(conc < 0)
  unsorted <- df |>
    dplyr::group_by(id) |>
    dplyr::summarise(ok = !is.unsorted(time), .groups = "drop") |>
    dplyr::filter(!ok)
  baseline <- df |>
    dplyr::group_by(id) |>
    dplyr::summarise(t0 = min(time), c0 = conc[which.min(time)], .groups = "drop")
  wt_rng   <- range(df$wt)
  dose_rng <- range(df$dose)

  dplyr::bind_rows(
    add("Subject count",
        n_subj == 12,
        glue::glue("{n_subj} unique subjects (expected 12)")),
    add("Samples per subject",
        all(per_subj$n_obs == 11),
        glue::glue("range {min(per_subj$n_obs)}–{max(per_subj$n_obs)} ",
                   "(Theoph has 11 samples/subject)")),
    add("Duplicate id×time rows",
        nrow(dup) == 0,
        glue::glue("{nrow(dup)} duplicate time points")),
    add("Missing values",
        sum(miss) == 0,
        glue::glue("{sum(miss)} NA cells across all columns")),
    add("Negative concentrations",
        nrow(neg_conc) == 0,
        glue::glue("{nrow(neg_conc)} negative conc values")),
    add("Time monotonic within subject",
        nrow(unsorted) == 0,
        glue::glue("{nrow(unsorted)} subjects with non-increasing time")),
    add("Baseline (t=0) concentrations",
        all(baseline$c0 >= 0),
        glue::glue("first-sample conc range ",
                   "{round(min(baseline$c0),2)}–{round(max(baseline$c0),2)} mg/L")),
    add("Weight range plausible",
        wt_rng[1] > 30 & wt_rng[2] < 120,
        glue::glue("{wt_rng[1]}–{wt_rng[2]} kg")),
    add("Dose range plausible",
        dose_rng[1] > 1 & dose_rng[2] < 12,
        glue::glue("{dose_rng[1]}–{dose_rng[2]} mg/kg"))
  )
}

# -----------------------------------------------------------------------------
# 3. subject_summary(): per-subject descriptive summary
# -----------------------------------------------------------------------------
#' Build a per-subject summary (demographics + observed PK landmarks).
#'
#' @param df Output of prepare_theoph().
#' @return A tibble, one row per subject.
subject_summary <- function(df) {
  df |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      wt        = dplyr::first(wt),
      dose_mgkg = dplyr::first(dose),
      amt_mg    = dplyr::first(amt),
      n_obs     = dplyr::n(),
      tmax_obs  = time[which.max(conc)],
      cmax_obs  = max(conc),
      tlast     = max(time),
      clast     = conc[which.max(time)],
      .groups   = "drop"
    ) |>
    dplyr::arrange(id)
}

# -----------------------------------------------------------------------------
# 4. EDA plot builders
# -----------------------------------------------------------------------------

#' Individual concentration–time profiles, faceted by subject (linear scale).
eda_plot_individual <- function(df, semilog = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(time, conc, group = id)) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 0.5) +
    ggplot2::geom_point(colour = "#E74C3C", size = 1.1) +
    ggplot2::facet_wrap(~ id, ncol = 4, labeller = ggplot2::label_both) +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = if (semilog) "Individual profiles (semi-log)" else "Individual profiles (linear)",
      subtitle = "Theophylline — single oral dose, 12 subjects"
    )
  if (semilog) {
    p <- p + ggplot2::scale_y_log10()
  }
  p
}

#' Spaghetti plot: all subjects overlaid, coloured by subject.
eda_plot_spaghetti <- function(df, semilog = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(time, conc, colour = factor(id), group = id)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    ggplot2::geom_point(size = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo", name = "Subject") +
    ggplot2::labs(
      x = "Time (h)", y = "Concentration (mg/L)",
      title = if (semilog) "All subjects overlaid (semi-log)" else "All subjects overlaid (linear)",
      subtitle = "Spaghetti plot of individual concentration–time data"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
  if (semilog) p <- p + ggplot2::scale_y_log10()
  p
}

#' Mean (± SD) concentration–time profile across subjects.
eda_plot_mean <- function(df) {
  mean_df <- df |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      mean_conc = mean(conc),
      sd_conc   = stats::sd(conc),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  ggplot2::ggplot(mean_df, ggplot2::aes(time, mean_conc)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(mean_conc - sd_conc, 0),
                                      ymax = mean_conc + sd_conc),
                         fill = "#3498DB", alpha = 0.25) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 0.9) +
    ggplot2::geom_point(colour = "#2C3E50", size = 1.6) +
    ggplot2::labs(
      x = "Time (h)", y = "Mean concentration (mg/L)",
      title = "Mean concentration–time profile (± SD)",
      subtitle = "Pooled across all 12 subjects at nominal sampling times"
    )
}

#' Body-weight distribution across subjects.
eda_plot_weight <- function(df) {
  wt_df <- df |> dplyr::distinct(id, wt)
  ggplot2::ggplot(wt_df, ggplot2::aes(wt)) +
    ggplot2::geom_histogram(binwidth = 5, fill = "#3498DB",
                            colour = "white", boundary = 0) +
    ggplot2::geom_rug(colour = "#E74C3C") +
    ggplot2::labs(
      x = "Body weight (kg)", y = "Number of subjects",
      title = "Body-weight distribution",
      subtitle = glue::glue("n = {nrow(wt_df)} subjects; ",
                            "range {min(wt_df$wt)}–{max(wt_df$wt)} kg")
    )
}

#' Observed Cmax vs body weight, with subject labels.
eda_plot_cmax_weight <- function(ss) {
  ggplot2::ggplot(ss, ggplot2::aes(wt, cmax_obs)) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#3498DB",
                         fill = "#3498DB", alpha = 0.15, formula = y ~ x) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(id)), size = 2.6) +
    ggrepel::geom_text_repel(ggplot2::aes(label = id), size = 3, seed = 1) +
    viridis::scale_colour_viridis(discrete = TRUE, option = "turbo", guide = "none") +
    ggplot2::labs(
      x = "Body weight (kg)", y = "Observed Cmax (mg/L)",
      title = "Observed Cmax vs body weight",
      subtitle = "Heavier subjects received more total drug (dose is mg/kg)"
    )
}

# -----------------------------------------------------------------------------
# 5. Demo / pipeline block (runs when the file is sourced directly)
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  theoph    <- prepare_theoph()
  dq        <- run_dq_checks(theoph)
  ss        <- subject_summary(theoph)

  # Persist the prepared dataset for downstream phases
  readr::write_csv(theoph, fs::path(PATH_DATA_PROC, "theoph_clean.csv"))
  readr::write_csv(ss,     fs::path(PATH_TABLES,    "subject_summary.csv"))

  # Build and save figures
  save_fig(eda_plot_individual(theoph, semilog = FALSE), "01_individual_linear")
  save_fig(eda_plot_individual(theoph, semilog = TRUE),  "01_individual_semilog")
  save_fig(eda_plot_mean(theoph),                        "01_mean_profile")
  save_fig(eda_plot_spaghetti(theoph),                   "01_spaghetti")
  save_fig(eda_plot_weight(theoph),                      "01_weight_dist")
  save_fig(eda_plot_cmax_weight(ss),                     "01_cmax_vs_weight")

  print(dq)
  message("✓ Phase 1 EDA pipeline complete. Figures in outputs/figures/, ",
          "tables in outputs/tables/.")
}
