# =============================================================================
# 00_setup.R
# Project : Population PK Portfolio — Theophylline
# Purpose : Load packages, set global options, define the ggplot2 theme,
#           save_fig() helper, and here::here()-based project paths.
# Author  : Pramod BR
# Notes   : Source this at the top of every script / report:
#               source(here::here("R", "00_setup.R"))
# =============================================================================

# ── 0. Guarded package loader ────────────────────────────────────────────────
# Loads a package if installed; otherwise warns instead of hard-erroring.
# This lets early phases run before the compiled modelling stack
# (rxode2 / nlmixr2) is confirmed working.
.need <- function(pkg, critical = TRUE) {
  ok <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
  if (ok) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  } else if (critical) {
    stop(sprintf("Required package '%s' is not installed. Run R/00_renv_init.R.", pkg),
         call. = FALSE)
  } else {
    message(sprintf("• Optional package '%s' not available yet — skipping.", pkg))
  }
  invisible(ok)
}

# ── 1. Package loading ───────────────────────────────────────────────────────
# Core data wrangling + reporting (required for all phases)
.need("tidyverse")     # dplyr, ggplot2, tidyr, readr, purrr, tibble, stringr
.need("here")          # project-relative paths
.need("fs")            # filesystem operations
.need("janitor")       # clean_names(), tabyl()
.need("glue")          # string interpolation
.need("scales")        # axis / label formatting

# Visualization
.need("ggplot2")
.need("viridis")       # colour-blind-safe palettes
.need("ggrepel")       # non-overlapping text labels
.need("patchwork")     # compose multiple ggplots

# Tables
.need("knitr")
.need("kableExtra")

# Pharmacometrics (loaded but NON-critical until the toolchain is verified)
.need("PKNCA",   critical = FALSE)   # noncompartmental analysis (Phase 2)
.need("rxode2",  critical = FALSE)   # ODE engine (Phases 3-7)
.need("nlmixr2", critical = FALSE)   # population PK estimation (Phases 3-7)

# ── 2. Global knitr chunk options ────────────────────────────────────────────
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(
    echo       = TRUE,
    warning    = FALSE,
    message    = FALSE,
    fig.width  = 8,
    fig.height = 6,
    dpi        = 300,
    out.width  = "100%"
  )
}

# ── 3. Global ggplot2 theme: theme_pk() ──────────────────────────────────────
#' A clean publication theme for PK figures.
#' @param base_size Base font size in points.
theme_pk <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.4),
      strip.background = ggplot2::element_rect(fill = "#2C3E50", colour = NA),
      strip.text       = ggplot2::element_text(colour = "white", face = "bold",
                                               size = ggplot2::rel(0.9)),
      axis.title       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1)),
      plot.subtitle    = ggplot2::element_text(colour = "grey40"),
      plot.caption     = ggplot2::element_text(colour = "grey50", size = ggplot2::rel(0.8)),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.8, "lines")
    )
}
ggplot2::theme_set(theme_pk())

# Consistent fill/colour scales for up to 12 Theoph subjects
SUBJECT_PAL <- viridisLite::viridis(12, option = "turbo")

# ── 4. Project paths (here::here) ────────────────────────────────────────────
PATH_DATA_RAW  <- here::here("data", "raw")
PATH_DATA_PROC <- here::here("data", "processed")
PATH_FIGS      <- here::here("outputs", "figures")
PATH_TABLES    <- here::here("outputs", "tables")
PATH_MODELS    <- here::here("outputs", "models")

fs::dir_create(c(PATH_DATA_RAW, PATH_DATA_PROC, PATH_FIGS, PATH_TABLES, PATH_MODELS))

# ── 5. save_fig(): write a ggplot to outputs/figures in PNG (+ optional PDF) ──
#' @param plot A ggplot object.
#' @param name Base filename without extension.
#' @param w,h  Width / height in inches.
#' @param pdf  Also write a vector PDF (via cairo_pdf).
save_fig <- function(plot, name, w = 8, h = 6, pdf = TRUE) {
  png_path <- fs::path(PATH_FIGS, paste0(name, ".png"))
  ggplot2::ggsave(png_path, plot = plot, width = w, height = h, dpi = 300, bg = "white")
  if (isTRUE(pdf)) {
    pdf_path <- fs::path(PATH_FIGS, paste0(name, ".pdf"))
    # Use the base pdf device (no cairo dependency — cairo fails to load on
    # this build even though capabilities() reports it as available).
    tryCatch(
      suppressWarnings(
        ggplot2::ggsave(pdf_path, plot = plot, width = w, height = h,
                        device = grDevices::pdf)
      ),
      error = function(e)
        message(sprintf("• PDF skipped for '%s' (%s). PNG was still written.",
                        name, conditionMessage(e)))
    )
  }
  invisible(plot)
}

# ── 6. Reproducibility ───────────────────────────────────────────────────────
set.seed(20260627)

message("✓ 00_setup.R loaded: theme_pk(), save_fig(), and project paths ready.")
