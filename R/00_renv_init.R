# =============================================================================
# 00_renv_init.R
# Project : Population PK Portfolio — Theophylline
# Purpose : Bootstrap renv and install project packages in TWO stages.
#
#   STAGE A  — pure-R / lightweight stack (Phases 0-2). Safe to snapshot.
#   STAGE B  — compiled modelling stack rxode2 + nlmixr2 (Phases 3-7).
#              Requires a working gfortran toolchain. Do NOT snapshot until
#              a trivial rxode2 model compiles successfully (see Phase 3).
#
# Run ONCE for first-time setup. Afterwards: renv::restore().
# =============================================================================

# ── Step 1: renv itself ──────────────────────────────────────────────────────
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# ── Step 2: Initialise renv in this project (creates renv/activate.R) ─────────
# Skip if already initialised.
if (!file.exists(here::here("renv", "activate.R"))) {
  renv::init(bare = TRUE)
}

# ── Step 3a: STAGE A — install + snapshot the lightweight stack ───────────────
packages_stage_a <- c(
  # Core data wrangling
  "tidyverse",
  # Utilities
  "here", "fs", "glue", "janitor",
  # Visualization
  "ggplot2", "viridis", "ggrepel", "patchwork", "scales",
  # Tables & reporting
  "knitr", "kableExtra", "rmarkdown", "quarto",
  # Phase 2 NCA (pure R, no Fortran)
  "PKNCA"
)
renv::install(packages_stage_a)

# Snapshot now — this is the reproducible baseline committed to Git.
renv::snapshot(prompt = FALSE)
message("✓ STAGE A installed and snapshotted (renv.lock updated).")

# ── Step 3b: STAGE B — compiled modelling stack (run AFTER gfortran check) ────
# rxode2 / nlmixr2 compile C++/Fortran. On macOS Apple Silicon, rxode2 expects
# gfortran at /opt/gfortran/bin/gfortran. If you used Homebrew gfortran, create
# the symlink (run once in Terminal):
#
#   sudo mkdir -p /opt/gfortran/bin
#   sudo ln -sf /opt/homebrew/bin/gfortran /opt/gfortran/bin/gfortran
#
# Then uncomment and run the lines below:
#
# renv::install(c("rxode2", "nlmixr2"))
#
# Verify a trivial compile works (Phase 3 will reuse this check):
#   library(rxode2)
#   m <- rxode2({ d/dt(x) <- -k * x }); print(m)
#
# ONLY after that prints without error:
#   renv::snapshot(prompt = FALSE)
#
message("⚠ STAGE B (rxode2/nlmixr2) is intentionally NOT installed here.\n",
        "  Confirm the gfortran toolchain first, then run it in Phase 3.")
