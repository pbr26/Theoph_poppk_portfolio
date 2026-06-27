# Population PK Portfolio — Theophylline

A complete, reproducible end-to-end **Population Pharmacokinetics (PopPK)** analysis of the
classic `Theoph` dataset, built in R with `nlmixr2` / `rxode2`, `PKNCA`, and Quarto.

> **Live site:** published from `docs/` via GitHub Pages once Phase 8 is complete.

## Why Theophylline?

Theophylline is a narrow-therapeutic-index bronchodilator with well-characterised
one-compartment oral PK. The built-in `Theoph` dataset has rich serial sampling from 12
subjects after a single oral dose, with body weight as a covariate of interest — a benchmark
dataset for PopPK methodology.

## Project structure

```
theoph_poppk_portfolio/
├── theoph_poppk_portfolio.Rproj   # RStudio project (Website build type)
├── _quarto.yml                    # Quarto website config (output-dir: docs)
├── renv.lock                      # Locked package environment
├── index.qmd                      # Portfolio landing page
├── R/                             # Analysis scripts (one per phase)
│   ├── 00_setup.R                 # packages, theme_pk(), save_fig(), paths
│   ├── 00_renv_init.R             # one-time environment bootstrap
│   └── 01_eda.R … 07_simulation.R
├── reports/                       # Quarto reports (.qmd) → rendered to docs/
├── data/{raw,processed}/          # data artifacts
├── outputs/{figures,tables,models}/
└── docs/                          # rendered website (GitHub Pages root)
```

## Analysis phases

| Phase | Topic | Tooling |
|------:|-------|---------|
| 0 | Project setup | renv, Quarto, here, fs |
| 1 | Exploratory data analysis | tidyverse, ggplot2, patchwork |
| 2 | Noncompartmental analysis | PKNCA |
| 3 | Base PopPK model (1-cmt oral) | nlmixr2, rxode2 |
| 4 | Inter-individual variability | nlmixr2 |
| 5 | Covariate model (weight) | nlmixr2 |
| 6 | Diagnostics (GOF, VPC, CWRES) | nlmixr2, ggplot2 |
| 7 | Clinical simulation | rxode2 |
| 8 | Deploy to GitHub Pages | docs/ |

## Getting started

```r
# 1. Open theoph_poppk_portfolio.Rproj in RStudio
# 2. First-time setup (installs packages, writes renv.lock):
source("R/00_renv_init.R")
# 3. Load helpers in any session:
source(here::here("R", "00_setup.R"))
```

### macOS Apple Silicon note (modelling stack)

`rxode2` expects gfortran at `/opt/gfortran/bin/gfortran`. If you installed it via Homebrew:

```bash
sudo mkdir -p /opt/gfortran/bin
sudo ln -sf /opt/homebrew/bin/gfortran /opt/gfortran/bin/gfortran
```

## Reproducibility

- **renv** locks all package versions (`renv.lock`); restore with `renv::restore()`.
- Snapshots are taken after each phase. The compiled stack (`rxode2`/`nlmixr2`) is excluded
  from the lockfile until a trivial model compiles successfully.
- **Quarto** renders every report to `docs/` for the published website.

---

*Educational / portfolio project using publicly available data. Not for clinical use.*
