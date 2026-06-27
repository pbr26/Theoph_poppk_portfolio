# =============================================================================
# 01_eda.R
# Purpose : Helper functions for Exploratory Data Analysis of Theoph dataset
# Author  : Pramod BR
# Date    : 2026-06-21
# =============================================================================

# ── 1. Load and prepare Theoph dataset ───────────────────────────────────────

#' Prepare the Theoph dataset for PK analysis
#'
#' Returns a cleaned tibble with renamed columns following NONMEM convention:
#'   ID, TIME, DV, AMT, WT, DOSE
#'
prepare_theoph <- function() {
  data(Theoph, package = "datasets")

  df <- Theoph |>
    as_tibble() |>
    rename(
      ID   = Subject,   # Subject identifier (factor → numeric)
      TIME = Time,      # Time after dose (hours)
      DV   = conc,      # Observed concentration (mg/L)
      WT   = Wt,        # Body weight (kg)
      DOSE = Dose       # Dose (mg/kg)
    ) |>
    mutate(
      ID     = as.numeric(as.character(ID)),
      AMT    = ifelse(TIME == 0, DOSE * WT, 0),  # Actual dose in mg
      EVID   = ifelse(TIME == 0, 1, 0),           # Event ID: 1=dose, 0=obs
      MDV    = ifelse(TIME == 0, 1, 0),           # Missing DV flag
      LLOQ   = 0.1                                # Assumed LLOQ (mg/L)
    ) |>
    arrange(ID, TIME)

  return(df)
}


# ── 2. Data quality checks ────────────────────────────────────────────────────

#' Run data quality checks and return a summary list
run_dq_checks <- function(df) {
  obs <- df |> filter(EVID == 0)   # observations only

  list(
    n_subjects     = n_distinct(df$ID),
    n_obs_total    = nrow(obs),
    n_obs_per_subj = obs |> count(ID) |> pull(n) |> summary(),
    time_range     = range(obs$TIME),
    conc_range     = range(obs$DV),
    n_below_lloq   = sum(obs$DV < obs$LLOQ),
    n_missing_conc = sum(is.na(obs$DV)),
    n_missing_time = sum(is.na(obs$TIME)),
    weight_range   = range(df$WT),
    dose_range     = range(df$DOSE),
    dose_mg_range  = range(df$AMT[df$EVID == 1])
  )
}


# ── 3. Subject-level summary ──────────────────────────────────────────────────

#' Compute per-subject descriptive PK summary
subject_summary <- function(df) {
  df |>
    filter(EVID == 0) |>
    group_by(ID, WT, DOSE) |>
    summarise(
      AMT_mg   = unique(DOSE * WT),
      Cmax     = max(DV, na.rm = TRUE),
      Tmax     = TIME[which.max(DV)],
      Clast    = DV[which.max(TIME)],
      Tlast    = max(TIME),
      n_obs    = n(),
      .groups  = "drop"
    ) |>
    mutate(across(where(is.numeric), \(x) round(x, 2)))
}


# ── 4. Plotting functions ─────────────────────────────────────────────────────

#' Individual concentration-time profiles (linear scale)
plot_individual_linear <- function(df) {
  obs <- df |> filter(EVID == 0)

  ggplot(obs, aes(x = TIME, y = DV, group = ID, colour = factor(ID))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ID, ncol = 4, labeller = label_both) +
    scale_colour_viridis_d(option = "turbo", guide = "none") +
    scale_x_continuous(breaks = seq(0, 25, 5)) +
    labs(
      title    = "Individual Concentration-Time Profiles",
      subtitle = "Theophylline – oral single dose, n = 12 subjects",
      x        = "Time after dose (h)",
      y        = "Theophylline concentration (mg/L)",
      caption  = "Each panel = one subject. Dose administered at time 0."
    )
}


#' Individual concentration-time profiles (semi-log scale)
plot_individual_semilog <- function(df) {
  obs <- df |> filter(EVID == 0, DV > 0)

  ggplot(obs, aes(x = TIME, y = DV, group = ID, colour = factor(ID))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ID, ncol = 4, labeller = label_both) +
    scale_y_log10(
      breaks = c(0.1, 0.5, 1, 2, 5, 10, 20),
      labels = scales::label_number()
    ) +
    scale_colour_viridis_d(option = "turbo", guide = "none") +
    scale_x_continuous(breaks = seq(0, 25, 5)) +
    annotation_logticks(sides = "l", colour = "grey60", linewidth = 0.3) +
    labs(
      title    = "Individual Concentration-Time Profiles (Semi-log)",
      subtitle = "Log-linear terminal phase reveals first-order elimination",
      x        = "Time after dose (h)",
      y        = "Theophylline concentration (mg/L) – log scale",
      caption  = "Straight line on semi-log = first-order elimination kinetics"
    )
}


#' Mean ± SD concentration-time profile (all subjects overlaid + mean)
plot_mean_profile <- function(df) {
  obs <- df |> filter(EVID == 0)

  mean_df <- obs |>
    group_by(TIME) |>
    summarise(
      mean_conc = mean(DV, na.rm = TRUE),
      sd_conc   = sd(DV, na.rm = TRUE),
      n         = n(),
      se_conc   = sd_conc / sqrt(n),
      .groups   = "drop"
    )

  ggplot() +
    # Individual lines (faded)
    geom_line(
      data    = obs,
      aes(x = TIME, y = DV, group = ID),
      colour  = "steelblue", alpha = 0.25, linewidth = 0.5
    ) +
    geom_point(
      data    = obs,
      aes(x = TIME, y = DV),
      colour  = "steelblue", alpha = 0.3, size = 1.5
    ) +
    # Mean ± SD ribbon
    geom_ribbon(
      data    = mean_df,
      aes(x = TIME, ymin = mean_conc - sd_conc, ymax = mean_conc + sd_conc),
      fill    = "#E74C3C", alpha = 0.15
    ) +
    # Mean line
    geom_line(
      data      = mean_df,
      aes(x = TIME, y = mean_conc),
      colour    = "#E74C3C", linewidth = 1.2
    ) +
    geom_point(
      data    = mean_df,
      aes(x = TIME, y = mean_conc),
      colour  = "#E74C3C", size = 3, shape = 18
    ) +
    scale_x_continuous(breaks = seq(0, 25, 5)) +
    labs(
      title    = "Mean (±SD) Concentration-Time Profile",
      subtitle = "Red = mean ± SD; blue = individual observations (n = 12)",
      x        = "Time after dose (h)",
      y        = "Theophylline concentration (mg/L)",
      caption  = "Absorption phase visible in first 2 h; elimination phase from ~2–25 h"
    )
}


#' All subjects overlaid on one panel (spaghetti plot)
plot_spaghetti <- function(df) {
  obs <- df |> filter(EVID == 0)

  ggplot(obs, aes(x = TIME, y = DV, group = ID, colour = factor(ID))) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    geom_point(size = 2, alpha = 0.85) +
    scale_colour_viridis_d(option = "turbo", name = "Subject ID") +
    scale_x_continuous(breaks = seq(0, 25, 5)) +
    labs(
      title    = "Spaghetti Plot – All Subjects",
      subtitle = "Between-subject variability in Cmax, Tmax, and elimination",
      x        = "Time after dose (h)",
      y        = "Theophylline concentration (mg/L)"
    ) +
    guides(colour = guide_legend(nrow = 2))
}


#' Weight distribution plot
plot_weight_dist <- function(df) {
  subj <- df |> distinct(ID, WT, DOSE) |>
    mutate(AMT_mg = DOSE * WT)

  p1 <- ggplot(subj, aes(x = WT)) +
    geom_histogram(bins = 8, fill = "#2980B9", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = mean(subj$WT), linetype = "dashed", colour = "#E74C3C") +
    labs(title = "Body Weight Distribution", x = "Weight (kg)", y = "Count") +
    annotate("text", x = mean(subj$WT) + 1.5, y = 2.8,
             label = paste0("Mean = ", round(mean(subj$WT), 1), " kg"),
             colour = "#E74C3C", size = 3.5)

  p2 <- ggplot(subj, aes(x = AMT_mg)) +
    geom_histogram(bins = 8, fill = "#27AE60", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = mean(subj$AMT_mg), linetype = "dashed", colour = "#E74C3C") +
    labs(title = "Total Dose Distribution", x = "Dose (mg)", y = "Count") +
    annotate("text", x = mean(subj$AMT_mg) + 15, y = 2.8,
             label = paste0("Mean = ", round(mean(subj$AMT_mg), 1), " mg"),
             colour = "#E74C3C", size = 3.5)

  p1 + p2 +
    plot_annotation(
      title    = "Subject Characteristics",
      subtitle = "Weight-normalised dosing (mg/kg) leads to variable absolute doses"
    )
}


#' Cmax vs Body Weight scatter
plot_cmax_vs_wt <- function(df) {
  ss <- subject_summary(df)

  ggplot(ss, aes(x = WT, y = Cmax, label = ID)) +
    geom_point(aes(size = AMT_mg), colour = "#8E44AD", alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, colour = "#E74C3C",
                linetype = "dashed", linewidth = 0.8) +
    ggrepel::geom_text_repel(size = 3, colour = "grey30") +
    scale_size_continuous(name = "Dose (mg)", range = c(3, 8)) +
    labs(
      title    = "Cmax vs Body Weight",
      subtitle = "Exploring covariate relationships before formal modeling",
      x        = "Body weight (kg)",
      y        = "Observed Cmax (mg/L)",
      caption  = "Point size ∝ total dose administered"
    )
}
