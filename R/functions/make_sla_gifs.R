# ============================================================
# Script: make_sla_gifs.R
# Purpose:
#   Create animated GIFs showing how SLA summaries change as the
#   number of replications R increases (e.g., 10 → 25 → 50 → 75 → 100).
#
# Inputs:
#   - Raw SLA outputs saved by your pipeline (outputs/tables/*_sla_raw.csv)
#     Expected columns (minimum):
#       step, rep, ATT, ps.flag
#     Optional (if you want blue dashed lines):
#       ci_low, ci_high
#
# Outputs:
#   - GIF files saved to outputs/figures/
#
# Notes:
#   - This script avoids forcing package upgrades.
#   - If gganimate is not available due to ggplot2 version constraints,
#     it will stop with a clear message.
# ============================================================

rm(list = ls())

# ---- Package setup (do NOT auto-upgrade ggplot2) ----
pkgs <- c("dplyr", "ggplot2", "gifski", "transformr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gganimate)
  library(gifski)
  library(transformr)
})
remotes::install_version( "gganimate", version = "1.0.8", repos = "https://cloud.r-project.org" )

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

.make_summary_for_R <- function(sla_res, R_use) {
  sub <- sla_res %>% dplyr::filter(rep <= R_use)
  
  att_sum <- sub %>%
    group_by(step) %>%
    summarise(
      mean_ATT = mean(ATT, na.rm = TRUE),
      sd_ATT   = sd(ATT,  na.rm = TRUE),
      .groups  = "drop"
    )
  
  # CI is optional (only if columns exist in raw)
  has_ci <- all(c("ci_low", "ci_high") %in% names(sub))
  if (has_ci) {
    ci_sum <- sub %>%
      group_by(step) %>%
      summarise(
        drdid_ci_low  = mean(ci_low,  na.rm = TRUE),
        drdid_ci_high = mean(ci_high, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    ci_sum <- att_sum %>%
      transmute(step, drdid_ci_low = NA_real_, drdid_ci_high = NA_real_)
  }
  
  ps_sum <- sub %>%
    group_by(step) %>%
    summarise(
      glm_pct = mean(ps.flag == 2, na.rm = TRUE),
      .groups = "drop"
    )
  
  att_sum %>%
    left_join(ci_sum, by = "step") %>%
    left_join(ps_sum, by = "step") %>%
    mutate(R_use = R_use)
}

make_sla_gif <- function(raw_path,
                         out_gif,
                         R_grid = c(10, 25, 50, 75, 100),
                         ylim = NULL,
                         title_prefix = "SLA sensitivity to number of replications",
                         nframes = 80,
                         fps = 10,
                         width = 800,
                         height = 500) {
  
  if (!file.exists(raw_path)) stop("File not found: ", raw_path)
  
  sla_res <- read.csv(raw_path)
  
  # Required columns check
  req <- c("step", "rep", "ATT", "ps.flag")
  miss <- setdiff(req, names(sla_res))
  if (length(miss) > 0) {
    stop("Missing required columns in ", raw_path, ": ", paste(miss, collapse = ", "))
  }
  
  # Ensure numeric where needed
  sla_res <- sla_res %>%
    mutate(
      step = as.integer(step),
      rep  = as.integer(rep),
      ATT  = as.numeric(ATT)
    )
  
  summary_all <- bind_rows(lapply(R_grid, function(Ru) .make_summary_for_R(sla_res, Ru)))
  
  # Scale factor to overlay glm_pct (0-1) on ATT axis
  scale_factor <- stats::quantile(abs(summary_all$mean_ATT), probs = 0.98, na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
  
  has_ci <- any(!is.na(summary_all$drdid_ci_low)) && any(!is.na(summary_all$drdid_ci_high))
  
  p <- ggplot(summary_all, aes(x = step)) +
    geom_ribbon(aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT), alpha = 0.15) +
    geom_line(aes(y = mean_ATT)) +
    { if (has_ci) geom_line(aes(y = drdid_ci_low),  linetype = "dashed", colour = "blue") } +
    { if (has_ci) geom_line(aes(y = drdid_ci_high), linetype = "dashed", colour = "blue") } +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(aes(y = glm_pct * scale_factor), colour = "red") +
    scale_y_continuous(
      name = "ATT",
      sec.axis = sec_axis(~ . / scale_factor, name = "GLM fallback rate (ps.flag = 2)")
    ) +
    labs(
      title = paste0(title_prefix, " (R = {closest_state})"),
      x = "k (balanced treated/controls)",
      y = "ATT"
    ) +
    theme_minimal() +
    transition_states(R_use, state_length = 1) +
    ease_aes("linear")
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  dir.create(dirname(out_gif), recursive = TRUE, showWarnings = FALSE)
  
  animate(
    p,
    nframes = nframes, fps = fps, width = width, height = height,
    renderer = gifski_renderer(out_gif)
  )
}

# ------------------------------------------------------------
# RUNS
# ------------------------------------------------------------

# Zoomed datasets
make_sla_gif(
  raw_path = "outputs/tables/adh_chinashock_sla_raw.csv",
  out_gif  = "outputs/figures/adh_chinashock_sweep.gif",
  ylim     = c(-8, 8),
  title_prefix = "ADH China Shock — SLA sensitivity to R"
)

make_sla_gif(
  raw_path = "outputs/tables/lalonde_sla_raw.csv",
  out_gif  = "outputs/figures/lalonde_sweep.gif",
  ylim     = c(-4000, 5000),
  title_prefix = "Lalonde (NSW+CPS) — SLA sensitivity to R"
)

# Others (no zoom)
make_sla_gif(
  raw_path = "outputs/tables/callaway_santanna_sla_raw.csv",
  out_gif  = "outputs/figures/callaway_santanna_sweep.gif",
  ylim     = NULL,
  title_prefix = "Callaway–Sant’Anna — SLA sensitivity to R"
)

make_sla_gif(
  raw_path = "outputs/tables/card_krueger_fastfood_sla_raw.csv",
  out_gif  = "outputs/figures/card_krueger_fastfood_sweep.gif",
  ylim     = NULL,
  title_prefix = "Card & Krueger (Fast Food) — SLA sensitivity to R"
)

make_sla_gif(
  raw_path = "outputs/tables/null_sla_raw.csv",
  out_gif  = "outputs/figures/null_sweep.gif",
  ylim     = NULL,
  title_prefix = "Null Effect — SLA sensitivity to R"
)

make_sla_gif(
  raw_path = "outputs/tables/strong_sla_raw.csv",
  out_gif  = "outputs/figures/strong_sweep.gif",
  ylim     = NULL,
  title_prefix = "Strong Effect — SLA sensitivity to R"
)
