# ------------------------------------------------------------
# 01_adh_chinashock_example.R
# SLA for DRDID – ADH (China Shock) via ShiftShareSE::ADH
# ------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(DRDID)
library(ShiftShareSE)

# ---- Source your existing SLA functions ----
source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")   # loads .xlab_k, .scale_colour2, .scale_ltype2,
# .legend_theme, save_sla_combined_plot, etc.
source("R/functions/sla_outputs.R")

# ------------------------------------------------------------
# 0) Load ADH data (CZ x 2 periods)
# ------------------------------------------------------------
data("ADH", package = "ShiftShareSE")
df <- ADH$reg  # 1444 rows = 722 CZ * 2 periods

df <- ADH$reg %>%
  mutate(
    time = as.integer(t2),   # FALSE -> 0, TRUE -> 1
    id   = as.integer(czone)
  )

# ------------------------------------------------------------
# 1) Construct time-invariant treatment at CZ level
#    (top tercile of mean shock across periods)
# ------------------------------------------------------------
cz_treat <- df %>%
  group_by(id) %>%
  summarise(shock_mean = mean(shock, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    treated = as.integer(shock_mean >= quantile(shock_mean, probs = 2/3, na.rm = TRUE))
  )

# ------------------------------------------------------------
# 2) Freeze covariates at baseline (time==0)
# ------------------------------------------------------------
baseline_covars <- df %>%
  filter(time == 0) %>%
  select(
    id,
    l_shind_manuf_cbp,
    l_sh_popedu_c,
    l_sh_popfborn,
    l_sh_empl_f,
    l_sh_routine33,
    l_task_outsource
  )

# Merge treatment + baseline covars back into both periods
df_long <- df %>%
  select(id, time, d_sh_empl_mfg) %>%
  left_join(cz_treat, by = "id") %>%
  left_join(baseline_covars, by = "id") %>%
  drop_na()

# Covariates vector (time-invariant)
covars_vec <- c(
  "l_shind_manuf_cbp",
  "l_sh_popedu_c",
  "l_sh_popfborn",
  "l_sh_empl_f",
  "l_sh_routine33",
  "l_task_outsource"
)

# Quick check: treated/control counts
table(df_long$treated)

# ------------------------------------------------------------
# 3) Run SLA (two-period panel DRDID)
# ------------------------------------------------------------
treated_ids <- df_long %>% filter(treated == 1) %>% pull(id) %>% unique()
control_ids <- df_long %>% filter(treated == 0) %>% pull(id) %>% unique()
Kmax_bal    <- min(length(treated_ids), length(control_ids))

sla_res <- run_sla_drdid_panel(
  df          = df_long,
  id_var      = "id",
  time_var    = "time",
  outcome_var = "d_sh_empl_mfg",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 0,
  post_time   = 1,
  max_steps   = NULL,
  min_step    = 10,
  R           = 100,
  seed        = 123,
  trim.level  = 0.995
)

# ------------------------------------------------------------
# 4) Diagnostics
# ------------------------------------------------------------
sla_summary <- summarise_sla(sla_res)

k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

si_grid <- compute_SI_grid(sla_summary)

# ------------------------------------------------------------
# 5) Zoomed ATT plot — consistent styling with sla_plots.R
#    Uses .xlab_k, .scale_colour2, .scale_ltype2, .legend_theme
#    sourced from sla_plots.R
# ------------------------------------------------------------

has_drdid_ci <- all(c("drdid_ci_low", "drdid_ci_high") %in% names(sla_summary))

p_adh_zoom <- ggplot2::ggplot(sla_summary, ggplot2::aes(x = step)) +
  # Dispersion ribbon — grey filled, not in legend
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_ATT - sd_ATT,
                 ymax = mean_ATT + sd_ATT),
    fill = "grey50", alpha = 0.15, colour = NA
  ) +
  # Zero reference — thin grey solid, not in legend
  ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                      colour = "grey60", linewidth = 0.35) +
  # Mean ATT — black solid
  ggplot2::geom_line(
    ggplot2::aes(y = mean_ATT,
                 colour   = "Mean ATT",
                 linetype = "Mean ATT"),
    linewidth = 0.7
  ) +
  # DRDID CI — blue dotdash
  { if (has_drdid_ci) ggplot2::geom_line(
    ggplot2::aes(y = drdid_ci_low,
                 colour   = "DRDID 95% CI",
                 linetype = "DRDID 95% CI"),
    linewidth = 0.6) } +
  { if (has_drdid_ci) ggplot2::geom_line(
    ggplot2::aes(y = drdid_ci_high,
                 colour   = "DRDID 95% CI",
                 linetype = "DRDID 95% CI"),
    linewidth = 0.6) } +
  .scale_colour2 +
  .scale_ltype2  +
  ggplot2::coord_cartesian(ylim = c(-8, 8)) +
  ggplot2::labs(
    x = .xlab_k,
    y = expression(widehat(ATT))
  ) +
  ggplot2::theme_minimal() +
  .legend_theme

print(p_adh_zoom)

ggplot2::ggsave(
  filename = "outputs/figures/adh_chinashock_att_zoom.jpeg",
  plot     = p_adh_zoom,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ------------------------------------------------------------
# 6) Combined zoom plot (ATT + DRDID CI + GLM fallback, zoomed)
#    save_sla_combined_plot() returns the ggplot object invisibly,
#    so we add coord_cartesian() on top and re-save.
# ------------------------------------------------------------

p_adh_combined_zoom <- save_sla_combined_plot(
  summary_df = sla_summary,
  prefix     = "adh_chinashock_zoom",   # saves an unzoomed copy as side effect
  out_dir    = "outputs/figures",
  width      = 8,
  height     = 5,
  dpi        = 300
) +
  ggplot2::coord_cartesian(ylim = c(-8, 8))

print(p_adh_combined_zoom)

ggplot2::ggsave(
  filename = "outputs/figures/adh_chinashock_att_glm_ci_zoom.jpeg",
  plot     = p_adh_combined_zoom,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ------------------------------------------------------------
# 7) Full standard plots + tables
# ------------------------------------------------------------
save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "adh_chinashock",
  out_dir    = "outputs/figures"
)

save_sla_tables(
  sla_res     = sla_res,
  sla_summary = sla_summary,
  si_grid     = si_grid,
  prefix      = "adh_chinashock",
  out_dir     = "outputs/tables",
  k_strict    = k_strict,
  k_tol       = k_tol_05,
  tol         = 0.05,
  meta        = list(seed = 123, R = 100, trim_level = 0.995)
)

cat("\nADH China Shock — Kmax_bal =", Kmax_bal,
    "\nStrict k*IPT  =", k_strict,
    "\nTolerant k*IPT =", k_tol_05, "\n")