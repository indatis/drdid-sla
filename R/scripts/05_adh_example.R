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
source("R/functions/sla_plots.R")
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
#    because your pivot_wider uses covariates as id_cols
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
  select(id, time, d_sh_empl_mfg) %>%   # choose outcome here
  left_join(cz_treat, by = "id") %>%
  left_join(baseline_covars, by = "id") %>%
  drop_na()

# Covariates vector (time-invariant now)
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
# Max feasible balanced step:
treated_ids <- df_long %>% filter(treated == 1) %>% pull(id) %>% unique()
control_ids <- df_long %>% filter(treated == 0) %>% pull(id) %>% unique()
Kmax_bal <- min(length(treated_ids), length(control_ids))

sla_res <- run_sla_drdid_panel(
  df          = df_long,
  id_var      = "id",
  time_var    = "time",
  outcome_var = "d_sh_empl_mfg",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 0,
  post_time   = 1,
  max_steps   = NULL,   # or set a smaller cap like 200 for speed
  min_step    = 10,         # match your default
  R           = 100,
  seed        = 123,
  trim.level  = 0.995
)

# ------------------------------------------------------------
# 4) Diagnostics + plots + outputs
# ------------------------------------------------------------
sla_summary <- summarise_sla(sla_res)
# ------------------------------------------------------------
# 5) Zoomed ATT plot (ADH only)
# ------------------------------------------------------------

p_adh_zoom <- ggplot(sla_summary, aes(x = step, y = mean_ATT)) +
  geom_ribbon(
    aes(ymin = mean_ATT - sd_ATT,
        ymax = mean_ATT + sd_ATT),
    alpha = 0.2
  ) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-8, 8)) +
  labs(
    title = "ADH China Shock — SLA (zoomed view)",
    subtitle = "ATT trajectory restricted to [−10, 10]",
    x = "Balanced subsample size k",
    y = "ATT"
  ) +
  theme_minimal()

print(p_adh_zoom)

ggsave(
  filename = "outputs/figures/adh_chinashock_att_zoom.png",
  plot = p_adh_zoom,
  width = 8,
  height = 5,
  dpi = 300
)
###
p_adh_combined_zoom <- save_sla_combined_plot(
  summary_df = sla_summary,
  prefix     = "adh_chinashock_zoom",
  out_dir    = "outputs/figures",
  width      = 8,
  height     = 5,
  dpi        = 300
) +
  coord_cartesian(ylim = c(-8, 8))

print(p_adh_combined_zoom)

ggsave(
  filename = "outputs/figures/adh_chinashock_att_glm_ci_zoom.png",
  plot = p_adh_combined_zoom,
  width = 8,
  height = 5,
  dpi = 300
)

###
k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

si_grid <- compute_SI_grid(sla_summary)

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
    "\nStrict k*IPT =", k_strict,
    "\nTolerant k*IPT =", k_tol_05, "\n")
