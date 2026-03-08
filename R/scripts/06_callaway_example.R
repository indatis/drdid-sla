# ------------------------------------------------------------
# 05_callaway_santanna.R
# Callaway & Sant'Anna benchmark — SLA with DRDID
# ------------------------------------------------------------

rm(list = ls())

# ---- Packages ----
library(did)
library(DRDID)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Source repo functions ----
source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")   # loads .xlab_k, .scale_colour*, .legend_theme, etc.
source("R/functions/sla_outputs.R")

# ------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------
data("mpdta", package = "did")

# idname:  countyreal
# tname:   year
# gname:   first treatment year
# treated: I(year >= g & g > 0)
summary(mpdta)

df_long <- mpdta %>%
  mutate(treated = treat)

# Covariates (light to avoid overfitting)
covars_vec <- c("lpop")

# ------------------------------------------------------------
# 2) Filter to two-period panel (pre = 2003, post = 2007)
#    run_sla_drdid_panel() requires exactly two time points;
#    mpdta has multiple years so we subset before passing it in.
# ------------------------------------------------------------
df_long <- df_long %>%
  filter(year %in% c(2003, 2007))

# Verify exactly two periods per unit after filtering
stopifnot(all(table(df_long$countyreal) == 2))

# ------------------------------------------------------------
# 3) Run SLA (two-period subset: 2003 vs 2007)
# ------------------------------------------------------------
sla_res <- run_sla_drdid_panel(
  df          = df_long,
  id_var      = "countyreal",
  time_var    = "year",
  outcome_var = "lemp",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 2003,
  post_time   = 2007,
  min_step    = 10,
  max_steps   = NULL,
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

cat("\n[CS] Minimum IPT-stable step (strict):          ", k_strict, "\n")
cat("[CS] Minimum IPT-stable step (tolerant tol=0.05):", k_tol_05, "\n")

si_grid <- compute_SI_grid(sla_summary)

# ------------------------------------------------------------
# 5) Standard plots — all styling inherited from sla_plots.R
# ------------------------------------------------------------
save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "callaway_santanna",
  out_dir    = "outputs/figures"
)

# ------------------------------------------------------------
# 6) Combined plot — ATT + DRDID CI + GLM fallback
# ------------------------------------------------------------
p_cs_combined <- save_sla_combined_plot(
  summary_df = sla_summary,
  prefix     = "callaway_santanna",
  out_dir    = "outputs/figures",
  width      = 8,
  height     = 5,
  dpi        = 300
)

print(p_cs_combined)

# ------------------------------------------------------------
# 7) Tables
# ------------------------------------------------------------
save_sla_tables(
  sla_res     = sla_res,
  sla_summary = sla_summary,
  si_grid     = si_grid,
  prefix      = "callaway_santanna",
  out_dir     = "outputs/tables",
  k_strict    = k_strict,
  k_tol       = k_tol_05,
  tol         = 0.05,
  meta        = list(seed = 123, R = 100, trim_level = 0.995)
)