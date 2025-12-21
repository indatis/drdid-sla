# ------------------------------------------------------------
# 02_simulation_strong.R
# SLA for DRDID – simulation with strong treatment effect
# ------------------------------------------------------------

rm(list = ls())

library(DRDID)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Source functions ----
source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")
source("R/functions/simulate_panel.R")
source("R/functions/sla_outputs.R")


# ---- Generate data (strong effect) ----
df_long <- simulate_two_period_panel(
  N = 100,
  N_T = 50,
  true_ATT = 20,
  noise_sd = 2,
  seed = 123
)

# ---- Benchmark DRDID ----
drdid_res <- DRDID::drdid(
  yname    = "Y",
  tname    = "time",
  idname   = "id",
  dname    = "treated",
  xformla  = ~ X1 + X2,
  data     = df_long,
  panel    = TRUE,
  estMethod = "imp"
)
print(drdid_res)

# ---- SLA ----
covars_vec <- c("X1", "X2")

sla_res <- run_sla_drdid_panel(
  df          = df_long,
  id_var      = "id",
  time_var    = "time",
  outcome_var = "Y",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 0,
  post_time   = 1,
  max_steps   = 50,
  min_step    = 1,
  R           = 20,
  seed        = 123,
  trim.level  = 0.995
)

# ---- Diagnostics ----
sla_summary <- summarise_sla(sla_res)

k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

cat("\n[Strong] Minimum IPT-stable step (strict): ", k_strict, "\n")
cat("[Strong] Minimum IPT-stable step (tolerant tol=0.05): ", k_tol_05, "\n")

si_grid <- compute_SI_grid(sla_summary)

# ---- Save plots ----
save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "sim_strong",
  out_dir    = "outputs/figures"
)
