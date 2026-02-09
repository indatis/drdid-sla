# ------------------------------------------------------------
# Callaway–Sant’Anna benchmark — SLA with DRDID
# ------------------------------------------------------------

library(did)
library(DRDID)
library(dplyr)
library(tidyr)

# Load data
data("mpdta", package = "did")

# Basic structure
# idname: id
# tname: year
# gname: first treatment year
# treated: I(year >= g & g > 0)
summary(mpdta)
df_long <- mpdta %>%
  mutate(
    treated = treat
  )

# Covariates (keep it light to avoid overfitting)
covars_vec <- c("lpop")

# ---- Run SLA ----
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

# ---- Diagnostics ----
sla_summary <- summarise_sla(sla_res)

k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

si_grid <- compute_SI_grid(sla_summary)

print(k_strict)
print(k_tol_05)

# ---- Save outputs ----
save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "callaway_santanna",
  out_dir    = "outputs/figures"
)

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
