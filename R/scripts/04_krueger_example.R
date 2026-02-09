############################################################
# 04_card_krueger_fastfood.R
# Card & Krueger (1994)-style DiD using Fastfood data (loedata)
# Subsample Ladder Analysis (SLA) for DRDID – Panel Data
#
# Design choice:
# - Use a parsimonious, standard covariate set commonly used in replications:
#   chain indicators + ownershiors.
# - Exclude "prices", "hours", and other operational variables that are often
#   time-varying and/or missing, and are not necessary for the canonical DiD.
############################################################

rm(list = ls())
set.seed(123)

# ---- Packages ----
pkgs <- c("DRDID", "dplyr", "tidyr", "ggplot2")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

library(DRDID)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Source repo functions ----
source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")
source("R/functions/simulate_panel.R")
source("R/functions/sla_outputs.R")

############################################################
# 1) Load data (preferred: loedata::Fastfood; fallback: evident::ck)
############################################################

dat <- NULL
dat_name <- NULL

if (requireNamespace("loedata", quietly = TRUE)) {
  data("Fastfood", package = "loedata")
  dat <- get("Fastfood")
  dat_name <- "Fastfood (loedata)"
} else if (requireNamespace("evident", quietly = TRUE)) {
  data("ck", package = "evident")
  dat <- get("ck")
  dat_name <- "ck (evident)"
} else {
  stop(
    "Could not find a packaged Card–Krueger dataset.\n",
    "Install one of:\n",
    "  install.packages('loedata')  # preferred\n",
    "  install.packages('evident')  # fallback\n"
  )
}

message("Using dataset: ", dat_name)
dat <- as.data.frame(dat)

stopifnot(all(c("id","after","nj","fte") %in% names(dat)))

############################################################
# 2) Balanced two-period panel
############################################################

dat <- dat %>% filter(after %in% c(0, 1))

ids_bal <- dat %>%
  distinct(id, after) %>%
  count(id) %>%
  filter(n == 2) %>%
  pull(id)

dat_bal <- dat %>% filter(id %in% ids_bal)

############################################################
# 3) Build long panel and baseline-fix covariates (panel requirement)
############################################################

# Parsimonious, standard covariate set:
# - chain: Burger King / KFC / Roy Rogers / Wendy’s (categorical)
# - co_owned: ownership structure
covars_raw <- c("chain","co_owned")

df_long_raw <- dat_bal %>%
  transmute(
    id      = as.integer(id),
    time    = as.integer(after),  # 0=pre, 1=post
    treated = as.integer(nj),     # 1=NJ, 0=PA
    Y       = as.numeric(fte),
    
    chain    = chain,
    co_owned = co_owned,
    southj   = southj,
    centralj = centralj,
    northj   = northj,
    shore    = shore,
    pa1      = pa1,
    pa2      = pa2
  )

# Drop missing rows for required fields (Y + covariates)
needed <- c("id","time","treated","Y", covars_raw)
df_long_cc <- df_long_raw %>%
  select(all_of(needed)) %>%
  drop_na()

# Re-enforce balanced panel after dropping NAs
ids_bal2 <- df_long_cc %>%
  distinct(id, time) %>%
  count(id) %>%
  filter(n == 2) %>%
  pull(id)

df_long_cc <- df_long_cc %>% filter(id %in% ids_bal2)

# Baseline covariates from pre only, attached to both periods
base_cov <- df_long_cc %>%
  filter(time == 0) %>%
  select(id, all_of(covars_raw))

df_long_fixed <- df_long_cc %>%
  select(id, time, treated, Y) %>%
  left_join(base_cov, by = "id")

stopifnot(!anyNA(df_long_fixed[, covars_raw]))

############################################################
# 4) Dummy encode covariates (numeric design matrix)
############################################################

X_mm <- model.matrix(~ . - 1, data = df_long_fixed %>% select(all_of(covars_raw)))
X_mm <- as.data.frame(X_mm)
names(X_mm) <- make.names(names(X_mm))

df_long <- df_long_fixed %>%
  select(id, time, treated, Y) %>%
  bind_cols(X_mm)

covars_vec <- names(X_mm)
xformla <- as.formula(paste("~", paste(covars_vec, collapse = " + ")))

# Basic checks
stopifnot(all(covars_vec %in% names(df_long)))
stopifnot(all(df_long$time %in% c(0, 1)))
stopifnot(range(table(df_long$id))[1] == 2 && range(table(df_long$id))[2] == 2)

############################################################
# 5) Baseline DRDID estimate (panel = TRUE)
############################################################

drdid_res <- DRDID::drdid(
  yname   = "Y",
  tname   = "time",
  idname  = "id",
  dname   = "treated",
  xformla = xformla,
  data    = df_long,
  panel   = TRUE,
  estMethod = "imp"
)

print(drdid_res)

############################################################
# 6) Run SLA (balanced k treated + k controls at each step)
############################################################

sla_res <- run_sla_drdid_panel(
  df          = df_long,
  id_var      = "id",
  time_var    = "time",
  outcome_var = "Y",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 0,
  post_time   = 1,
  min_step    = 10,#1
  max_steps   = NULL,  # min(nT, nC)
  R           = 100, #20
  seed        = 123,
  trim.level  = 0.995
)

############################################################
# 7) Summaries + plots (ATT band, GLM fallback, SI grid, combined)
############################################################

sla_summary <- summarise_sla(sla_res)   # <- this is the key line

max_step <- max(sla_summary$step, na.rm = TRUE)

si_grid <- data.frame(
  k0 = 1:max_step,
  SI = sapply(1:max_step, function(k) compute_SI(sla_summary, k0 = k, K = max_step))
)

save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "card_krueger_fastfood",
  out_dir    = "outputs/figures"
)

k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

cat("\n[Null] Minimum IPT-stable step (strict): ", k_strict, "\n")
cat("[Null] Minimum IPT-stable step (tolerant tol=0.05): ", k_tol_05, "\n")

si_grid <- compute_SI_grid(sla_summary)


save_sla_tables(
  sla_res     = sla_res,
  sla_summary = sla_summary,
  si_grid     = si_grid,
  prefix      = "card_krueger_fastfood",
  out_dir     = "outputs/tables",
  k_strict = k_strict,
  k_tol       = k_tol_05,
  tol         = 0.05,
  meta = list(seed = 123, R = max(sla_res$rep), trim_level = 0.995)
)

