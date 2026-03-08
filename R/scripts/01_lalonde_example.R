# # ------------------------------------------------------------
# # 01_lalonde_example.R
# # Subsample Ladder Analysis (SLA) for DRDID - Lalonde example
# # ------------------------------------------------------------
# 
# rm(list = ls())
# 
# # ---- Packages ----
# library(DRDID)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# # ---- Source functions ----
# source("R/functions/drdid_att_panel.R")
# source("R/functions/sla_runner.R")
# source("R/functions/sla_diagnostics.R")
# source("R/functions/sla_plots.R")
# source("R/functions/sla_outputs.R")
# 
# 
# # ------------------------------------------------------------
# # 0) Build panel-long Lalonde dataset (NSW + CPS)
# # ------------------------------------------------------------
# data("nsw")
# 
# lalonde_wide <- nsw %>%
#   filter(sample %in% c(1, 2)) %>%
#   mutate(id = row_number())
# 
# lalonde_long <- lalonde_wide %>%
#   select(id, treated, re75, re78,
#          age, educ, black, married, nodegree, hisp, re74) %>%
#   pivot_longer(
#     cols      = c(re75, re78),
#     names_to  = "year",
#     values_to = "re"
#   ) %>%
#   mutate(year = ifelse(year == "re75", 1975, 1978))
# 
# covars_vec <- c("age", "educ", "black", "married", "nodegree", "hisp", "re74")
# 
# # ------------------------------------------------------------
# # 1) Run SLA
# # ------------------------------------------------------------
# sla_res <- run_sla_drdid_panel(
#   df          = lalonde_long,
#   id_var      = "id",
#   time_var    = "year",
#   outcome_var = "re",
#   treat_var   = "treated",
#   covariates  = covars_vec,
#   pre_time    = 1975,
#   post_time   = 1978,
#   max_steps   = NULL,
#   min_step    = 10,
#   R           = 100,
#   seed        = 123,
#   trim.level  = 0.995
# )
# 
# print(head(sla_res))
# print(table(sla_res$ps.method, useNA = "ifany"))
# 
# # ------------------------------------------------------------
# # 2) Diagnostics
# # ------------------------------------------------------------
# sla_summary <- summarise_sla(sla_res)
# print(head(sla_summary))
# 
# k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
# k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)
# 
# cat("\nMinimum IPT-stable step (strict): ", k_strict, "\n")
# cat("Minimum IPT-stable step (tolerant tol=0.05): ", k_tol_05, "\n")
# 
# si_grid <- compute_SI_grid(sla_summary)
# print(head(si_grid))
# 
# # ------------------------------------------------------------
# # 3) Plots (basic)
# # ------------------------------------------------------------
# 
# # Mean ATT +/- 1 SD band
# yl <- quantile(sla_summary$mean_ATT - sla_summary$sd_ATT, probs = 0.02, na.rm = TRUE)
# yu <- quantile(sla_summary$mean_ATT + sla_summary$sd_ATT, probs = 0.98, na.rm = TRUE)
# 
# p1 <- ggplot(sla_summary, aes(x = step, y = mean_ATT)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT), alpha = 0.2) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   coord_cartesian(ylim = c(yl, yu)) +
#   # labs(
#   #   title = "SLA for DRDID – Lalonde (NSW + CPS)",
#   #   x = "Subsample size (treated = controls = step)",
#   #   y = "Mean ATT across replications"
#   # ) +
#   theme_minimal()
# 
# 
# # GLM fallback rate
# p2 <- ggplot(sla_summary, aes(x = step, y = glm_pct)) +
#   geom_line() +
#   labs(
#     title = "GLM fallback rate by ladder step (ps.flag = 2)",
#     x = "Subsample size (step)",
#     y = "GLM fallback rate"
#   ) +
#   theme_minimal()
# 
# print(p2)
# 
# # Stability Index grid
# p3 <- ggplot(si_grid, aes(x = k0, y = SI)) +
#   geom_line() +
#   geom_hline(yintercept = 0.9, linetype = "dashed") +
#   labs(
#     title = "Stability Index SI(k0, K)",
#     x = "Initial k0",
#     y = "SI"
#   ) +
#   theme_minimal()
# 
# print(p3)
# 
# save_sla_plots(
#   summary_df = sla_summary,
#   si_grid    = si_grid,
#   prefix     = "lalonde",
#   out_dir    = "outputs/figures"
# )
# 
# save_sla_tables(
#   sla_res     = sla_res,
#   sla_summary = sla_summary,
#   si_grid     = si_grid,
#   prefix      = "lalonde",           # change per script
#   out_dir     = "outputs/tables",
#   k_strict    = k_strict,
#   k_tol       = k_tol_05,
#   tol         = 0.05,
#   meta = list(seed = 123, R = max(sla_res$rep), trim_level = 0.995)  # change per script
# )
# 
# p_lalonde <- save_sla_combined_plot(
#   summary_df = sla_summary,
#   prefix     = "lalonde",
#   out_dir    = "outputs/figures"
# )
# 
# p_lalonde <- p_lalonde +
#   coord_cartesian(ylim = c(-3000, 5000))
# 
# ggplot2::ggsave(
#   filename = "outputs/figures/lalonde_att_glm_ci_zoom.png",
#   plot = p_lalonde,
#   width = 8, height = 5, dpi = 300
# )
# 
# ------------------------------------------------------------
# 01_lalonde_example.R
# Subsample Ladder Analysis (SLA) for DRDID - Lalonde example
# ------------------------------------------------------------

rm(list = ls())

# ---- Packages ----
library(DRDID)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Source functions ----
source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")   # loads .xlab_k, .scale_colour2, .scale_ltype2,
# .legend_theme, save_sla_combined_plot, etc.
source("R/functions/sla_outputs.R")

# ------------------------------------------------------------
# 0) Build panel-long Lalonde dataset (NSW + CPS)
# ------------------------------------------------------------
data("nsw")

lalonde_wide <- nsw %>%
  filter(sample %in% c(1, 2)) %>%
  mutate(id = row_number())

lalonde_long <- lalonde_wide %>%
  select(id, treated, re75, re78,
         age, educ, black, married, nodegree, hisp, re74) %>%
  pivot_longer(
    cols      = c(re75, re78),
    names_to  = "year",
    values_to = "re"
  ) %>%
  mutate(year = ifelse(year == "re75", 1975, 1978))

covars_vec <- c("age", "educ", "black", "married", "nodegree", "hisp", "re74")

# ------------------------------------------------------------
# 1) Run SLA
# ------------------------------------------------------------
sla_res <- run_sla_drdid_panel(
  df          = lalonde_long,
  id_var      = "id",
  time_var    = "year",
  outcome_var = "re",
  treat_var   = "treated",
  covariates  = covars_vec,
  pre_time    = 1975,
  post_time   = 1978,
  max_steps   = NULL,
  min_step    = 10,
  R           = 100,
  seed        = 123,
  trim.level  = 0.995
)

print(head(sla_res))
print(table(sla_res$ps.method, useNA = "ifany"))

# ------------------------------------------------------------
# 2) Diagnostics
# ------------------------------------------------------------
sla_summary <- summarise_sla(sla_res)
print(head(sla_summary))

k_strict <- min_ipt_stable_step(sla_summary, method = "strict")
k_tol_05 <- min_ipt_stable_step(sla_summary, method = "tolerant", tol = 0.05)

cat("\nMinimum IPT-stable step (strict): ", k_strict, "\n")
cat("Minimum IPT-stable step (tolerant tol=0.05): ", k_tol_05, "\n")

si_grid <- compute_SI_grid(sla_summary)
print(head(si_grid))

# ------------------------------------------------------------
# 3) Standard plots — all styling inherited from sla_plots.R
# ------------------------------------------------------------
save_sla_plots(
  summary_df = sla_summary,
  si_grid    = si_grid,
  prefix     = "lalonde",
  out_dir    = "outputs/figures"
)

# ------------------------------------------------------------
# 4) Combined zoom plot — consistent styling with sla_plots.R
#    save_sla_combined_plot() returns the ggplot object invisibly;
#    coord_cartesian() is added on top for the zoomed version.
# ------------------------------------------------------------
p_lalonde <- save_sla_combined_plot(
  summary_df = sla_summary,
  prefix     = "lalonde",
  out_dir    = "outputs/figures"
) +
  ggplot2::coord_cartesian(ylim = c(-3000, 5000))

print(p_lalonde)

ggplot2::ggsave(
  filename = "outputs/figures/lalonde_att_glm_ci_zoom.jpeg",
  plot     = p_lalonde,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ------------------------------------------------------------
# 5) Tables
# ------------------------------------------------------------
save_sla_tables(
  sla_res     = sla_res,
  sla_summary = sla_summary,
  si_grid     = si_grid,
  prefix      = "lalonde",
  out_dir     = "outputs/tables",
  k_strict    = k_strict,
  k_tol       = k_tol_05,
  tol         = 0.05,
  meta        = list(seed = 123, R = max(sla_res$rep), trim_level = 0.995)
)