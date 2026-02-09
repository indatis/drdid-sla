# ------------------------------------------------------------
# SLA diagnostics: summaries, IPT-stability, Stability Index
# ------------------------------------------------------------

# ---- Step-level summary table ----
summarise_sla <- function(sla_res, z = 1.96) {
  
  # 1) ATT stability summaries across replications
  att_sum <- sla_res |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      mean_ATT = mean(ATT, na.rm = TRUE),
      sd_ATT   = stats::sd(ATT, na.rm = TRUE),
      n_ok     = sum(!is.na(ATT)),
      
      # 2) Aggregated DRDID analytical 95% CI across replications
      # (robust aggregation: medians of the CI endpoints)
      drdid_ci_low  = stats::median(ci_low,  na.rm = TRUE),
      drdid_ci_high = stats::median(ci_high, na.rm = TRUE),
      
      # Optional: Monte Carlo CI for the mean trajectory across replications
      # (NOT a DRDID inferential CI; this is about precision of mean_ATT)
      se_mean = dplyr::if_else(n_ok > 0, sd_ATT / sqrt(n_ok), NA_real_),
      mc_ci_low  = mean_ATT - z * se_mean,
      mc_ci_high = mean_ATT + z * se_mean,
      
      .groups  = "drop"
    )
  
  # 3) Significance rate across replications (CI excludes 0)
  sig_sum <- sla_res |>
    dplyr::mutate(
      sig_95 = dplyr::if_else(
        !is.na(ci_low) & !is.na(ci_high),
        (ci_low > 0) | (ci_high < 0),
        NA
      )
    ) |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      sig_pct = mean(sig_95, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4) Propensity score diagnostic rates
  ps_sum <- sla_res |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      ipt_pct = mean(ps.flag == 0, na.rm = TRUE),
      ipw_pct = mean(ps.flag == 1, na.rm = TRUE),
      glm_pct = mean(ps.flag == 2, na.rm = TRUE),
      .groups = "drop"
    )
  
  att_sum |>
    dplyr::left_join(sig_sum, by = "step") |>
    dplyr::left_join(ps_sum,  by = "step")
}



# ---- Minimum IPT-stable step ----
# strict: glm_pct == 0 for all steps >= k
# tolerant: glm_pct <= tol for all steps >= k
min_ipt_stable_step <- function(summary_df,
                                method = c("strict", "tolerant"),
                                tol = 0.05) {

  method <- match.arg(method)
  if (!("glm_pct" %in% names(summary_df))) {
    stop("summary_df must contain glm_pct. Use summarise_sla() first.")
  }

  summary_df <- summary_df |>
    dplyr::arrange(step)

  steps <- summary_df$step

  cond_fun <- function(k) {
    tail_vals <- summary_df$glm_pct[summary_df$step >= k]
    if (length(tail_vals) == 0) return(FALSE)
    if (method == "strict") return(all(tail_vals == 0))
    if (method == "tolerant") return(all(tail_vals <= tol))
    FALSE
  }

  k_star <- steps[sapply(steps, cond_fun)]
  if (length(k_star) == 0) return(NA_integer_)
  k_star[1]
}

# ---- Stability Index (SI) ----
compute_SI <- function(summary_df,
                       k0,
                       K = max(summary_df$step),
                       eps = 1e-8) {

  sub <- summary_df |>
    dplyr::filter(step >= k0, step <= K, !is.na(mean_ATT))

  if (nrow(sub) == 0) return(NA_real_)

  ATT_bar <- mean(sub$mean_ATT)
  dev_rel <- abs(sub$mean_ATT - ATT_bar) / (abs(ATT_bar) + eps)

  1 - mean(dev_rel)
}

compute_SI_grid <- function(summary_df,
                            K = max(summary_df$step),
                            eps = 1e-8) {
  
  min_step <- min(summary_df$step, na.rm = TRUE)
  max_step <- max(summary_df$step, na.rm = TRUE)
  
  data.frame(
    k0 = min_step:max_step,
    SI = sapply(min_step:max_step, function(k) compute_SI(summary_df, k0 = k, K = K, eps = eps))
  )
}
