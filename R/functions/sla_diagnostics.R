# ------------------------------------------------------------
# SLA diagnostics: summaries, IPT-stability, Stability Index
# ------------------------------------------------------------

# ---- Step-level summary table ----
summarise_sla <- function(sla_res) {

  # Basic ATT summary
  att_sum <- sla_res |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      mean_ATT = mean(ATT, na.rm = TRUE),
      sd_ATT   = stats::sd(ATT, na.rm = TRUE),
      n_ok     = sum(!is.na(ATT)),
      .groups  = "drop"
    )

  # Propensity score diagnostic rates
  ps_sum <- sla_res |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      ipt_pct = mean(ps.flag == 0, na.rm = TRUE),
      ipw_pct = mean(ps.flag == 1, na.rm = TRUE),
      glm_pct = mean(ps.flag == 2, na.rm = TRUE),
      .groups = "drop"
    )

  dplyr::left_join(att_sum, ps_sum, by = "step")
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

  max_step <- max(summary_df$step)
  data.frame(
    k0 = 1:max_step,
    SI = sapply(1:max_step, function(k) compute_SI(summary_df, k0 = k, K = K, eps = eps))
  )
}
