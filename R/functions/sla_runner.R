# ------------------------------------------------------------
# run_sla_drdid_panel()
# ------------------------------------------------------------
# Subsample Ladder Analysis (SLA) for two-period panel DRDID.
#
# Returns: data.frame with columns:
#   rep, step, nT, nC, ATT, se, ps.flag, ps.method
# ------------------------------------------------------------

decode_ps_flag <- function(flag) {
  if (is.na(flag)) return(NA_character_)
  if (flag == 0) return("IPT (trust)")
  if (flag == 1) return("IPW backup")
  if (flag == 2) return("GLM logit")
  "Unknown"
}

run_sla_drdid_panel <- function(df,
                                id_var,
                                time_var,
                                outcome_var,
                                treat_var,
                                covariates,
                                pre_time,
                                post_time,
                                treated_ids = NULL,
                                control_ids = NULL,
                                max_steps = NULL,
                                min_step = 10,
                                R = 100,
                                seed = 123,
                                trim.level = 0.995) {

  set.seed(seed)

  # ---- Infer treated / control ids if not provided ----
  if (is.null(treated_ids)) {
    treated_ids <- df |>
      dplyr::filter(.data[[treat_var]] == 1) |>
      dplyr::pull(.data[[id_var]]) |>
      unique()
  }
  if (is.null(control_ids)) {
    control_ids <- df |>
      dplyr::filter(.data[[treat_var]] == 0) |>
      dplyr::pull(.data[[id_var]]) |>
      unique()
  }

  nT <- length(treated_ids)
  nC <- length(control_ids)
  if (nT == 0 || nC == 0) stop("Need both treated and control units.")

  if (is.null(max_steps)) max_steps <- min(nT, nC)
  if (min_step > max_steps) stop("min_step cannot be larger than max_steps.")

  steps <- min_step:max_steps

  # Preallocate results list
  out_list <- vector("list", length = R * length(steps))
  idx <- 1L

  for (r in seq_len(R)) {

    treated_shuffled <- sample(treated_ids)
    control_shuffled <- sample(control_ids)

    for (k in steps) {

      t_sub <- treated_shuffled[1:k]
      c_sub <- control_shuffled[1:k]

      est <- compute_att_drdid_panel(
        df          = df,
        id_var      = id_var,
        time_var    = time_var,
        outcome_var = outcome_var,
        treat_var   = treat_var,
        covariates  = covariates,
        pre_time    = pre_time,
        post_time   = post_time,
        treated_sub = t_sub,
        control_sub = c_sub,
        trim.level  = trim.level
      )

      out_list[[idx]] <- data.frame(
        rep       = r,
        step      = k,
        nT        = k,
        nC        = k,
        ATT       = est$ATT,
        se        = est$se,
        ps.flag   = est$ps.flag,
        ps.method = decode_ps_flag(est$ps.flag),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(out_list)
}
