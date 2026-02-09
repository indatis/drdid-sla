# ------------------------------------------------------------
# compute_att_drdid_panel()
# ------------------------------------------------------------
# Computes ATT using DRDID (improved estimator) for
# two-period panel data on a selected subsample.
#
# Returns:
#   list(ATT, se, ci_low, ci_high, ps.flag)
# ------------------------------------------------------------

compute_att_drdid_panel <- function(df,
                                    id_var,
                                    time_var,
                                    outcome_var,
                                    treat_var,
                                    covariates,
                                    pre_time,
                                    post_time,
                                    treated_sub,
                                    control_sub,
                                    trim.level = 0.995,
                                    z = 1.96) {
  
  # ---- Input checks ----
  if (length(treated_sub) == 0 || length(control_sub) == 0) {
    return(list(
      ATT = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_,
      ps.flag = NA_integer_
    ))
  }
  
  # ---- Subset data ----
  subdf <- df |>
    dplyr::filter(.data[[id_var]] %in% c(treated_sub, control_sub)) |>
    dplyr::select(dplyr::all_of(c(id_var, time_var, treat_var, outcome_var, covariates)))
  
  # ---- Convert to wide format ----
  wide <- subdf |>
    tidyr::pivot_wider(
      id_cols      = c(dplyr::all_of(id_var),
                       dplyr::all_of(treat_var),
                       dplyr::all_of(covariates)),
      names_from   = dplyr::all_of(time_var),
      values_from  = dplyr::all_of(outcome_var),
      names_prefix = "y"
    )
  
  y0_name <- paste0("y", pre_time)
  y1_name <- paste0("y", post_time)
  
  if (!all(c(y0_name, y1_name) %in% names(wide))) {
    return(list(
      ATT = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_,
      ps.flag = NA_integer_
    ))
  }
  
  # ---- Extract variables ----
  y0 <- wide[[y0_name]]
  y1 <- wide[[y1_name]]
  D  <- wide[[treat_var]]
  
  X <- wide |>
    dplyr::select(dplyr::all_of(covariates)) |>
    as.matrix()
  
  # ---- DRDID estimation ----
  out <- tryCatch(
    DRDID::drdid_imp_panel(
      y1 = y1,
      y0 = y0,
      D  = D,
      covariates = X,
      trim.level = trim.level
    ),
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(list(
      ATT = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_,
      ps.flag = NA_integer_
    ))
  }
  
  # ---- Extract + compute CI ----
  att <- as.numeric(out$ATT)
  se  <- as.numeric(out$se)
  
  # If se is missing/non-finite, CI is NA
  if (!is.finite(se) || is.na(se) || se <= 0) {
    ci_low  <- NA_real_
    ci_high <- NA_real_
  } else {
    ci_low  <- att - z * se
    ci_high <- att + z * se
  }
  
  list(
    ATT     = att,
    se      = se,
    ci_low  = ci_low,
    ci_high = ci_high,
    ps.flag = out$ps.flag
  )
}
