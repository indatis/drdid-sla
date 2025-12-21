# ------------------------------------------------------------
# Output utilities: save metrics and summaries to CSV
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

save_sla_tables <- function(sla_res,
                            sla_summary,
                            si_grid,
                            prefix,
                            out_dir = "outputs/tables",
                            k_strict = NA_integer_,
                            k_tol = NA_integer_,
                            tol = 0.05,
                            meta = list()) {

  ensure_dir(out_dir)

  # 1) Save full step-level summary
  utils::write.csv(
    sla_summary,
    file = file.path(out_dir, paste0(prefix, "_sla_summary.csv")),
    row.names = FALSE
  )

  # 2) Save raw replication-level results (optional but valuable)
  utils::write.csv(
    sla_res,
    file = file.path(out_dir, paste0(prefix, "_sla_raw.csv")),
    row.names = FALSE
  )

  # 3) Save SI grid
  utils::write.csv(
    si_grid,
    file = file.path(out_dir, paste0(prefix, "_si_grid.csv")),
    row.names = FALSE
  )

  # 4) Save a compact “paper metrics” table (single row)
  Kmax <- max(sla_summary$step, na.rm = TRUE)

  # Compute SI at k_strict / k_tol if available
  SI_strict <- if (!is.na(k_strict)) compute_SI(sla_summary, k0 = k_strict, K = Kmax) else NA_real_
  SI_tol    <- if (!is.na(k_tol))    compute_SI(sla_summary, k0 = k_tol,    K = Kmax) else NA_real_

  # Extract last-step diagnostics
  last_row <- sla_summary[sla_summary$step == Kmax, , drop = FALSE]
  glm_last <- if (nrow(last_row) == 1) last_row$glm_pct else NA_real_
  ipt_last <- if (nrow(last_row) == 1) last_row$ipt_pct else NA_real_

  metrics <- data.frame(
    experiment = prefix,
    Kmax = Kmax,
    k_strict = k_strict,
    k_tol = k_tol,
    tol = tol,
    SI_strict = SI_strict,
    SI_tol = SI_tol,
    glm_pct_at_Kmax = glm_last,
    ipt_pct_at_Kmax = ipt_last,
    stringsAsFactors = FALSE
  )

  # Add meta fields if provided (seed, R, trim.level, etc.)
  if (length(meta) > 0) {
    for (nm in names(meta)) metrics[[nm]] <- meta[[nm]]
  }

  utils::write.csv(
    metrics,
    file = file.path(out_dir, paste0(prefix, "_metrics.csv")),
    row.names = FALSE
  )

  invisible(metrics)
}
