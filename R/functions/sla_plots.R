# ------------------------------------------------------------
# Plot helpers for SLA outputs
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ---- Combined plot: ATT band + DRDID CI (blue dashed) + GLM fallback (secondary axis) ----
save_sla_combined_plot <- function(summary_df,
                                   prefix,
                                   out_dir = "outputs/figures",
                                   width = 8,
                                   height = 5,
                                   dpi = 300) {
  
  ensure_dir(out_dir)
  
  # Scale factor to overlay glm_pct on ATT axis
  scale_factor <- stats::quantile(abs(summary_df$mean_ATT), probs = 0.98, na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
  
  # Check that DRDID CI columns exist (added by updated summarise_sla)
  has_drdid_ci <- all(c("drdid_ci_low", "drdid_ci_high") %in% names(summary_df))
  
  p_combined <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step)) +
    # Dispersion band across replications (stability)
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      alpha = 0.15
    ) +
    # Mean trajectory
    ggplot2::geom_line(ggplot2::aes(y = mean_ATT)) +
    # Optional: DRDID analytical CI aggregated across replications (blue dashed)
    { if (has_drdid_ci) ggplot2::geom_line(ggplot2::aes(y = drdid_ci_low),
                                           linetype = "dashed", colour = "blue") } +
    { if (has_drdid_ci) ggplot2::geom_line(ggplot2::aes(y = drdid_ci_high),
                                           linetype = "dashed", colour = "blue") } +
    # Zero reference
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    # GLM fallback overlay (scaled)
    ggplot2::geom_line(
      ggplot2::aes(y = glm_pct * scale_factor),
      colour = "red"
    ) +
    ggplot2::scale_y_continuous(
      name = "Mean ATT across replications",
      sec.axis = ggplot2::sec_axis(~ . / scale_factor,
                                   name = "GLM fallback rate (ps.flag = 2)")
    ) +
    ggplot2::labs(
      title = paste0(prefix, " — ATT stability, DRDID CI, and IPT convergence"),
      x = "Balanced subsample size k (k treated and k control units; k ≥ k_min)"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_att_glm_ci.png")),
    plot = p_combined,
    width = width, height = height, dpi = dpi
  )
  
  invisible(p_combined)
}

# ---- Main saver: creates all standard plots + the combined one ----
save_sla_plots <- function(summary_df,
                           si_grid,
                           prefix,
                           out_dir = "outputs/figures",
                           width = 8,
                           height = 5,
                           dpi = 300) {
  
  ensure_dir(out_dir)
  
  has_drdid_ci <- all(c("drdid_ci_low", "drdid_ci_high") %in% names(summary_df))
  
  # 1) ATT: mean +/- 1 SD + DRDID analytical CI (blue dashed)
  p_att <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      alpha = 0.15
    ) +
    ggplot2::geom_line(ggplot2::aes(y = mean_ATT)) +
    { if (has_drdid_ci) ggplot2::geom_line(ggplot2::aes(y = drdid_ci_low),
                                           linetype = "dashed", colour = "blue") } +
    { if (has_drdid_ci) ggplot2::geom_line(ggplot2::aes(y = drdid_ci_high),
                                           linetype = "dashed", colour = "blue") } +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = paste0(prefix, " — ATT trajectory (dispersion + DRDID 95% CI)"),
      x = "Balanced subsample size k (k treated and k control units; k ≥ k_min)",
      y = "ATT"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_att_band.png")),
    plot = p_att, width = width, height = height, dpi = dpi
  )
  
  # 2) GLM fallback rate
  p_glm <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step, y = glm_pct)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste0(prefix, " — GLM fallback (ps.flag = 2)"),
      x = "Balanced subsample size k (k treated and k control units; k ≥ k_min)",
      y = "GLM fallback rate"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_glm_fallback.png")),
    plot = p_glm, width = width, height = height, dpi = dpi
  )
  
  # 3) Stability Index
  p_si <- ggplot2::ggplot(si_grid, ggplot2::aes(x = k0, y = SI)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed") +
    ggplot2::labs(
      title = paste0(prefix, " — Stability Index SI(k0, K)"),
      x = "Starting ladder step k0",
      y = "SI"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_si_grid.png")),
    plot = p_si, width = width, height = height, dpi = dpi
  )
  
  # 4) Combined plot: ATT band + DRDID CI + GLM fallback
  p_combined <- save_sla_combined_plot(
    summary_df = summary_df,
    prefix     = prefix,
    out_dir    = out_dir,
    width      = width,
    height     = height,
    dpi        = dpi
  )
  
  invisible(list(p_att = p_att, p_glm = p_glm, p_si = p_si, p_combined = p_combined))
}
