# ------------------------------------------------------------
# Plot helpers for SLA outputs
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ---- Combined plot: ATT band + GLM fallback (secondary axis) ----
save_sla_combined_plot <- function(summary_df,
                                   prefix,
                                   out_dir = "outputs/figures",
                                   width = 8,
                                   height = 5,
                                   dpi = 300) {
  
  ensure_dir(out_dir)
  
  # Scale factor to overlay glm_pct on ATT axis
  scale_factor <- max(abs(summary_df$mean_ATT), na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
  
  p_combined <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      alpha = 0.15
    ) +
    ggplot2::geom_line(ggplot2::aes(y = mean_ATT)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
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
      title = paste0(prefix, " — ATT stability and IPT convergence"),
      x = "Balanced subsample size k (k treated and k control units; k≥10)"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_att_and_glm.png")),
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
  
  # 1) ATT mean +/- 1 SD
  p_att <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step, y = mean_ATT)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      alpha = 0.2
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = paste0(prefix, " — ATT stability"),
      x = "Balanced subsample size k (k treated and k control units; k≥10)",
      y = "Mean ATT across replications"
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
      x = "Balanced subsample size k (k treated and k control units; k≥10)",
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
      x = "Ladder step",
      y = "SI"
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_si_grid.png")),
    plot = p_si, width = width, height = height, dpi = dpi
  )
  
  # 4) Combined plot: ATT band + GLM fallback
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
