# ------------------------------------------------------------
# Plot helpers for SLA outputs
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

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
      x = "Subsample size (treated = controls = step)",
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
      x = "Subsample size (step)",
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
    ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed") +
    ggplot2::labs(
      title = paste0(prefix, " — Stability Index SI(k0, K)"),
      x = "Initial k0",
      y = "SI"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_si_grid.png")),
    plot = p_si, width = width, height = height, dpi = dpi
  )

  invisible(list(p_att = p_att, p_glm = p_glm, p_si = p_si))
}
