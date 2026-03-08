# ------------------------------------------------------------
# Plot helpers for SLA outputs
# Double encoding: color + linetype (works in color AND B&W print)
#
# Element          Color     Linetype
# -------          -----     --------
# Mean ATT         black     solid
# DRDID 95% CI     blue      dotdash
# GLM fallback     red       longdash
# Zero reference   grey60    solid, thin     (not in legend)
# 0.95 threshold   black     dotted          (not in legend)
# Ribbon           grey50    filled          (not in legend)
#
# x-axis (ATT/GLM plots): Ladder step k (k treated, k control)
# x-axis (SI plot):        k_0 (stability window start)
# ------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Shared scale definitions — identical name/breaks forces legend merge
.line_labels  <- c("Mean ATT", "DRDID 95% CI", "GLM fallback rate")
.line_colours <- c("Mean ATT" = "black", "DRDID 95% CI" = "blue",
                   "GLM fallback rate" = "red")
.line_types   <- c("Mean ATT" = "solid", "DRDID 95% CI" = "dotdash",
                   "GLM fallback rate" = "longdash")

.scale_colour <- ggplot2::scale_colour_manual(
  name   = "",
  values = .line_colours,
  breaks = .line_labels
)
.scale_ltype <- ggplot2::scale_linetype_manual(
  name   = "",
  values = .line_types,
  breaks = .line_labels
)

# Two-element version (no GLM) for att_band
.line_labels2  <- c("Mean ATT", "DRDID 95% CI")
.scale_colour2 <- ggplot2::scale_colour_manual(
  name   = "",
  values = .line_colours[.line_labels2],
  breaks = .line_labels2
)
.scale_ltype2 <- ggplot2::scale_linetype_manual(
  name   = "",
  values = .line_types[.line_labels2],
  breaks = .line_labels2
)

# Shared x-axis labels
.xlab_k <- expression("Ladder step" ~ k ~ "(" * k ~ "treated," ~ k ~ "control)")
.xlab_k0 <- expression(k[0] ~ "(ladder starting point)")

# Shared legend theme — compact box anchored bottom-right inside the panel
.legend_theme <- ggplot2::theme(
  legend.position       = c(0.98, 0.02),
  legend.justification  = c("right", "bottom"),
  legend.direction      = "vertical",
  legend.box.background = ggplot2::element_rect(colour = "grey70", fill = "white",
                                                linewidth = 0.4),
  legend.margin         = ggplot2::margin(4, 6, 4, 6),
  legend.key.width      = ggplot2::unit(1.6, "cm"),
  legend.key.height     = ggplot2::unit(0.45, "cm"),
  legend.text           = ggplot2::element_text(size = 8),
  legend.title          = ggplot2::element_blank()
)

# ---- Combined plot: ATT band + DRDID CI + GLM fallback ----
save_sla_combined_plot <- function(summary_df,
                                   prefix,
                                   out_dir = "outputs/figures",
                                   width   = 8,
                                   height  = 5,
                                   dpi     = 300) {
  
  ensure_dir(out_dir)
  
  scale_factor <- stats::quantile(abs(summary_df$mean_ATT), probs = 0.98, na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
  
  has_drdid_ci <- all(c("drdid_ci_low", "drdid_ci_high") %in% names(summary_df))
  
  p_combined <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step)) +
    # Dispersion ribbon — grey filled, not in legend
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      fill = "grey50", alpha = 0.15, colour = NA
    ) +
    # Zero reference — thin grey solid, not in legend
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "grey60", linewidth = 0.35) +
    # Mean ATT — black solid
    ggplot2::geom_line(
      ggplot2::aes(y = mean_ATT,
                   colour   = "Mean ATT",
                   linetype = "Mean ATT"),
      linewidth = 0.7
    ) +
    # DRDID CI — blue dotdash
    { if (has_drdid_ci) ggplot2::geom_line(
      ggplot2::aes(y = drdid_ci_low,
                   colour   = "DRDID 95% CI",
                   linetype = "DRDID 95% CI"),
      linewidth = 0.6) } +
    { if (has_drdid_ci) ggplot2::geom_line(
      ggplot2::aes(y = drdid_ci_high,
                   colour   = "DRDID 95% CI",
                   linetype = "DRDID 95% CI"),
      linewidth = 0.6) } +
    # GLM fallback — red longdash
    ggplot2::geom_line(
      ggplot2::aes(y = glm_pct * scale_factor,
                   colour   = "GLM fallback rate",
                   linetype = "GLM fallback rate"),
      linewidth = 0.6
    ) +
    .scale_colour +
    .scale_ltype  +
    ggplot2::scale_y_continuous(
      name = expression(widehat(ATT) ~ "across replications"),
      sec.axis = ggplot2::sec_axis(
        ~ . / scale_factor,
        name = expression("GLM fallback rate" ~ (widehat(ps)[k]^{(r)} == 2))
      )
    ) +
    ggplot2::labs(x = .xlab_k) +
    ggplot2::theme_minimal() +
    .legend_theme
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_att_glm_ci.jpeg")),
    plot = p_combined,
    width = width, height = height, dpi = dpi
  )
  
  invisible(p_combined)
}

# ---- Main saver ----
save_sla_plots <- function(summary_df,
                           si_grid,
                           prefix,
                           out_dir = "outputs/figures",
                           width   = 8,
                           height  = 5,
                           dpi     = 300) {
  
  ensure_dir(out_dir)
  
  has_drdid_ci <- all(c("drdid_ci_low", "drdid_ci_high") %in% names(summary_df))
  
  # 1) ATT band + DRDID CI — legend bottom-right inside panel
  p_att <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_ATT - sd_ATT, ymax = mean_ATT + sd_ATT),
      fill = "grey50", alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "grey60", linewidth = 0.35) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean_ATT,
                   colour   = "Mean ATT",
                   linetype = "Mean ATT"),
      linewidth = 0.7
    ) +
    { if (has_drdid_ci) ggplot2::geom_line(
      ggplot2::aes(y = drdid_ci_low,
                   colour   = "DRDID 95% CI",
                   linetype = "DRDID 95% CI"),
      linewidth = 0.6) } +
    { if (has_drdid_ci) ggplot2::geom_line(
      ggplot2::aes(y = drdid_ci_high,
                   colour   = "DRDID 95% CI",
                   linetype = "DRDID 95% CI"),
      linewidth = 0.6) } +
    .scale_colour2 +
    .scale_ltype2  +
    ggplot2::labs(
      x = .xlab_k,
      y = expression(widehat(ATT))
    ) +
    ggplot2::theme_minimal() +
    .legend_theme
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_att_band.jpeg")),
    plot = p_att, width = width, height = height, dpi = dpi
  )
  
  # 2) GLM fallback rate — no legend (single line)
  p_glm <- ggplot2::ggplot(summary_df, ggplot2::aes(x = step, y = glm_pct)) +
    ggplot2::geom_line(linetype = "longdash", colour = "red", linewidth = 0.6) +
    ggplot2::labs(
      x = .xlab_k,
      y = expression("GLM fallback rate" ~ (widehat(ps)[k]^{(r)} == 2))
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_glm_fallback.jpeg")),
    plot = p_glm, width = width, height = height, dpi = dpi
  )
  
  # 3) Stability Index — no legend (single line)
  p_si <- ggplot2::ggplot(si_grid, ggplot2::aes(x = k0, y = SI)) +
    ggplot2::geom_line(colour = "black", linewidth = 0.7) +
    # 0.95 threshold — black dotted
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dotted",
                        colour = "black", linewidth = 0.5) +
    ggplot2::labs(
      x = .xlab_k0,
      y = expression(SI(k[0], K))
    ) +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_si_grid.jpeg")),
    plot = p_si, width = width, height = height, dpi = dpi
  )
  
  # 4) Combined plot
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
