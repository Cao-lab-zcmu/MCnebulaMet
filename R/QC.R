# ==========================================================================
# QC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Title
#'
#' @param strings
#' @param patterns
#' @param target
#'
#' @return
#' @export get_metadata
get_metadata <- function(strings, patterns = NULL, target = "sample") {
  n <- length(strings)

  group_assign <- rep(NA_character_, n)
  group_names <- names(patterns)
  if (is.null(group_names)) {
    group_names <- paste0("Group", seq_along(patterns))
  }
 
  for (i in seq_along(patterns)) {
    idx <- which(is.na(group_assign) & grepl(patterns[i], strings, perl = TRUE))
    group_assign[idx] <- group_names[i]
  }
 
  df <- data.frame(
    target = strings,
    group = group_assign,
    stringsAsFactors = FALSE
  )
 
  df <- df[!is.na(df$group), ]
  colnames(df)[1] <- target
  tibble::as_tibble(df)
}

#' Title
#'
#' @param df
#' @param group
#'
#' @return
#' @export compute_rsd
compute_rsd <- function(intensity_df) {

  rsd <- apply(intensity_df, 1, cal_rsd)
  rsd[is.na(rsd)] <- Inf

  rsd_df <- data.frame(
    .features_id = rownames(intensity_df),
    RSD = rsd
  )

  rsd_sorted <- sort(rsd)
  cum_df <- data.frame(
   RSD = rsd_sorted,
   CumFraction = seq_along(rsd_sorted) / length(rsd_sorted)
  )

  qc_cols <- grep("^QC", colnames(intensity_df), value = TRUE)
  sample_cols <- setdiff(colnames(intensity_df), qc_cols)
  violin_df <- NULL
  
  if (length(qc_cols) > 0 && length(sample_cols) > 0) {
    rsd_qc <- apply(as.matrix(intensity_df[, qc_cols, drop = FALSE]), 1, cal_rsd)
    rsd_sample <- apply(as.matrix(intensity_df[, sample_cols, drop = FALSE]), 1, cal_rsd)

    violin_df <- rbind(
      data.frame(type = "QC", value = rsd_qc),
      data.frame(type = "sample", value = rsd_sample)
    )
    violin_df <- violin_df[is.finite(violin_df$value), ]

  }
  n_low <- sum(rsd < 30, na.rm = TRUE)
  percent_low <- round(n_low / length(rsd) * 100, 2)

  message(paste0(
    n_low, " features (", percent_low, "%) have RSD < 30%"
  ))

  res <- list(
    rsd_df = rsd_df,
    cum_df = cum_df
  )
  if (!is.null(violin_df)) {
    res$violin_df <- violin_df
  }

  return(res)
}

cal_rsd <- function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.na(m) || m == 0) return(Inf)
    sd(x, na.rm = TRUE) / m * 100
}

draw_rsd_scat <- function(rsd_df) {

  rsd_scat_p <- ggplot(rsd_df, aes(x = .features_id, y = RSD)) +
    geom_point(color = "steelblue") +
    geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
    #theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(title = "RSD Scatter Plot", y = "RSD (%)", x = "Features")

  return(rsd_scat_p)
}

draw_rsd_cumu <- function(cum_df) {

  rsd_cum_p <- ggplot(cum_df, aes(x = RSD, y = CumFraction)) +
    geom_line(
      color = "#4C78A8",
      linewidth = 0.9
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      linewidth = 0.6,
      color = "#E45756"
    ) +
    annotate(
      "text",
      x = 30,
      y = 0.05,
      label = "30%",
      hjust = -0.2,
      size = 3
    ) +
    coord_cartesian(xlim = c(0, max(cum_df$RSD, na.rm = TRUE)),
                    ylim = c(0, 1)) +
    labs(
      x = "RSD (%)",
      y = "Cumulative fraction"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
  return(rsd_cum_p)
}

draw_rsd_violin <- function(violin_df) {

  has_qc <- any(violin_df$type == "QC")
  if (!has_qc) {
    message("No QC samples found, plotting only Sample")
  }

  n_qc <- sum(violin_df$type == "QC")
  n_sample <- sum(violin_df$type == "Sample")

  label_list <- c(
    QC = paste0("QC (n = ", n_qc, ")"),
    Sample = paste0("Sample (n = ", n_sample, ")")
  )

  ggplot(violin_df, aes(x = type, y = value)) +

    # violin（更干净）
    geom_violin(
      aes(fill = type),
      trim = FALSE,
      alpha = 0.8,
      color = NA
    ) +

    # boxplot（细一点，更精致）
    geom_boxplot(
      width = 0.10,
      outlier.size = 0.5,
      fill = "white",
      color = "black",
      linewidth = 0.4
    ) +

    # cutoff（更克制）
    geom_hline(
      yintercept = 30,
      linetype = "dashed",
      linewidth = 0.6,
      color = "#E45756"
    ) +

    # cutoff标注（关键）
    annotate(
      "text",
      x = 1.5,
      y = 32,
      label = "30%",
      size = 3
    ) +

    # 配色（顶刊常用）
    scale_fill_manual(values = c(
      QC = "#4C78A8",
      Sample = "#E45756"
    )) +

    scale_x_discrete(labels = label_list) +

    coord_cartesian(ylim = c(0, max(violin_df$value, na.rm = TRUE))) +

    labs(
      x = NULL,
      y = "RSD (%)"
    ) +

    theme_classic(base_size = 12) +

    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),

      legend.position = "none",

      panel.grid = element_blank(),

      plot.margin = margin(8, 12, 8, 8)
    )
}

#draw_rsd_violin <- function(violin_df) {
#
#  has_qc <- any(violin_df$type == "QC")
#  if (!has_qc) {
#    message("No QC samples found, plotting only Sample")
#  }
#  n_qc <- sum(violin_df$type == "QC")
#  n_sample <- sum(violin_df$type == "Sample")
#
#  label_list <- c(
#    QC = paste0("QC (n = ", n_qc, ")"),
#    Sample = paste0("Sample (n = ", n_sample, ")")
#  )
#
#  rsd_violin_p <- ggplot(violin_df, aes(x = type, y = value, fill = type)) +
#    geom_violin(trim = FALSE, alpha = 0.8) +
#    geom_boxplot(width = 0.12, outlier.size = 0.6,
#                 fill = "white", color = "black") +
#    geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
#    labs(x = "", y = "RSD (%)") +
#    scale_x_discrete(labels = label_list) +
#    theme_bw() +
#    theme(legend.position = "none")
#
#  return(rsd_violin_p)
#}

#draw_intentity <- function(df_plot) {
#  p <-
#    ggplot(df_plot, aes(x = sample, y = log10(intensity), fill = group)) +
#    geom_boxplot() +
#    scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(df_plot$group)), "Set1")) +
#    labs(x = "", y = "Intensity") +
#    theme_bw()
#  return(p)
#}

draw_intentity <- function(df_plot) {

  ggplot(df_plot, aes(x = sample, y = intensity)) +

    # 箱线图（细一点）
    geom_boxplot(
      aes(fill = group),
      width = 0.6,
      outlier.size = 0.4,
      linewidth = 0.4
    ) +
    #geom_jitter(
    #  width = 0.2,
    #  size = 0.5,
    #  alpha = 0.3,
    #  color = "black"
    #) +
    scale_y_log10() +
    scale_fill_manual(values = c(
      "#4C78A8",
      "#E45756",
      "#54A24B"
    )) +

    labs(
      x = NULL,
      y = expression(Log[10]~"intensity"),
      fill = NULL
    ) +

    theme_classic(base_size = 12) +

    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11),

      legend.position = "top",
      legend.text = element_text(size = 10),

      panel.grid = element_blank(),

      plot.margin = margin(8, 12, 8, 8)
    )
}

draw_corr <- function(quant, sample_info) {
  corr <- stats::cor(quant, use = "pairwise.complete.obs")
  corr_df <- reshape2::melt(corr)
  sample_order <- sample_info$sample[order(sample_info$group)]
  corr_df$Var1 <- factor(corr_df$Var1, levels = sample_order)
  corr_df$Var2 <- factor(corr_df$Var2, levels = sample_order)
  ggplot2::ggplot(corr_df, ggplot2::aes(Var1, Var2, fill = value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756", midpoint = 0.8, limits = c(0, 1), name = "Correlation") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7), axis.text.y = ggplot2::element_text(size = 7), legend.position = "right", legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8), panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(8, 8, 8, 8))
}

draw_pca <- function(quant, sample_info, label = FALSE) {
  intensity_scaled_df <- scale(log2(quant + 1), center = TRUE, scale = TRUE)
  quant_t <- t(intensity_scaled_df)
  res_pca <- FactoMineR::PCA(quant_t, graph = FALSE)
  pca_df <- as.data.frame(res_pca$ind$coord)
  pca_df$group <- sample_info$group
  pca_df$sample <- rownames(pca_df)
  var_explained <- res_pca$eig[,2]
  pca_df$group <- factor(pca_df$group, levels = unique(sample_info$group))
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(Dim.1, Dim.2)) +
    ggplot2::stat_ellipse(ggplot2::aes(color = group), level = 0.7, linewidth = 0.6, alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(color = group), size = 2.6, alpha = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste0("PC1 (", round(var_explained[1],1), "%)"), y = paste0("PC2 (", round(var_explained[2],1), "%)"), color = NULL) +
    ggplot2::scale_color_manual(values = c("#4C78A8","#E45756","#54A24B")) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8,12,8,8)
    )
  if (label) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3, max.overlaps = Inf)
  }
  p
}
