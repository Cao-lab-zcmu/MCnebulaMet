run_univariate <- function(mat, group, method = c("t.test", "wilcox.test"),
                          p_adjust_method = "BH") {
  
  method <- match.arg(method)
  
  g1 <- unique(group)[1]
  g2 <- unique(group)[2]
  
  idx1 <- which(group == g1)
  idx2 <- which(group == g2)
  
  res_list <- apply(mat, 1, function(x) {
    
    x1 <- x[idx1]
    x2 <- x[idx2]
    
    # fold change（log2 scale）
    fc <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)
    
    # 统计检验
    pval <- tryCatch({
      if (method == "t.test") {
        t.test(x1, x2)$p.value
      } else {
        wilcox.test(x1, x2)$p.value
      }
    }, error = function(e) NA)
    
    c(log2FC = fc, pvalue = pval)
  })
  
  res <- as.data.frame(t(res_list))
  
  res$padj <- p.adjust(res$pvalue, method = p_adjust_method)
  res$.features_id <- rownames(res)
  
  res <- dplyr::as_tibble(res) |>
    dplyr::relocate(.features_id)
  
  return(res)
}

#' Volcano plot
#'
#' @param df data.frame, differential analysis result
#' @param logfc_col column name for logFC
#' @param p_col column name for adjusted p value
#' @param fc_cutoff numeric, logFC threshold
#' @param p_cutoff numeric, adjusted p-value threshold
#' @param title plot title
#'
#' @return ggplot object
#' @export
plot_volcano <- function(df,
                         logfc_col = "logFC",
                         p_col = "adj.P.Val",
                         fc_cutoff = 1.5,
                         p_cutoff = 0.05,
                         title = "Volcano Plot",
                         subtitle = NULL) {
  df_plot <- df |>
    dplyr::mutate(
      .features_id = as.character(.features_id),
      neg_log10_q = -log10(.data[[p_col]]),
      sig = dplyr::case_when(
        .data[[p_col]] < p_cutoff & .data[[logfc_col]] >= fc_cutoff ~ "Up",
        .data[[p_col]] < p_cutoff & .data[[logfc_col]] <= -fc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data[[logfc_col]],
      y = neg_log10_q
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = sig),
      size = 1.8,
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(
      values = c(
        Down = "#4C6A92",
        NS   = "#CFCFCF",
        Up   = "#C77C6B"
      )
    ) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "#8F8F8F"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "#8F8F8F"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~adjusted~italic(P))
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 11.5,
        color = "grey30",
        margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(
        size = 13,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 11,
        face = "bold"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.7
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.4
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        linewidth = 0.8
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
  return(p)
}

#' Heatmap plot
#'
#' @param mat numeric matrix (features × samples)
#' @param scale_rows logical, whether to z-score by row
#' @param cluster_rows logical
#' @param cluster_cols logical
#' @param show_rownames logical
#' @param show_colnames logical
#'
#' @return pheatmap object
#' @export
plot_heatmap <- function(mat,
                         scale_rows = TRUE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = FALSE,
                         show_colnames = FALSE) {

  # 行标准化（热图常规操作）
  if (scale_rows) {
    mat <- t(scale(t(mat)))
  }

  # 颜色：和你 volcano 的低饱和风格一致
  heat_colors <- colorRampPalette(c(
    "#4C6A92",  # down (蓝)
    "#F5F5F5",  # 中性灰
    "#C77C6B"   # up (红)
  ))(100)

  pheatmap::pheatmap(
    mat,
    color = heat_colors,
    border_color = NA,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 9,
    treeheight_row = 40,
    treeheight_col = 40,
    angle_col = 45
  )
}
