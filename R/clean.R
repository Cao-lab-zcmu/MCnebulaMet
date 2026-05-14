#' Remove noisy features based on missing value and RSD
#'
#' @param df data.frame / matrix (features × samples)
#' @param qc_mv_cutoff numeric, QC missing value ratio cutoff (e.g. 0.2)
#' @param sample_mv_cutoff numeric, sample missing value ratio cutoff (e.g. 0.5)
#' @param rsd_cutoff numeric, QC RSD cutoff (e.g. 0.3)
#'
#' @return filtered data.frame
#' @export
remove_noisy_feature <- function(df,
                         qc_mv_cutoff = NULL,
                         sample_mv_cutoff = NULL,
                         rsd_cutoff = NULL) {

  qc_col <- grep("QC", colnames(df), value = TRUE)
  sample_col <- setdiff(colnames(df), qc_col)

  has_qc <- length(qc_col) > 0

  if (!has_qc && (!is.null(qc_mv_cutoff) || !is.null(rsd_cutoff))) {
    message("⚠️ No QC columns found → QC-based filtering will be skipped.")
  }

  if (!is.null(qc_mv_cutoff) && has_qc) {
    df_qc <- df[, qc_col, drop = FALSE]
    qc_mv_ratio <- rowSums(is.na(df_qc)) / ncol(df_qc)
  } else {
    qc_mv_ratio <- rep(0, nrow(df))
  }

  if (!is.null(sample_mv_cutoff)) {
    df_sample <- df[, sample_col, drop = FALSE]
    sample_mv_ratio <- rowSums(is.na(df_sample)) / ncol(df_sample)
  } else {
    sample_mv_ratio <- rep(0, nrow(df))
  }

  if (!is.null(rsd_cutoff) && has_qc) {
    df_qc <- df[, qc_col, drop = FALSE]

    rsd <- apply(df_qc, 1, function(x) {
      if (all(is.na(x))) return(NA)
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    })
  } else {
    rsd <- rep(0, nrow(df))
  }

  remove_idx <- rep(FALSE, nrow(df))

  if (!is.null(qc_mv_cutoff) && has_qc) {
    remove_idx <- remove_idx | (qc_mv_ratio > qc_mv_cutoff)
  }

  if (!is.null(sample_mv_cutoff)) {
    remove_idx <- remove_idx | (sample_mv_ratio > sample_mv_cutoff)
  }

  if (!is.null(rsd_cutoff) && has_qc) {
    remove_idx <- remove_idx | (rsd > rsd_cutoff)
  }

  high_mv_idx <- which(remove_idx)

  removed_pct <- length(high_mv_idx) / nrow(df) * 100

  message(sprintf("🧹 Removed %d features (%.2f%%)",
                  length(high_mv_idx), removed_pct))

  if (!is.null(qc_mv_cutoff) && has_qc) {
    message(sprintf("   QC missing > %.2f", qc_mv_cutoff))
  }

  if (!is.null(sample_mv_cutoff)) {
    message(sprintf("   Sample missing > %.2f", sample_mv_cutoff))
  }

  if (!is.null(rsd_cutoff) && has_qc) {
    message(sprintf("   RSD > %.2f", rsd_cutoff))
  }

  if (length(high_mv_idx) > 0) {
    message("🚩 Removed feature indices:")
    message(paste(high_mv_idx, collapse = ", "))
  }

  df_filtered <- df[!remove_idx, , drop = FALSE]

  return(df_filtered)
}

#' Title
#'
#' @param quant
#' @param method
#'
#' @return
#' @export impute_mv
impute_mv <- function(quant, method = c("knn", "median", "min")) {

  .check_feature_table(quant)
  method <- match.arg(method)

  if (sum(is.na(quant)) == 0) {
    message("✅ No missing values detected. Skipping imputation.")
    return(quant)
  }

  if (method == "median") {
    quant_imputed <- apply(quant, 2, function(x) {
      ifelse(is.na(x), median(x, na.rm = TRUE), x)
    })
  } else if (method == "min") {
    quant_imputed <- apply(quant, 2, function(x) {
      ifelse(is.na(x), min(x, na.rm = TRUE), x)
    })
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE)) {
      BiocManager::install("impute")
    }
    quant_imputed <- impute::impute.knn(as.matrix(quant))$data
  }

  return(quant_imputed)
}

#' Title
#'
#' @param df
#' @param metadata
#'
#' @return
#' @export detect_mv
detect_mv <- function(df, metadata) {

  .check_feature_table(df)

  df_mv <- data.frame(
    sample = names(df),
    mv_rate = colSums(is.na(df)) / nrow(df) * 100
  )
  if (sum(is.na(df_mv)) == 0) {
    message("✅ No missing values detected. Skipping imputation.")
  } else {
    message(paste0("⚠ ", sum(is.na(df_mv)), " missing values detected."))
  }
  df_mv <- dplyr::left_join(df_mv, metadata, by = "sample")

  return(df_mv)
}

draw_mv <- function(df_mv) {
  mv_p <- ggplot(df_mv, aes(x = sample, y = mv_rate)) +
    geom_col(
      aes(fill = group),
      width = 0.7,
      alpha = 0.9
    ) +
    scale_fill_manual(values = c(
      "#4C78A8",  # blue
      "#E45756",  # red
      "#54A24B"   # green
    )) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      x = NULL,
      y = "Missing values (%)",
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

  return(mv_p)
}
