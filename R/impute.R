#' Title
#'
#' @param df
#'
#' @return
#' @export
remove_noisy <- function(df, qc_mv_cutoff = 0.2, sample_mv_cutoff = 0.5) {
  qc_col <- grep("QC", colnames(df), value = TRUE)
  sample_col <- setdiff(colnames(df), qc_col)

  if (length(qc_col) == 0) {
    stop("❌ No QC columns found (column name must contain 'QC').")
  }


  df_qc     <- df[, qc_col, drop = FALSE]
  df_sample <- df[, sample_col, drop = FALSE]

  qc_mv_ratio <- rowSums(is.na(df_qc)) / ncol(df_qc)
  sample_mv_ratio <- rowSums(is.na(df_sample)) / ncol(df_sample)

  high_mv_idx <- which(qc_mv_ratio > qc_mv_cutoff | sample_mv_ratio > sample_mv_cutoff)
  removed_pct <- length(high_mv_idx) / nrow(df) * 100
  message(sprintf("🧹 Removed %d features (%.2f%%) due to high missingness. (QC > %.2f, Sample > %.2f)",
                  length(high_mv_idx), removed_pct, qc_mv_cutoff, sample_mv_cutoff))

  if (length(high_mv_idx) > 0) {
    message("🚩 Removed feature indices:")
    message(paste(high_mv_idx, collapse = ", "))
  }

  df_filtered <- df[-high_mv_idx, , drop = FALSE]
  return(df_filtered)
}


#' Title
#'
#' @param df
#' @param method
#'
#' @return
#' @export impute_mv
impute_mv <- function(df, method = c("knn", "median", "min")) {

  df <- dplyr::select(df, dplyr::where(is.double))

  if (sum(is.na(df)) == 0) {
    message("✅ No missing values detected. Skipping imputation.")
    return(df)
  }

  if (method == "median") {
    df_imputed <- apply(df, 2, function(x) {
      ifelse(is.na(x), median(x, na.rm = TRUE), x)
    })
  } else if (method == "min") {
    df_imputed <- apply(df, 2, function(x) {
      ifelse(is.na(x), min(x, na.rm = TRUE), x)
    })
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE)) {
      BiocManager::install("impute")
    }
    df_imputed <- impute::impute.knn(as.matrix(df))$data
  }

  return(df_imputed)
}

#' Title
#'
#' @param df
#' @param metadata
#'
#' @return
#' @export detect_mv
detect_mv <- function(df, metadata) {

  numeric_df <- dplyr::select(df, dplyr::where(is.numeric))

  df_mv <- data.frame(
    sample = names(numeric_df),
    mv_rate = colSums(is.na(numeric_df)) / nrow(numeric_df) * 100
  )
  df_mv <- dplyr::left_join(df_mv, metadata, by = "sample")

  p <- draw_mv(df_mv)
  ggsave("missing_values_rate.png", p, width = 10, height = 10, dpi = 300)

  if (sum(is.na(df)) == 0) {
    message("✅ No missing values detected. Skipping imputation.")
  } else {
    message(paste0("⚠ ", sum(is.na(df)), " missing values detected."))
  }

}

draw_mv <- function(df_mv) {

  mv_p <- ggplot(df_mv, aes(x = sample, y = mv_rate, fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(df_mv$group)), "Set1")) +
    labs(x = "", y = "Number of missing values", title = "Missing Values Rate") +
    ylim(0, 100) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(mv_p)
}
