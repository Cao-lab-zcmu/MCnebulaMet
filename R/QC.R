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
    geom_line(color = "darkgreen", size = 1) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Cumulative RSD Distribution", x = "RSD (%)", y = "Cumulative Fraction")

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

  rsd_violin_p <- ggplot(violin_df, aes(x = type, y = value, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.size = 0.6,
                 fill = "white", color = "black") +
    geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
    labs(x = "", y = "RSD (%)") +
    scale_x_discrete(labels = label_list) +
    theme_bw() +
    theme(legend.position = "none")

  return(rsd_violin_p)
}

draw_intentity <- function(df_plot) {
  p <-
    ggplot(df_plot, aes(x = sample, y = log10(intensity), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(df_plot$group)), "Set1")) +
    labs(x = "", y = "Intensity") +
    theme_bw()
  return(p)
}

