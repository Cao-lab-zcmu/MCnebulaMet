get_mirror_data <- function(annotation, sample, library) {
  # Check
  if (nrow(annotation) != 1) {
    stop("annotation must contain exactly one row.")
  }

  # Prepare
  feature_id <- sub("_.*", "", annotation$ID)
  library_id <- annotation$`Library ID`

  spec_sample <- sample$peaks |>
    dplyr::filter(.features_id == feature_id)

  spec_library <- library$peaks |>
    dplyr::filter(id == library_id)

  list(
    sample = spec_sample,
    library = spec_library,
    meta = annotation
  )
}

plot_mirror <- function(
    spec_sample,
    spec_library,
    title = NULL,
    subtitle = NULL
) {
  xmax <- max(c(spec_sample$mz, spec_library$mz), na.rm = TRUE)

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = spec_sample,
      ggplot2::aes(
        x = mz,
        xend = mz,
        y = 0,
        yend = rel.int.
      ),
      linewidth = 0.9,
      colour = "#1B4F72"
    ) +
    ggplot2::geom_segment(
      data = spec_library,
      ggplot2::aes(
        x = mz,
        xend = mz,
        y = 0,
        yend = -rel_intensity
      ),
      linewidth = 0.9,
      colour = "#7B241C"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.8,
      colour = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 105,
      hjust = 0,
      size = 6,
      fontface = "bold",
      label = title
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = -105,
      hjust = 0,
      size = 6,
      fontface = "bold",
      label = subtitle
    ) +
    ggplot2::annotate(
      "text",
      x = xmax,
      y = 105,
      hjust = 1,
      size = 6,
      fontface = "bold",
      colour = "#2C3E50",
      label = "Sample"
    ) +
    ggplot2::annotate(
      "text",
      x = xmax,
      y = -105,
      hjust = 1,
      size = 6,
      fontface = "bold",
      colour = "#8B1E3F",
      label = "Library"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-110, 110)
    ) +
    ggplot2::labs(
      x = "m/z",
      y = "Relative Intensity (%)"
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text = ggplot2::element_text(size = 16),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}
plot_intensity_std <- function(i, stat, origin, sample_info) {
  feature_id <- stat$.features_id[i]
  dat <- origin |>
    dplyr::filter(.features_id == feature_id) |>
    tidyr::pivot_longer(
      cols = -.features_id,
      names_to = "sample",
      values_to = "intensity"
    ) |>
    dplyr::left_join(
      sample_info,
      by = "sample"
    ) |>
    dplyr::filter(!is.na(group))
  if (nrow(dat) == 0) return(NULL)
  dat$group <- factor(dat$group, levels = unique(dat$group))
  meta <- stat[i, ]
  group_levels <- levels(dat$group)
  if (length(group_levels) == 2) {
    group_colors <- c(
      "#1B4F72",
      "#7B241C"
    )
    names(group_colors) <- group_levels
  } else {
    group_colors <- c(
      "#1B4F72",
      "#7B241C",
      "#2E7D32",
      "#7B1FA2",
      "#EF6C00",
      "#00838F"
    )
    group_colors <- group_colors[
      seq_along(group_levels)
    ]
    names(group_colors) <- group_levels
  }
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = group,
      y = intensity,
      fill = group
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.55,
      linewidth = 0.9,
      outlier.shape = NA,
      alpha = 0.55,
      colour = "black"
    ) +
    ggplot2::geom_jitter(
      width = 0.10,
      height = 0,
      shape = 21,
      size = 3.5,
      stroke = 0.7,
      colour = "black",
      alpha = 0.95
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_scientific(
        digits = 2
      ),
      expand = ggplot2::expansion(
        mult = c(0.05, 0.15)
      )
    ) +
    ggplot2::labs(
      title = paste0(
        "Feature ",
        feature_id
      ),
      subtitle = paste0(
        "logFC = ",
        sprintf("%.2f", meta$logFC),
        "   |   FDR = ",
        format(
          meta$adj.P.Val,
          scientific = TRUE,
          digits = 2
        ),
        "   |   VIP = ",
        sprintf("%.2f", meta$VIP)
      ),
      x = NULL,
      y = "Intensity"
    ) +
    ggplot2::theme_classic(
      base_size = 18
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(
        size = 16,
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.x = ggplot2::element_text(
        size = 16,
        face = "bold",
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 14,
        colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.8,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.8,
        colour = "black"
      ),
      plot.title = ggplot2::element_text(
        size = 18,
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 12.5,
        colour = "#4D4D4D",
        margin = ggplot2::margin(b = 10)
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        10, 15, 10, 10
      )
    )
  return(p)
}
