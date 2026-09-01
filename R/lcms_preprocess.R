extractChromatograms <- function(
    raw_data,
    metadata,
    msLevel = 1L
) {
  
  stopifnot(
    nrow(metadata) == length(raw_data)
  )
  
  tic <- xcms::chromatogram(
    raw_data,
    aggregationFun = "sum",
    msLevel = msLevel
  )
  
  bpc <- xcms::chromatogram(
    raw_data,
    aggregationFun = "max",
    msLevel = msLevel
  )
  
  extract_data <- function(chrom, type) {
    
    data_list <- lapply(seq_along(chrom), function(i) {
      
      x <- chrom[[i]]
      
      if (length(x@rtime) == 0) {
        data.frame(
          rt = NA_real_,
          intensity = NA_real_,
          type = type,
          group = metadata$group[i],
          sample = metadata$sample[i]
        )
      } else {
        data.frame(
          rt = x@rtime / 60,
          intensity = x@intensity,
          type = type,
          group = metadata$group[i],
          sample = metadata$sample[i]
        )
      }
    })
    
    do.call(rbind, data_list)
  }
  
  list(
    tic = extract_data(tic, "TIC"),
    bpc = extract_data(bpc, "BPC")
  )
}

plotChromatograms <- function(
    chrom_data,
    rt_range = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Check input
  # ---------------------------------------------------------------------------
  
  if (!is.list(chrom_data) ||
      !all(c("tic", "bpc") %in% names(chrom_data))) {
    stop(
      "`chrom_data` must be a list containing `tic` and `bpc`."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Combine TIC and BPC
  # ---------------------------------------------------------------------------
  
  plot_data <- dplyr::bind_rows(
    chrom_data$tic,
    chrom_data$bpc
  )
  
  # ---------------------------------------------------------------------------
  # Retention time range
  # ---------------------------------------------------------------------------
  
  if (!is.null(rt_range)) {
    
    if (length(rt_range) != 2 ||
        !is.numeric(rt_range) ||
        any(!is.finite(rt_range)) ||
        rt_range[1] >= rt_range[2]) {
      stop(
        "`rt_range` must be a numeric vector of length 2, ",
        "e.g. `c(0, 10)`."
      )
    }
    
    plot_data <- plot_data |>
      dplyr::filter(
        rt >= rt_range[1],
        rt <= rt_range[2]
      )
    
    rt_label <- paste0(
      "RT: ",
      format(rt_range[1], trim = TRUE),
      "–",
      format(rt_range[2], trim = TRUE),
      " min"
    )
    
    # Automatically determine x-axis breaks
    rt_width <- diff(rt_range)
    
    if (rt_width <= 2) {
      x_breaks <- 0.25
    } else if (rt_width <= 5) {
      x_breaks <- 0.5
    } else if (rt_width <= 10) {
      x_breaks <- 1
    } else if (rt_width <= 20) {
      x_breaks <- 2
    } else if (rt_width <= 50) {
      x_breaks <- 5
    } else {
      x_breaks <- 10
    }
    
    x_scale <- ggplot2::scale_x_continuous(
      limits = rt_range,
      breaks = scales::breaks_width(x_breaks),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    )
    
  } else {
    
    # NULL = full retention-time range
    rt_label <- "RT: -Inf–Inf min"
    
    x_scale <- ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.005, 0.01))
    )
  }
  
  # ---------------------------------------------------------------------------
  # Y-axis
  # ---------------------------------------------------------------------------
  
  y_max <- max(
    plot_data$intensity,
    na.rm = TRUE
  )
  
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }
  
  y_scale <- ggplot2::scale_y_continuous(
    labels = scales::label_scientific(),
    limits = c(0, y_max * 1.08),
    expand = ggplot2::expansion(mult = c(0, 0))
  )
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = rt,
      y = intensity,
      group = interaction(sample, type),
      colour = type
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.55,
      alpha = 0.75
    ) +
    
    ggplot2::labs(
      title = NULL,
      x = "Retention time (min)",
      y = "Intensity",
      colour = NULL
    ) +
    
    x_scale +
    y_scale +
    
    # Top-journal style palette
    ggplot2::scale_colour_manual(
      values = c(
        TIC = "#3C5488",
        BPC = "#E64B35"
      )
    ) +
    
    ggplot2::theme_classic(
      base_size = 14
    ) +
    
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 16
      ),
      axis.text = ggplot2::element_text(
        size = 14,
        colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.5
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.5
      ),
      
      # TIC / BPC legend on the right
      legend.position = "right",
      legend.direction = "vertical",
      legend.text = ggplot2::element_text(
        size = 13
      ),
      legend.key.height = grid::unit(
        0.8,
        "cm"
      ),
      legend.key.width = grid::unit(
        0.8,
        "cm"
      )
    )
  
  # ---------------------------------------------------------------------------
  # Retention time annotation
  # ---------------------------------------------------------------------------
  
  p +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = rt_label,
      hjust = -0.05,
      vjust = 1.5,
      size = 4.5,
      fontface = "bold"
    )
}
