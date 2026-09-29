plot_id_venn <- function(
    df1,
    df2,
    name1 = "MS/MS Annotation",
    name2 = "GNPS Library"
) {
  # Prepare
  ids1 <- unique(stats::na.omit(df1$ID))
  ids2 <- unique(stats::na.omit(df2$ID))

  n1 <- length(ids1)
  n2 <- length(ids2)
  n12 <- length(intersect(ids1, ids2))
  n1_only <- length(setdiff(ids1, ids2))
  n2_only <- length(setdiff(ids2, ids1))

  # Plot
  circle1 <- data.frame(
    x = 0.42,
    y = 0.5,
    r = 0.30
  )

  circle2 <- data.frame(
    x = 0.58,
    y = 0.5,
    r = 0.30
  )

  theta <- seq(0, 2 * pi, length.out = 500)

  circle1 <- data.frame(
    x = 0.42 + 0.30 * cos(theta),
    y = 0.50 + 0.30 * sin(theta)
  )

  circle2 <- data.frame(
    x = 0.58 + 0.30 * cos(theta),
    y = 0.50 + 0.30 * sin(theta)
  )

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = circle1,
      ggplot2::aes(x, y),
      fill = "#D9E2F3",
      colour = "#465C88",
      alpha = 0.65,
      linewidth = 0.8
    ) +
    ggplot2::geom_polygon(
      data = circle2,
      ggplot2::aes(x, y),
      fill = "#E8D8D5",
      colour = "#C65D3A",
      alpha = 0.65,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = 0.30,
      y = 0.84,
      label = name1,
      family = "Times New Roman",
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::annotate(
      "text",
      x = 0.70,
      y = 0.84,
      label = name2,
      family = "Times New Roman",
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::annotate(
      "text",
      x = 0.30,
      y = 0.50,
      label = n1_only,
      family = "Times New Roman",
      fontface = "bold",
      size = 7
    ) +
    ggplot2::annotate(
      "text",
      x = 0.50,
      y = 0.50,
      label = n12,
      family = "Times New Roman",
      fontface = "bold",
      size = 7
    ) +
    ggplot2::annotate(
      "text",
      x = 0.70,
      y = 0.50,
      label = n2_only,
      family = "Times New Roman",
      fontface = "bold",
      size = 7
    ) +
    ggplot2::coord_fixed(
      xlim = c(0.05, 0.95),
      ylim = c(0.08, 0.92)
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}
