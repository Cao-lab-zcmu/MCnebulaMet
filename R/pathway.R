plot_pathway_network <- function(
    cytos_plot_data,
    nodes_info,
    layout = c("kk", "stress"),
    label = c("name", "synonym"),
    seed = 1
) {
  
  layout <- match.arg(layout)
  label  <- match.arg(label)
  
  # ─────────────────────────────
  # Check data structure
  # ─────────────────────────────
  
  required_edge_cols <- c(
    "source",
    "target",
    "value",
    "type"
  )
  
  required_node_cols <- c(
    ".features_id",
    "synonym",
    "logFC",
    "adj.P.Val",
    "VIP"
  )
  
  missing_edge_cols <- setdiff(
    required_edge_cols,
    colnames(cytos_plot_data)
  )
  
  missing_node_cols <- setdiff(
    required_node_cols,
    colnames(nodes_info)
  )
  
  if (length(missing_edge_cols) > 0) {
    stop(
      "cytos_plot_data 缺少必要列: ",
      paste(missing_edge_cols, collapse = ", ")
    )
  }
  
  if (length(missing_node_cols) > 0) {
    stop(
      "nodes_info 缺少必要列: ",
      paste(missing_node_cols, collapse = ", ")
    )
  }
  
  
  # ─────────────────────────────
  # Check label column
  # ─────────────────────────────
  
  if (label == "name") {
    label_col <- "name"
  } else {
    label_col <- "synonym"
  }
  
  if (label_col == "synonym" &&
      !"synonym" %in% colnames(nodes_info)) {
    stop("nodes_info 中不存在 synonym 列。")
  }
  
  
  # ─────────────────────────────
  # Build graph
  # ─────────────────────────────
  
  g_tbl <- tidygraph::as_tbl_graph(
    cytos_plot_data,
    directed = FALSE
  )
  
  g_tbl <- g_tbl |>
    tidygraph::activate(nodes) |>
    dplyr::left_join(
      nodes_info,
      by = c("name" = ".features_id")
    ) |>
    dplyr::mutate(
      node_fill = factor(
        dplyr::case_when(
          logFC > 0.5  & adj.P.Val < 0.05 ~ "up",
          logFC < -0.5 & adj.P.Val < 0.05 ~ "down",
          VIP > 1      & adj.P.Val < 0.05 ~ "vip_high",
          TRUE                          ~ "no significance"
        )
      ),
      node_label = .data[[label_col]]
    )
  
  
  # ─────────────────────────────
  # Layout
  # ─────────────────────────────
  
  set.seed(seed)
  
  layout_data <- ggraph::create_layout(
    g_tbl,
    layout = layout
  )
  
  
  # ─────────────────────────────
  # Plot
  # ─────────────────────────────
  
  p <- ggraph::ggraph(layout_data) +
    
    # Edges
    ggraph::geom_edge_link(
      ggplot2::aes(
        color = type,
        width = value
      ),
      alpha = 0.45
    ) +
    
    ggraph::scale_edge_width(
      range = c(0.3, 1.5)
    ) +
    
    # Nodes
    ggraph::geom_node_point(
      ggplot2::aes(
        fill = node_fill,
        size = VIP
      ),
      shape = 21,
      color = "white",
      stroke = 0.7,
      alpha = 0.95
    ) +
    
    ggplot2::scale_size_continuous(
      range = c(3, 8),
      breaks = c(1, 1.5, 2, 3),
      name = "VIP"
    ) +
    
    # Labels
    ggraph::geom_node_text(
      ggplot2::aes(
        label = node_label
      ),
      repel = TRUE,
      size = 3.2,
      color = "grey20",
      family = "sans",
      fontface = "plain",
      max.overlaps = Inf
    ) +
    
    # Node colors
    ggplot2::scale_fill_manual(
      values = c(
        "up" = "#D73027",
        "down" = "#4575B4",
        "vip_high" = "#7B61A8",
        "no significance" = "#BDBDBD"
      ),
      name = "Metabolite status"
    ) +
    
    # Edge colors
    ggraph::scale_edge_color_manual(
      values = c(
        "spectral"  = "#8EC6D9",
        "structure" = "#496A88",
        "KEGG"      = "#E6A15C"
      ),
      name = "Association"
    ) +
    
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        15, 15, 15, 15
      )
    )
  
  return(p)
}
