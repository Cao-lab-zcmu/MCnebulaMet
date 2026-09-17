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

plot_kegg_pathway <- function(
    pathway_id,
    nodes_info,
    layout = c("kk", "stress"),
    label = c("name", "synonym"),
    seed = 1,
    species = "hsa"
) {
  
  layout <- match.arg(layout)
  label  <- match.arg(label)
  
  # ─────────────────────────────
  # Check packages
  # ─────────────────────────────
  
  required_pkgs <- c(
    "xml2",
    "dplyr",
    "tidyr",
    "tibble",
    "tidygraph",
    "ggraph",
    "ggplot2"
  )
  
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "缺少 R 包: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }
  
  
  # ─────────────────────────────
  # Check nodes_info
  # ─────────────────────────────
  
  required_node_cols <- c(
    ".features_id",
    "synonym",
    "logFC",
    "adj.P.Val",
    "VIP"
  )
  
  missing_node_cols <- setdiff(
    required_node_cols,
    colnames(nodes_info)
  )
  
  if (length(missing_node_cols) > 0) {
    stop(
      "nodes_info 缺少必要列: ",
      paste(missing_node_cols, collapse = ", ")
    )
  }
  
  
  # ─────────────────────────────
  # KEGG pathway ID
  # ─────────────────────────────
  
  if (!grepl("^\\d{5}$", pathway_id)) {
    pathway_id <- sub(
      "^path:",
      "",
      pathway_id
    )
    
    if (!grepl("^\\d{5}$", pathway_id)) {
      stop(
        "pathway_id 应该类似 '00010' 或 'path:hsa00010'"
      )
    }
  }
  
  
  pathway_name <- paste0(
    species,
    pathway_id
  )
  
  
  # ─────────────────────────────
  # Download KGML
  # ─────────────────────────────
  
  kgml_url <- paste0(
    "https://rest.kegg.jp/get/",
    pathway_name,
    "/kgml"
  )
  
  kgml <- xml2::read_xml(kgml_url)
  
  
  # ─────────────────────────────
  # Extract entries
  # ─────────────────────────────
  
  entries <- xml2::xml_find_all(
    kgml,
    ".//entry"
  )
  
  entry_tbl <- tibble::tibble(
    id = as.integer(
      xml2::xml_attr(entries, "id")
    ),
    type = xml2::xml_attr(
      entries,
      "type"
    ),
    name = xml2::xml_attr(
      entries,
      "name"
    ),
    entry_name = xml2::xml_attr(
      entries,
      "name"
    )
  )
  
  
  # ─────────────────────────────
  # Keep compounds
  # ─────────────────────────────
  
  compound_tbl <- entry_tbl |>
    dplyr::filter(
      type == "compound"
    ) |>
    dplyr::mutate(
      kegg = sub(
        "^cpd:",
        "",
        name
      )
    )
  
  
  # ─────────────────────────────
  # Extract relations
  # ─────────────────────────────
  
  relations <- xml2::xml_find_all(
    kgml,
    ".//relation"
  )
  
  relation_tbl <- tibble::tibble(
    source_id = as.integer(
      xml2::xml_attr(
        relations,
        "entry1"
      )
    ),
    target_id = as.integer(
      xml2::xml_attr(
        relations,
        "entry2"
      )
    ),
    type = xml2::xml_attr(
      relations,
      "type"
    )
  )
  
  
  # ─────────────────────────────
  # Compound → Compound edges
  # ─────────────────────────────
  
  edge_tbl <- relation_tbl |>
    dplyr::inner_join(
      compound_tbl |>
        dplyr::select(
          source_id = id,
          source = kegg
        ),
      by = "source_id"
    ) |>
    dplyr::inner_join(
      compound_tbl |>
        dplyr::select(
          target_id = id,
          target = kegg
        ),
      by = "target_id"
    ) |>
    dplyr::select(
      source,
      target,
      type
    ) |>
    dplyr::distinct()
  
  
  # ─────────────────────────────
  # Build nodes
  # ─────────────────────────────
  
  node_tbl <- compound_tbl |>
    dplyr::select(
      name = kegg
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(
      nodes_info,
      by = c(
        "name" = ".features_id"
      )
    ) |>
    dplyr::mutate(
      logFC = dplyr::coalesce(
        logFC,
        0
      ),
      adj.P.Val = dplyr::coalesce(
        adj.P.Val,
        1
      ),
      VIP = dplyr::coalesce(
        VIP,
        0
      ),
      node_fill = factor(
        dplyr::case_when(
          logFC > 0.5 &
            adj.P.Val < 0.05 ~ "up",
          
          logFC < -0.5 &
            adj.P.Val < 0.05 ~ "down",
          
          VIP > 1 &
            adj.P.Val < 0.05 ~ "vip_high",
          
          TRUE ~ "no significance"
        ),
        levels = c(
          "up",
          "down",
          "vip_high",
          "no significance"
        )
      )
    )
  
  
  # ─────────────────────────────
  # Label
  # ─────────────────────────────
  
  if (label == "name") {
    
    node_tbl <- node_tbl |>
      dplyr::mutate(
        node_label = name
      )
    
  } else {
    
    node_tbl <- node_tbl |>
      dplyr::mutate(
        node_label = dplyr::if_else(
          is.na(synonym) |
            synonym == "",
          name,
          synonym
        )
      )
  }
  
  
  # ─────────────────────────────
  # Build graph
  # ─────────────────────────────
  
  g_tbl <- tidygraph::tbl_graph(
    nodes = node_tbl,
    edges = edge_tbl,
    directed = FALSE
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
  
  p <- ggraph::ggraph(
    layout_data
  ) +
    
    # Edges
    ggraph::geom_edge_link(
      ggplot2::aes(
        color = type
      ),
      alpha = 0.45,
      width = 0.7
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
      breaks = c(
        1,
        1.5,
        2,
        3
      ),
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
    ggraph::scale_edge_color_discrete(
      name = "KEGG relation"
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
        15,
        15,
        15,
        15
      )
    )
  
  
  return(
    list(
      plot = p,
      nodes = node_tbl,
      edges = edge_tbl,
      kgml = kgml
    )
  )
}

plot_kegg_pathway2 <- function(
    pathway_id,
    nodes_info,
    layout = c("kk", "stress"),
    label = c("name", "synonym"),
    seed = 1,
    species = "hsa"
) {

  layout <- match.arg(layout)
  label <- match.arg(label)

  # ─────────────────────────────
  # Check packages
  # ─────────────────────────────

  required_pkgs <- c(
    "xml2",
    "dplyr",
    "tibble",
    "tidygraph",
    "ggraph",
    "ggplot2"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "缺少 R 包: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }


  # ─────────────────────────────
  # Check nodes_info
  # ─────────────────────────────

  required_node_cols <- c(
    ".features_id",
    "inchikey2d",
    "synonym",
    "logFC",
    "adj.P.Val",
    "VIP"
  )

  missing_node_cols <- setdiff(
    required_node_cols,
    colnames(nodes_info)
  )

  if (length(missing_node_cols) > 0) {
    stop(
      "nodes_info 缺少必要列: ",
      paste(missing_node_cols, collapse = ", ")
    )
  }


  # ─────────────────────────────
  # Clean nodes_info
  # ─────────────────────────────

  nodes_info <- nodes_info |>
    dplyr::mutate(
      inchikey2d = toupper(inchikey2d)
    ) |>
    dplyr::filter(
      !is.na(inchikey2d),
      inchikey2d != ""
    )


  # ─────────────────────────────
  # KEGG pathway ID
  # ─────────────────────────────

  pathway_id <- sub(
    "^path:",
    "",
    pathway_id
  )

  pathway_id <- sub(
    paste0("^", species),
    "",
    pathway_id
  )

  if (!grepl("^\\d{5}$", pathway_id)) {
    stop(
      "pathway_id 应该类似 '00330'、'hsa00330' 或 'path:hsa00330'"
    )
  }

  pathway_name <- paste0(
    species,
    pathway_id
  )


  # ─────────────────────────────
  # Download KGML
  # ─────────────────────────────

  kgml_url <- paste0(
    "https://rest.kegg.jp/get/",
    pathway_name,
    "/kgml"
  )

  kgml <- xml2::read_xml(
    kgml_url
  )


  # ─────────────────────────────
  # Extract entries
  # ─────────────────────────────

  entries <- xml2::xml_find_all(
    kgml,
    ".//entry"
  )

  entry_tbl <- tibble::tibble(
    id = as.integer(
      xml2::xml_attr(
        entries,
        "id"
      )
    ),
    type = xml2::xml_attr(
      entries,
      "type"
    ),
    name = xml2::xml_attr(
      entries,
      "name"
    )
  )


  # ─────────────────────────────
  # KEGG compounds
  # ─────────────────────────────

  compound_tbl <- entry_tbl |>
    dplyr::filter(
      type == "compound"
    ) |>
    dplyr::mutate(
      kegg = sub(
        "^cpd:",
        "",
        name
      )
    ) |>
    dplyr::select(
      id,
      kegg
    ) |>
    dplyr::distinct()


  # ─────────────────────────────
  # Get KEGG InChIKey
  # ─────────────────────────────

  kegg_compounds <- unique(
    compound_tbl$kegg
  )

  get_kegg_inchikey <- function(kegg_id) {

    url <- paste0(
      "https://rest.kegg.jp/get/",
      kegg_id
    )

    txt <- tryCatch(
      readLines(
        url,
        warn = FALSE
      ),
      error = function(e) {
        character()
      }
    )

    if (length(txt) == 0) {
      return(
        tibble::tibble(
          kegg = kegg_id,
          inchikey2d = NA_character_
        )
      )
    }

    inchikey <- txt[
      grepl(
        "^\\s*DBLINKS",
        txt
      ) |
        grepl(
          "^\\s+InChIKey:",
          txt
        )
    ]

    inchikey <- paste(
      inchikey,
      collapse = " "
    )

    inchikey <- sub(
      ".*InChIKey:\\s*",
      "",
      inchikey
    )

    inchikey <- trimws(
      inchikey
    )

    if (
      identical(inchikey, "") ||
      !grepl(
        "^[A-Z]{14}-",
        inchikey
      )
    ) {
      inchikey <- NA_character_
    }

    inchikey2d <- ifelse(
      is.na(inchikey),
      NA_character_,
      sub(
        "-.*$",
        "",
        inchikey
      )
    )

    tibble::tibble(
      kegg = kegg_id,
      inchikey2d = toupper(
        inchikey2d
      )
    )
  }


  kegg_info <- dplyr::bind_rows(
    lapply(
      kegg_compounds,
      get_kegg_inchikey
    )
  )


  # ─────────────────────────────
  # KEGG compound → InChIKey
  # ─────────────────────────────

  compound_tbl <- compound_tbl |>
    dplyr::left_join(
      kegg_info,
      by = "kegg"
    )


  # ─────────────────────────────
  # Extract relations
  # ─────────────────────────────

  relations <- xml2::xml_find_all(
    kgml,
    ".//relation"
  )

  relation_tbl <- tibble::tibble(
    source_id = as.integer(
      xml2::xml_attr(
        relations,
        "entry1"
      )
    ),
    target_id = as.integer(
      xml2::xml_attr(
        relations,
        "entry2"
      )
    ),
    type = xml2::xml_attr(
      relations,
      "type"
    )
  )


  # ─────────────────────────────
  # Compound → Compound edges
  # ─────────────────────────────

  edge_tbl <- relation_tbl |>
    dplyr::inner_join(
      compound_tbl |>
        dplyr::select(
          source_id = id,
          source = kegg
        ),
      by = "source_id"
    ) |>
    dplyr::inner_join(
      compound_tbl |>
        dplyr::select(
          target_id = id,
          target = kegg
        ),
      by = "target_id"
    ) |>
    dplyr::select(
      source,
      target,
      type
    ) |>
    dplyr::distinct()


  # ─────────────────────────────
  # Build nodes
  # ─────────────────────────────

  node_tbl <- compound_tbl |>
    dplyr::select(
      name = kegg,
      inchikey2d
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(
      nodes_info,
      by = "inchikey2d"
    ) |>
    dplyr::mutate(
      logFC = dplyr::coalesce(
        logFC,
        0
      ),
      adj.P.Val = dplyr::coalesce(
        adj.P.Val,
        1
      ),
      VIP = dplyr::coalesce(
        VIP,
        0
      ),
      node_fill = factor(
        dplyr::case_when(

          logFC > 0.5 &
            adj.P.Val < 0.05 ~ "up",

          logFC < -0.5 &
            adj.P.Val < 0.05 ~ "down",

          VIP > 1 &
            adj.P.Val < 0.05 ~ "vip_high",

          TRUE ~ "no significance"
        ),
        levels = c(
          "up",
          "down",
          "vip_high",
          "no significance"
        )
      )
    )


  # ─────────────────────────────
  # Label
  # ─────────────────────────────

  if (label == "name") {

    node_tbl <- node_tbl |>
      dplyr::mutate(
        node_label = name
      )

  } else {

    node_tbl <- node_tbl |>
      dplyr::mutate(
        node_label = dplyr::if_else(
          is.na(synonym) |
            synonym == "",
          name,
          synonym
        )
      )
  }


  # ─────────────────────────────
  # Build graph
  # ─────────────────────────────

  g_tbl <- tidygraph::tbl_graph(
    nodes = node_tbl,
    edges = edge_tbl,
    directed = FALSE
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

  p <- ggraph::ggraph(
    layout_data
  ) +

    ggraph::geom_edge_link(
      ggplot2::aes(
        color = type
      ),
      alpha = 0.45,
      width = 0.7
    ) +

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
      breaks = c(
        1,
        1.5,
        2,
        3
      ),
      name = "VIP"
    ) +

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

    ggplot2::scale_fill_manual(
      values = c(
        "up" = "#D73027",
        "down" = "#4575B4",
        "vip_high" = "#7B61A8",
        "no significance" = "#BDBDBD"
      ),
      name = "Metabolite status"
    ) +

    ggraph::scale_edge_color_discrete(
      name = "KEGG relation"
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
        15,
        15,
        15,
        15
      )
    )


  # ─────────────────────────────
  # Return
  # ─────────────────────────────

  return(
    list(
      plot = p,
      nodes = node_tbl,
      edges = edge_tbl,
      kgml = kgml,
      kegg_info = kegg_info
    )
  )
}
