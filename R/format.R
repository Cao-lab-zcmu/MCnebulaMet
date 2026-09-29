#' Create formatted metabolite table
#'
#' @param df data.frame, metabolite annotation result
#' @param fontname character, table font
#' @param fontsize numeric, table font size
#'
#' @return flextable object
#' @export
format_metabolite_table <- function(df, fontname = "Times New Roman", fontsize = 9) {
  # Prepare
  df <- df |>
    dplyr::mutate(
      mz = round(mz, 4),
      rt.min = round(rt.min, 2),
      tani.score = round(tani.score, 2),
      logFC = round(logFC, 2),
      P.Value = formatC(P.Value, format = "e", digits = 2),
      adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2),
      VIP = round(VIP, 2)
    )
  # Plot
  flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::font(fontname = fontname, part = "all") |>
    flextable::fontsize(size = fontsize, part = "all") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(
      j = c("synonym", "mol.formula"),
      align = "left",
      part = "all"
    ) |>
    flextable::set_table_properties(layout = "autofit") |>
    flextable::autofit()
}

format_metabolite_table_ppt_met <- function(
    df,
    fontname = "Times New Roman",
    fontsize = 9
) {
  # Prepare
  n_total <- nrow(df)

  df <- df |>
    dplyr::slice_head(n = 10) |>
    dplyr::mutate(
      mz = round(mz, 4),
      rt.min = round(rt.min, 2),
      tani.score = round(tani.score, 2),
      logFC = round(logFC, 2),
      P.Value = formatC(P.Value, format = "e", digits = 2),
      adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2),
      VIP = round(VIP, 2)
    )

  summary_row <- df[1, , drop = FALSE]

  summary_row[,] <- NA

  summary_row[[1]] <- paste0(
    "... (Total: ",
    n_total,
    " metabolites)"
  )

  df <- dplyr::bind_rows(
    df,
    summary_row
  )

  # Plot
  flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::font(
      fontname = fontname,
      part = "all"
    ) |>
    flextable::fontsize(
      size = fontsize,
      part = "all"
    ) |>
    flextable::align(
      align = "center",
      part = "all"
    ) |>
    flextable::align(
      j = c("synonym", "mol.formula"),
      align = "left",
      part = "all"
    ) |>
    flextable::align(
      i = nrow(df),
      align = "left",
      part = "body"
    ) |>
    flextable::italic(
      i = nrow(df),
      part = "body"
    ) |>
    flextable::set_table_properties(
      layout = "autofit"
    ) |>
    flextable::autofit()
}

format_metabolite_table2 <- function(
    df,
    fontname = "Times New Roman",
    fontsize = 9
) {
  # Prepare
  df <- df |>
    dplyr::mutate(
      `Precursor m/z` = round(`Precursor m/z`, 4),
      `RT(s)` = round(`RT(s)`, 2),
      `Cosine score` = round(`Cosine score`, 2)
    )

  # Plot
  flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::font(
      fontname = fontname,
      part = "all"
    ) |>
    flextable::fontsize(
      size = fontsize,
      part = "all"
    ) |>
    flextable::align(
      align = "center",
      part = "all"
    ) |>
    flextable::align(
      j = c(
        "Library synonym",
        "Formula",
        "InChIKey",
        "InChIKey planar"
      ),
      align = "left",
      part = "all"
    ) |>
    flextable::set_table_properties(
      layout = "autofit"
    ) |>
    flextable::autofit()
}

format_metabolite_table_ppt <- function(
    df,
    fontname = "Times New Roman",
    fontsize = 9
) {
  # Prepare
  n_total <- nrow(df)

  df <- df |>
    dplyr::slice_head(n = 10) |>
    dplyr::mutate(
      `Precursor m/z` = round(`Precursor m/z`, 4),
      `RT(s)` = round(`RT(s)`, 2),
      `Cosine score` = round(`Cosine score`, 2)
    )

  summary_row <- df[1, , drop = FALSE]

  summary_row$`Library ID` <- paste0(
    "... (Total: ",
    n_total,
    " metabolites)"
  )
  summary_row$ID <- NA_character_
  summary_row$`Library synonym` <- NA_character_
  summary_row$`Precursor m/z` <- NA_real_
  summary_row$`RT(s)` <- NA_real_
  summary_row$Adduct <- NA_character_
  summary_row$Formula <- NA_character_
  summary_row$InChIKey <- NA_character_
  summary_row$`InChIKey planar` <- NA_character_
  summary_row$`Cosine score` <- NA_real_

  df <- dplyr::bind_rows(
    df,
    summary_row
  )

  # Plot
  flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::font(
      fontname = fontname,
      part = "all"
    ) |>
    flextable::fontsize(
      size = fontsize,
      part = "all"
    ) |>
    flextable::align(
      align = "center",
      part = "all"
    ) |>
    flextable::align(
      j = c(
        "Library synonym",
        "Formula",
        "InChIKey",
        "InChIKey planar"
      ),
      align = "left",
      part = "all"
    ) |>
    flextable::align(
      i = nrow(df),
      align = "left",
      part = "body"
    ) |>
    flextable::italic(
      i = nrow(df),
      part = "body"
    ) |>
    flextable::set_table_properties(
      layout = "autofit"
    ) |>
    flextable::autofit()
}

format_metabolite_table_10_met <- function(
    df,
    fontname = "Times New Roman",
    fontsize = 9
) {
  # Prepare
  n_total <- nrow(df)

  df <- df |>
    dplyr::slice_head(n = 10) |>
    dplyr::mutate(
      `Precursor m/z` = round(`Precursor m/z`, 4),
      `RT(s)` = round(`RT(s)`, 1),
      Tani.score = round(Tani.score, 2)
    )

  summary_row <- df[1, , drop = FALSE]
  summary_row[,] <- NA
  summary_row[[1]] <- paste0(
    "... (Total: ",
    n_total,
    " metabolites)"
  )

  df <- dplyr::bind_rows(
    df,
    summary_row
  )

  # Plot
  flextable::flextable(df) |>
    flextable::theme_booktabs() |>
    flextable::font(
      fontname = fontname,
      part = "all"
    ) |>
    flextable::fontsize(
      size = fontsize,
      part = "all"
    ) |>
    flextable::align(
      align = "center",
      part = "all"
    ) |>
    flextable::align(
      j = c("synonym", "Formula", "Class"),
      align = "left",
      part = "all"
    ) |>
    flextable::align(
      i = nrow(df),
      align = "left",
      part = "body"
    ) |>
    flextable::italic(
      i = nrow(df),
      part = "body"
    ) |>
    flextable::set_table_properties(
      layout = "autofit"
    ) |>
    flextable::autofit()
}
