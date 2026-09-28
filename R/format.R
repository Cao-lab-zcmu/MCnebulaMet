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
