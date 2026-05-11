.normalize_adduct <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- gsub("\\[|\\]", "", x)
  dplyr::case_when(
    x %in% c("M-H", "M-H-") ~ "M-H",
    x %in% c("M+H", "M+H+") ~ "M+H",
    x %in% c("M+Na", "M+Na+") ~ "M+Na",
    x %in% c("M+K", "M+K+") ~ "M+K",
    x %in% c("M+Cl", "M+Cl-") ~ "M+Cl",
    TRUE ~ x
  )
}
