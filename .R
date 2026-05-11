usethis::use_inst()
dir.create("inst/extdata", recursive = TRUE)

devtools::load_all(".")
library(ggplot2)
files <- list.files(
  "~/workspace/project//kidney_fibrosis_rats_20230518/serum_uuo_byn/data_mzmine/neg",
  pattern = "MSMS\\.csv$", full.names = TRUE
)

origin <- data.table::fread(files)
origin <- as.data.frame(origin)

intensity_df <- dplyr::select(
  origin, .features_id = 1, dplyr::contains("Peak area")
)
rownames(intensity_df) <- intensity_df$.features_id
intensity_df <- dplyr::select(intensity_df, -.features_id)
colnames(intensity_df) <- stringr::str_replace_all(
  colnames(intensity_df), c("\\.mzML Peak area" = "", "_NEG" = "")
)
intensity_df <- dplyr::select(intensity_df, sort(colnames(intensity_df)))
intensity_df <- lapply(intensity_df, function(x) {
  x[x == 0] <- NA
  x
})
intensity_df <- as.data.frame(intensity_df)
sample_info <- get_metadata(
  #colnames(intensity_df), c(Sham = "^Sham", Model = "^M", QC = "^QC")
  colnames(intensity_df), c(Sham = "^Sham", Model = "^M")
)

QC_report(intensity_df, sample_info)
