QC_report <- function(
  intensity_df,
  sample_info,
  output = "QC_report.pdf"
) {

  if (missing(intensity_df)) {
    stop("intensity_df is required")
  }
  if (missing(sample_info)) {
    stop("sample_info is required")
  }

  rmd_path <- system.file(
    "extdata",
    "QC.rmd",
    package = "MCnebulaMet"
  )

  if (rmd_path == "") {
    stop("Cannot find QC.rmd in package")
  }

  rmarkdown::render(
    input = rmd_path,
    output_file = output,
    params = list(
      intensity_df = intensity_df,
      sample_info = sample_info
    ),
    envir = new.env(parent = globalenv())
  )

  message("Report generated: ", output)
}
