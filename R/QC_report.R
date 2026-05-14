QC_report <- function(
  quant,
  sample_info,
  output = "QC_report.pdf"
) {

  if (missing(quant)) {
    stop("quant is required")
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

  out_dir <- dirname(output)

  rmarkdown::render(
    input = rmd_path,
    output_file = basename(output),
    output_dir = out_dir,
    params = list(
      quant = quant,
      sample_info = sample_info
    ),
    envir = new.env(parent = globalenv())
  )
  message("Report generated: ", output)
}
