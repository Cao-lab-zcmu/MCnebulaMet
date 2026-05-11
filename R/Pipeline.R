run_level2 <- function(sample, lib, ionmode) {
  sample_spectra <- sample$spectra
  lib_spectra <- lib$spectra
  sample_peaks <- split(
    sample$peaks[, c("mz", "rel.int.")],
    sample$peaks$.features_id
  )
  lib_peaks <- split(
    lib$peaks[, c("mz", "rel_intensity")],
    lib$peaks$id
  )

  candidates <- generate_candidates(
    sample_spectra, lib_spectra,
    tol_ppm = 20,
    tol_rt = NULL,
    use_adduct = TRUE,
    bin_size = 0.01
  )
  res <- BiocParallel::bplapply(seq_len(nrow(candidates)), function(i) {
    fid <- candidates$.features_id[i]
    lid <- candidates$id[i]
    
    q_raw <- sample_peaks[[as.character(fid)]]
    l_raw <- lib_peaks[[as.character(lid)]]
    q <- .preprocess_spectrum(q_raw, "rel.int.", top_n = 50)
    l <- .preprocess_spectrum(l_raw, "rel_intensity", top_n = 50)

    if (is.null(q) || is.null(l)) return(NULL)

    score <- calc_cosine_fast(q, l, ppm = 50, min_match = 3)

    if (!is.na(score) && score >= 0.7) {
      cbind(candidates[i, ], cosine = score)
    } else {
      NULL
    }
  })
  results <- dplyr::bind_rows(res)
  results_dedup <- results |>
    dplyr::mutate(
      dedup_key = ifelse(
        !is.na(inchikey) & inchikey != "", 
        paste0("IK", inchikey),
        paste0("NP", name, precursor_mz)
      )
    ) |>
    dplyr::group_by(.features_id, dedup_key) |>
    dplyr::slice_max(cosine, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-dedup_key)
}

run_level1 <- function(sample, std_mix, ionmode) {
  sample_spectra <- sample$spectra
  std_mix_spectra <- std_mix$spectra
  sample_peaks <- split(
    sample$peaks[, c("mz", "rel.int.")],
    sample$peaks$.features_id
  )
  std_mix_peaks <- split(
    std_mix$peaks[, c("mz", "rel_intensity")],
    std_mix$peaks$id
  )

  candidates <- generate_candidates(
    sample_spectra, std_mix_spectra,
    tol_ppm = 20,
    tol_rt = 30,
    use_adduct = FALSE,
    bin_size = 0.01
  )

  res <- BiocParallel::bplapply(seq_len(nrow(candidates)), function(i) {
    fid <- candidates$.features_id[i]
    lid <- candidates$id[i]
    q_raw <- sample_peaks[[as.character(fid)]]
    l_raw <- std_mix_peaks[[as.character(lid)]]
    q <- .preprocess_spectrum(q_raw, "rel.int.", top_n = 50)
    l <- .preprocess_spectrum(l_raw, "rel_intensity", top_n = 50)

    if (is.null(q) || is.null(l)) return(NULL)
    score <- calc_cosine_fast(q, l, ppm = 50, min_match = 3)
    cbind(candidates[i, ], cosine = score)
  })
  results <- dplyr::bind_rows(res)
}

build_stdmix_feature_mgf <- function(
  files,
  output_file,
  ion_mode = c("pos", "neg"),
  snthresh = 10,
  noise = 1000,
  peakwidth = c(10, 40),
  ppm = 20,
  bw = 5,
  expandRt = 3,
  expandMz = 0.01,
  centroided = TRUE
) {

  ion_mode <- match.arg(ion_mode)

  ## ---- metadata (all belong to same std mix) ----
  metadata <- data.frame(
    file = files,
    sample = paste0("stdmix_", seq_along(files)),
    group = "std_mix",
    sample_type = "standard",
    ion_mode = ion_mode,
    stringsAsFactors = FALSE
  )

  ## ---- read data ----
  raw_data <- MSnbase::readMSData(
    files = metadata$file,
    centroided. = centroided,
    mode = "onDisk",
    pdata = new("AnnotatedDataFrame", metadata)
  )

  ## ---- peak detection ----
  cwp <- xcms::CentWaveParam(
    snthresh = snthresh,
    noise = noise,
    peakwidth = peakwidth,
    ppm = ppm
  )

  pro_data <- xcms::findChromPeaks(raw_data, param = cwp)

  ## ---- grouping (ALL files same group) ----
  pdp <- xcms::PeakDensityParam(
    sampleGroups = rep(1, length(files)),
    bw = bw,
    minFraction = 1,      # 混标必须全部存在
    ppm = ppm
  )

  grouped_data <- xcms::groupChromPeaks(pro_data, param = pdp)

  ## ---- extract MS2 per feature ----
  spectra <- xcms::featureSpectra(
    grouped_data,
    msLevel = 2L,
    ppm = ppm,
    expandRt = expandRt,
    expandMz = expandMz,
    return.type = "MSpectra"
  )

  ## ---- clean ----
  spectra <- MSnbase::clean(spectra, all = TRUE)

  ## ---- consensus combine (technical replicate merge) ----
  spectra_consensus <- MSnbase::combineSpectra(
    spectra,
    fcol = "feature_id",
    method = MSnbase::consensusSpectrum,
    ppm = ppm,
    mzd = expandMz
  )

  ## ---- write mgf ----
  MSnbase::writeMgfData(
    spectra_consensus,
    output_file
  )

  invisible(spectra_consensus)
}

#files_pos <- list.files(
#  "data_raw_stdmix",
#  pattern = "3 ul pos\\.mzML$",
#  full.names = TRUE
#)
#
#build_stdmix_feature_mgf(
#  files = files_pos,
#  output_file = "data_processed_stdmix/stdmix_3ul_pos.mgf",
#  ion_mode = "pos"
#)
