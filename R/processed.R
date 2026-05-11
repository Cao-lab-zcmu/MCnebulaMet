assignAnnotationLevel4 <- function(sample_pos, sample_neg) {
  
  anno_level4 <- rbind(sample_pos, sample_neg) |> 
    dplyr::mutate(ID = paste(.features_id, ionmode, sep = "_")) |> 
    dplyr::arrange(ID, desc(tani.score)) |>
    dplyr::distinct(ID, .keep_all = TRUE) |>
    dplyr::select(
      ID = ID,
      Synonym = synonym,
      `Precursor m/z` = mz,
      `Mass error(ppm)` = error.mass,
      `Fragment m/z` = frag_mz,
      `RT(s)` = rt.secound,
      Formula = mol.formula,
      Adduct = adduct,
      `Tanimoto similarity` = tani.score,
      `InChIKey planar` = inchikey2d,
      `log2(FC)` = logFC,
      `P-Value` = P.Value,
      `FDR` = adj.P.Val
    )
  return(anno_level4)
}

assignAnnotationLevel3 <- function(sample_pos, sample_neg) {

  anno_level3 <- rbind(sample_pos, sample_neg) |>
    #dplyr::filter(!is.na(inchikey2d), !is.na(tani.score)) |>
    #dplyr::arrange(inchikey2d, desc(tani.score)) |>
    #dplyr::distinct(inchikey2d, .keep_all = TRUE) |>
    dplyr::filter(
      #rank.formula == 1, zodiac.score >= .9, 
      tani.score >= .7
    ) |> 
    dplyr::mutate(ID = paste(.features_id, ionmode, sep = "_")) |> 
    dplyr::select(
      ID = ID,
      Synonym = synonym,
      `Precursor m/z` = mz,
      `Mass error(ppm)` = error.mass,
      `Fragment m/z` = frag_mz,
      `RT(s)` = rt.secound,
      Formula = mol.formula,
      Adduct = adduct,
      `Tanimoto similarity` = tani.score,
      `InChIKey planar` = inchikey2d,
      Class = lipidClass,
      `log2(FC)` = logFC,
      `P-Value` = P.Value,
      `FDR` = adj.P.Val
    )
  return(anno_level3)
}

assignAnnotationLevel2 <- function(anno_level2_pos, anno_level2_neg) {
  anno_level2 <- rbind(anno_level2_neg, anno_level2_pos) |>
    dplyr::mutate(
      compound_key = dplyr::case_when(
        !is.na(inchikey) & inchikey != "" ~
          paste0("IK_", inchikey),
        TRUE ~
          paste0("NP_", name, "_", round(precursor_mz, 4))
      ),
      InChIKey2D = substr(inchikey, 1, 14),
      Adduct = toupper(adduct)
    ) |>
    dplyr::group_by(compound_key) |>
    dplyr::slice_max(cosine, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::group_by(ID) |>
    dplyr::slice_max(cosine, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      `InChIKey planar` = substr(inchikey, 1, 14)
    ) |> 
    dplyr::select(
      `Library ID` = id,
      ID,
      `Library synonym` = name,
      `Precursor m/z` = precursor_mz,
      `RT(s)` = rt.secound,
      Adduct,
      Formula = formula,
      `InChIKey`= inchikey,
      `InChIKey planar`,
      `Cosine score` = cosine
    )
  return(anno_level2)
}

assignAnnotationLevel1 <- function(
  anno_level1_pos,
  anno_level1_neg, std_mix_info_file
) {
  anno_level1 <- dplyr::bind_rows(
    anno_level1_pos,
    anno_level1_neg
  ) |>
    dplyr::mutate(
      ID = paste(.features_id, ionmode, sep = "_")
    ) |>
    dplyr::select(
      ID,
      ionmode, 
      `Standard m/z` = precursor_mz,
      `Sample m/z` = mz,
      `Mass error(ppm)` = error.mass,
      `Standard RT(s)` = rt,
      `Sample RT(s)` = rt.secound,
      `RT Error(s)` = rt_error,
      `Cosine score` = cosine
    )

  std_mix_meta <- readxl::read_xlsx(std_mix_info_file)
  std_mix_meta_long <- std_mix_meta |>
  tidyr::crossing(
    tibble::tibble(
      ionmode = c(
        "Positive","Positive","Positive",
        "Negative","Negative"
      ),
      adduct = c(
        "M+H","M+Na","M+K",
        "M-H","M+Cl"
      ),
      shift = c(
        1.007276,
        22.989218,
        38.963158,
        -1.007276,
        34.969402
      )
    )
  ) |>
  dplyr::mutate(
    `Theoretical m/z` = `Monoisotopic Mass` + shift
  )

  anno_level1 <- anno_level1 |>
    dplyr::inner_join(
      std_mix_meta_long,
      by = "ionmode"
    ) |>
    dplyr::mutate(
      `Mass error(Std-Theo, ppm)` = abs(`Standard m/z` - `Theoretical m/z`) / `Theoretical m/z` * 1e6
    ) |>
    dplyr::filter(`Mass error(Std-Theo, ppm)` <= 20) |> 
    dplyr::mutate(
      adduct_priority = dplyr::case_when(
        adduct == "M+H" ~ 1,
        adduct == "M-H" ~ 1,
        adduct == "M+Na" ~ 2,
        adduct == "M+K" ~ 3,
        adduct == "M+Cl" ~ 4,
        TRUE ~ 5
      )
    ) |>
    dplyr::arrange(
      Synonym,
      adduct_priority,
      dplyr::desc(`Cosine score`)
    ) |>
    dplyr::group_by(Synonym, ionmode) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::group_by(ID) |>
    dplyr::slice_max(`Cosine score`, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |> 
    dplyr::select(
      ID, `CAS number`, Synonym, `Standard m/z`,
      `Sample m/z`, `Mass error(Std-Sample, ppm)` = `Mass error(ppm)`,
      `Standard RT(s)`, `Sample RT(s)`,
      `RT Error(s)`, `InChIKey planar`, `Cosine score`,
      `Monoisotopic Mass`, Adduct = adduct, `Theoretical m/z`,
      `Mass error(Std-Theo, ppm)`
    )
  return(anno_level1)
}

generate_candidates <- function(
  sample_spectra,
  lib_spectra,
  tol_ppm = 20,
  tol_rt = NULL,
  use_adduct = FALSE,
  bin_size = 0.01
) {

  max_mz <- 1000
  max_abs_error <- max_mz * tol_ppm * 1e-6
  expand_n <- ceiling(max_abs_error / bin_size)

  sample_spectra <- sample_spectra |>
    dplyr::mutate(
      mz_bin = floor(mz / bin_size)
    )
  lib_spectra <- lib_spectra |>
    dplyr::mutate(
      mz_bin = floor(precursor_mz / bin_size)
    )

  by_cols <- c("mz_bin")
  if (use_adduct) {
    sample_spectra <- sample_spectra |> dplyr::mutate(adduct = tolower(adduct))
    lib_spectra <- lib_spectra |> dplyr::mutate(adduct = tolower(adduct))
    by_cols <- c(by_cols, "adduct", "ionmode")
  }

  sample_spectra_expanded <- purrr::map_dfr(
    -expand_n:expand_n,
    function(k) {
      sample_spectra |>
        dplyr::mutate(mz_bin = mz_bin + k)
    }
  )

  candidates <- dplyr::inner_join(
    sample_spectra_expanded, lib_spectra,
    by = by_cols, multiple = "all"

  )

  if (nrow(candidates) == 0) return(NULL)

  candidates <- candidates |>
    dplyr::mutate(
      ppm_diff = abs(mz - precursor_mz) / precursor_mz * 1e6
    ) |>
    dplyr::filter(ppm_diff <= tol_ppm)

  if (nrow(candidates) == 0) return(NULL)

  if (!is.null(tol_rt)) {
    candidates <- candidates |>
      dplyr::mutate(
        rt_error = rt.secound - rt
      ) |>
      dplyr::filter(abs(rt_error) <= tol_rt)

    if (nrow(candidates) == 0) return(NULL)
  }

  return(candidates)
}


.preprocess_spectrum <- function(df, intensity_col, top_n = 50) {

  if (is.null(df) || nrow(df) == 0) return(NULL)

  df <- df |>
    dplyr::mutate(
      intensity = sqrt(.data[[intensity_col]])
    ) |>
    dplyr::arrange(desc(intensity)) |>
    dplyr::slice_head(n = top_n)

  return(
    df[, c("mz", "intensity")]
  )
}

calc_cosine_fast <- function(q, l, ppm = 50, min_match = 6) {
  if (is.null(q) || is.null(l)) return(NA_real_)

  if (nrow(q) < min_match || nrow(l) < min_match) return(NA_real_)

  q <- q[order(q$mz), ]
  l <- l[order(l$mz), ]

  i <- j <- 1L
  score_q <- score_l <- numeric(0)

  while (i <= nrow(q) && j <= nrow(l)) {

    dmz <- q$mz[i] - l$mz[j]
    tol <- q$mz[i] * ppm * 1e-6

    if (abs(dmz) <= tol) {

      score_q <- c(score_q, q$intensity[i])
      score_l <- c(score_l, l$intensity[j])

      i <- i + 1L
      j <- j + 1L

    } else if (dmz < 0) {

      i <- i + 1L

    } else {

      j <- j + 1L
    }
  }

  if (length(score_q) < min_match) return(NA_real_)

  sum(score_q * score_l) /
    sqrt(sum(score_q^2) * sum(score_l^2))
}

calc_monoisotopic_mass <- function(formula) {

  parse_formula <- function(f) {

    matches <- gregexpr("([A-Z][a-z]*)([0-9]*)", f, perl = TRUE)
    tokens <- regmatches(f, matches)[[1]]

    elements <- sub("([A-Z][a-z]*)([0-9]*)", "\\1", tokens)
    counts   <- sub("([A-Z][a-z]*)([0-9]*)", "\\2", tokens)

    counts[counts == ""] <- 1
    counts <- as.numeric(counts)

    data.frame(
      element = elements,
      count = counts,
      stringsAsFactors = FALSE
    )
  }

  comp <- parse_formula(formula)

  mass <- sum(
    mono_mass[comp$element] * comp$count,
    na.rm = TRUE
  )

  return(mass)
}

mono_mass <- c(
  H  = 1.00782503223,
  C  = 12.0000000,
  N  = 14.00307400443,
  O  = 15.99491461957,
  P  = 30.97376199842,
  S  = 31.9720711744,
  Cl = 34.968852682,
  Na = 22.98976928,
  K = 38.96370668,
  Br = 78.9183376,
  F  = 18.998403163,
  I  = 126.9044719
)

match_level <- function(df, ref) {

  hit <- mapply(function(id, ik) {

    ref_rows <- ref[ref$ID == id, ]

    if (nrow(ref_rows) == 0) return(FALSE)

    any(
      is.na(ik) | ik == "" |
      is.na(ref_rows$`InChIKey planar`) |
      ref_rows$`InChIKey planar` == "" |
      ref_rows$`InChIKey planar` == ik
    )

  }, df$ID, df$`InChIKey planar`)

  unique(df$ID[hit])
}

get_mirror_data_lib <- function(i){

  row <- anno_sub[i,]

  feature_full <- row$ID
  lib_id <- row$`Library ID`

  feature_id <- sub("_.*","",feature_full)
  ionmode <- sub(".*_","",feature_full)

  # sample spectrum
  if(ionmode=="Positive"){
    spec_sample <- sample_pos$peaks |>
      dplyr::filter(.features_id==feature_id)
  }else{
    spec_sample <- sample_neg$peaks |>
      dplyr::filter(.features_id==feature_id)
  }

  # library spectrum
  spec_lib <- gnps_lib$peaks |>
    dplyr::filter(id==lib_id)

  list(
    sample = spec_sample,
    library = spec_lib,
    meta = row
  )
}

plot_mirror_lib <- function(i){

  dat <- get_mirror_data_lib(i)

  spec_sample <- dat$sample
  spec_lib <- dat$library
  meta <- dat$meta

  ggplot() +
    geom_segment(
      data = spec_sample,
      aes(x=mz,xend=mz,y=0,yend=rel.int.),
      linewidth=0.3,
      colour="#1f78b4"
    ) +
    geom_segment(
      data = spec_lib,
      aes(x=mz,xend=mz,y=0,yend=-rel_intensity),
      linewidth=0.3,
      colour="#e31a1c"
    ) +
    theme_classic(base_size = 11) +
    labs(
      title = meta[["Library synonym"]][1],
      subtitle = paste(
        meta[["ID"]][1],
        "Cosine =", round(meta[["Cosine score"]][1],3)
      ),
      x="m/z",
      y="Intensity"
    )
}

get_mirror_data_std <- function(i, anno, sample, lib){

  meta <- anno |>
    dplyr::filter(ID == i) |>
    dplyr::slice_max(cosine)

  if(nrow(meta) == 0) return(NULL)

  lib_id <- meta$id

  if(length(lib_id) == 0) return(NULL)

  spec_lib <- lib$peaks |>
    dplyr::filter(id %in% lib_id)

  feat <- sub("_.*","", meta$ID)

  spec_sample <- sample$peaks |>
    dplyr::filter(.features_id %in% feat)

  list(
    sample = spec_sample,
    library = spec_lib,
    meta = meta
  )
}

plot_mirror_std <- function(i, anno, sample, lib){

  dat <- get_mirror_data_std(i, anno, sample, lib)

  if(is.null(dat)) return(NULL)

  spec_sample <- dat$sample
  spec_lib <- dat$library
  meta <- dat$meta

  ggplot2::ggplot() +

    ggplot2::geom_segment(
      data = spec_sample,
      ggplot2::aes(x=mz,xend=mz,y=0,yend=rel.int.),
      linewidth=0.3,
      colour="#1f78b4"
    ) +

    ggplot2::geom_segment(
      data = spec_lib,
      ggplot2::aes(x=mz,xend=mz,y=0,yend=-rel_intensity),
      linewidth=0.3,
      colour="#e31a1c"
    ) +

    ggplot2::theme_classic(base_size = 11) +

    ggplot2::labs(
      title = meta$synonym,
      subtitle = paste0(meta$ID," Cosine=",round(meta$cosine,3)),
      x = "m/z",
      y = "Intensity"
    )
}
