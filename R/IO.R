read_mgf <- function(
  file,
  id_prefix = "MGF",
  calc_rel_int = TRUE,
  ms_level = NULL,
  progress = interactive()
) {

  ## ---- argument check ----
  stopifnot(is.character(file), length(file) == 1)
  stopifnot(is.character(id_prefix), length(id_prefix) == 1)
  stopifnot(is.logical(calc_rel_int), length(calc_rel_int) == 1)
  stopifnot(is.logical(progress), length(progress) == 1)

  if (!is.null(ms_level)) {
    stopifnot(is.numeric(ms_level), length(ms_level) == 1)
  }

  lines <- readLines(file, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  rec_starts <- which(lines == "BEGIN IONS")
  rec_ends   <- which(lines == "END IONS")

  if (length(rec_starts) == 0) {
    stop("No MGF records found.")
  }

  if (length(rec_starts) != length(rec_ends)) {
    stop("Unmatched BEGIN IONS / END IONS.")
  }

  n <- length(rec_starts)

  spectra_list <- vector("list", n)
  peaks_list   <- vector("list", n)

  if (progress) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  keep_index <- 0

  for (i in seq_len(n)) {

    if (progress) setTxtProgressBar(pb, i)

    block <- lines[(rec_starts[i] + 1):(rec_ends[i] - 1)]

    ## ---- key=value ----
    kv_lines <- block[grepl("=", block)]

    if (length(kv_lines)) {
      kv <- data.table::tstrsplit(kv_lines, "=", fixed = TRUE, keep = 1:2)
      keys   <- trimws(kv[[1]])
      values <- trimws(kv[[2]])
      fields <- setNames(values, keys)
    } else {
      fields <- character()
    }

    this_ms_level <- suppressWarnings(as.integer(fields["MSLEVEL"]))

    ## ---- filter by ms_level ----
    if (!is.null(ms_level) && this_ms_level != ms_level) {
      next
    }

    keep_index <- keep_index + 1
    id <- sprintf("%s_%05d", id_prefix, keep_index)

    spectra_list[[keep_index]] <- data.table::data.table(
      id           = id,
      feature_id   = fields["FEATURE_ID"],
      precursor_mz = suppressWarnings(as.numeric(fields["PEPMASS"])),
      charge       = suppressWarnings(as.numeric(gsub("\\+", "", fields["CHARGE"]))),
      rt           = suppressWarnings(as.numeric(fields["RTINSECONDS"])),
      ms_level     = this_ms_level,
      filename     = fields["FILENAME"],
      scans        = fields["SCANS"],
      spectype     = fields["SPECTYPE"]
    )

    ## ---- peaks ----
    peak_lines <- block[!grepl("=", block)]
    peak_lines <- peak_lines[grepl("^\\d", peak_lines)]

    if (length(peak_lines)) {

      dt <- data.table::fread(
        text = paste(peak_lines, collapse = "\n"),
        col.names = c("mz", "intensity"),
        showProgress = FALSE
      )

      dt[, id := id]

      if (calc_rel_int && nrow(dt) > 0) {
        max_int <- max(dt$intensity, na.rm = TRUE)
        if (is.finite(max_int) && max_int > 0) {
          dt[, rel_intensity := 100 * intensity / max_int]
        } else {
          dt[, rel_intensity := NA_real_]
        }
      }

      peaks_list[[keep_index]] <- dt
    }
  }

  spectra_dt <- data.table::rbindlist(
    spectra_list[seq_len(keep_index)],
    fill = TRUE
  )

  peaks_dt <- data.table::rbindlist(
    peaks_list[seq_len(keep_index)],
    fill = TRUE
  )

  list(
    spectra = spectra_dt,
    peaks   = peaks_dt
  )
}

read_msp <- function(
  file,
  id_prefix = "MSP",
  calc_rel_int = TRUE,
  progress = interactive()
) {

  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(lines)]

  rec_starts <- grep("^(NAME|Name):", lines)
  rec_starts <- c(rec_starts, length(lines) + 1)

  n <- length(rec_starts) - 1

  spectra_list <- vector("list", n)
  peaks_list   <- vector("list", n)

  if (progress) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {

    if (progress) setTxtProgressBar(pb, i)

    block <- lines[rec_starts[i]:(rec_starts[i + 1] - 1)]
    id <- sprintf("%s_%05d", id_prefix, i)

    ## ---- fields ----
    kv_lines <- block[grepl(":", block) & !grepl("^\\d", block)]

    kv <- data.table::tstrsplit(kv_lines, ":", fixed = TRUE, keep = 1:2)

    fields <- setNames(trimws(kv[[2]]), toupper(kv[[1]]))

    ## ---- parse comments ----
    comments <- fields["COMMENTS"]

    smiles <- fields["SMILES"]

    if (is.na(smiles) && !is.na(comments)) {
      m <- regmatches(comments, regexpr('SMILES=[^"]+', comments))
      if (length(m)) smiles <- sub("SMILES=", "", m)
    }

    cas <- NA
    if (!is.na(comments)) {
      m <- regmatches(comments, regexpr('cas number=[^"]+', comments))
      if (length(m)) cas <- sub("cas number=", "", m)
    }

    ## ---- ion mode normalize ----
    ionmode <- fields["IONMODE"]
    if (is.na(ionmode)) ionmode <- fields["ION_MODE"]

    if (!is.na(ionmode)) {
      ionmode <- ifelse(
        toupper(substr(ionmode,1,1)) == "P",
        "Positive",
        "Negative"
      )
    }

    ## ---- precursor mz ----
    precursor <- fields["PRECURSORMZ"]
    if (is.na(precursor)) precursor <- fields["EXACTMASS"]

    spectra_list[[i]] <- data.table::data.table(
      id           = id,
      name         = fields[c("NAME","NAME")][1],
      precursor_mz = as.numeric(precursor),
      adduct       = fields["PRECURSORTYPE"],
      formula      = fields["FORMULA"],
      inchikey     = fields["INCHIKEY"],
      smiles       = smiles,
      cas          = cas,
      ionmode      = ionmode,
      instrument   = fields["INSTRUMENT"]
    )

    ## ---- peaks ----
    p_start <- grep("^Num Peaks:", block)

    if (length(p_start) == 1) {

      p_txt <- block[(p_start + 1):length(block)]
      p_txt <- p_txt[grepl("^\\d", p_txt)]

      if (length(p_txt)) {

        dt <- data.table::fread(
          text = paste(p_txt, collapse = "\n"),
          col.names = c("mz", "intensity"),
          showProgress = FALSE
        )

        dt[, id := id]

        if (calc_rel_int) {
          dt[, rel_intensity := 100 * intensity / max(intensity)]
        }

        peaks_list[[i]] <- dt
      }
    }
  }

  list(
    spectra = data.table::rbindlist(spectra_list, fill = TRUE),
    peaks   = data.table::rbindlist(peaks_list,   fill = TRUE)
  )
}
#read_msp <- function(
#  file,
#  id_prefix = "MSP",
#  calc_rel_int = TRUE,
#  progress = interactive()
#) {
#
#  lines <- readLines(file, warn = FALSE)
#  lines <- lines[nzchar(lines)]
#
#  rec_starts <- which(startsWith(lines, "NAME:"))
#  rec_starts <- c(rec_starts, length(lines) + 1)
#
#  n <- length(rec_starts) - 1
#
#  spectra_list <- vector("list", n)
#  peaks_list   <- vector("list", n)
#
#  ## ---- progress bar ----
#  if (progress) {
#    pb <- txtProgressBar(min = 0, max = n, style = 3)
#    on.exit(close(pb), add = TRUE)
#  }
#
#  for (i in seq_len(n)) {
#
#    if (progress) setTxtProgressBar(pb, i)
#
#    block <- lines[rec_starts[i]:(rec_starts[i + 1] - 1)]
#    id <- sprintf("%s_%05d", id_prefix, i)
#
#    ## ---- fields ----
#    kv_lines <- block[grepl(":", block) & !grepl("^\\d", block)]
#    kv <- data.table::tstrsplit(kv_lines, ":", fixed = TRUE, keep = 1:2)
#    fields <- setNames(trimws(kv[[2]]), kv[[1]])
#
#    spectra_list[[i]] <- data.table::data.table(
#      id           = id,
#      name         = fields["NAME"],
#      precursor_mz = as.numeric(fields["PRECURSORMZ"]),
#      adduct       = fields["PRECURSORTYPE"],
#      formula      = fields["FORMULA"],
#      inchikey     = fields["INCHIKEY"],
#      smiles       = fields["SMILES"],
#      ionmode      = fields["IONMODE"],
#      instrument   = fields["INSTRUMENT"]
#    )
#
#    ## ---- peaks ----
#    p_start <- which(startsWith(block, "Num Peaks:"))
#
#    if (length(p_start) == 1) {
#
#      p_txt <- block[(p_start + 1):length(block)]
#      p_txt <- p_txt[grepl("^\\d", p_txt)]
#
#      if (length(p_txt)) {
#
#        dt <- data.table::fread(
#          text = paste(p_txt, collapse = "\n"),
#          col.names = c("mz", "intensity"),
#          showProgress = FALSE
#        )
#
#        dt[, id := id]
#
#        if (calc_rel_int) {
#          dt[, rel_intensity := 100 * intensity / max(intensity)]
#        }
#
#        peaks_list[[i]] <- dt
#      }
#    }
#  }
#
#  list(
#    spectra = data.table::rbindlist(spectra_list, fill = TRUE),
#    peaks   = data.table::rbindlist(peaks_list,   fill = TRUE)
#  )
#}

read_mcn <- function(mcn, ionmode){
  anno_table <- features_annotation(mcn)
  stat_table <- top_table(statistic_set(mcn))[[1]]
  anno_table <- dplyr::mutate(
    anno_table,
    ionmode = ionmode,
    adduct = .normalize_adduct(adduct),
    ID = paste(.features_id, ionmode, sep = "_")
  )
  anno_table <- dplyr::left_join(
    anno_table,
    dplyr::select(stat_table, .features_id, logFC, P.Value, adj.P.Val),
    by = ".features_id"
  )

  sig_msms <- latest(mcn, "project_dataset", ".f3_spectra")
  sig_split <- split(
    sig_msms[, c("mz", "rel.int.")],
    sig_msms$.features_id
  )

  anno_table$frag_mz <- vapply(
  anno_table$.features_id,
    function(fid) {
        sp <- sig_split[[fid]]
        if (is.null(sp)) return(NA_character_)
        idx <- order(sp$rel.int., decreasing = TRUE)
        mz_top <- sp$mz[idx][1:min(5, length(idx))]
        paste(mz_top, collapse = ",")
      },
    character(1)
  )

  list(
    spectra = anno_table,
    peaks = sig_msms
  )
}


reformat_mgf_custom <- function(infile, outfile, force_charge = NULL) {
  lines <- readLines(infile)
  new_lines <- vector("list",  length = sum(grepl("^BEGIN IONS", lines)))
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  for (k in seq_along(begin_idx)) {
    block <- lines[begin_idx[k]:end_idx[k]]
    
    feature <- sub("^FEATURE_ID=", "", grep("^FEATURE_ID=", block, value = TRUE))
    pepmass <- grep("^PEPMASS=", block, value = TRUE)
    charge  <- grep("^CHARGE=", block, value = TRUE)
    charge_val <- if (length(charge) > 0) sub("^CHARGE=", "", charge) else ""
    rt <- grep("^RTINSECONDS", block, value = TRUE)
    
    if (!is.null(force_charge)) {
      if (force_charge == "positive") {
        charge_fmt <- "1+"
      } else if (force_charge == "negative") {
        charge_fmt <- "1-"
      } else {
        stop("force_charge must be 'positive', 'negative', or NULL")
      }
    } else {
      charge_fmt <- charge_val
      if (charge_fmt == "" || charge_fmt == "0") charge_fmt <- "1+"
    }
    
    title_line <- grep("^TITLE=", block, value = TRUE)
    mslevel <- if (length(title_line) > 0) sub(".*msLevel ([0-9]+).*", "\\1", title_line) else "NA"
    
    peaks <- block[!grepl("=", block) & !grepl("^BEGIN|^END", block)]
    scans <- grep("^SCANS", block, value = TRUE)
    
    new_lines[[k]] <- c(
      "BEGIN IONS",
      paste0("FEATURE_ID=", feature),
      pepmass,
      paste0("CHARGE=", charge_fmt),
      rt,
      paste0("MSLEVEL=", mslevel),
      scans,
      peaks,
      "END IONS",
      "",
      if (k < length(begin_idx)) "" else NULL
    )
  }
  
  writeLines(unlist(new_lines), outfile)
}

create_ST <- function(annotation_table, Level1, Level2, Level3) {
  

  Level1 <- Level1 |>
    dplyr::left_join(
      Level2 |> dplyr::select(ID, `InChIKey planar`) |> 
        dplyr::rename(`InChIKey planar level 2` = `InChIKey planar`),
      by = "ID"
    ) |>
    dplyr::left_join(
      Level3 |> dplyr::select(ID, `InChIKey planar`) |> 
        dplyr::rename(`InChIKey planar level 3` = `InChIKey planar`),
      by = "ID"
    ) |>
    dplyr::mutate(
      `In Level2` = ifelse(!is.na(`InChIKey planar level 2`), "Yes", "No"),
      `In Level3` = ifelse(!is.na(`InChIKey planar level 3`), "Yes", "No")
    )
  #Level1 <- Level1 |>
  #  dplyr::mutate(
  #    `In Level2` = ifelse(
  #      ID %in% Level2$ID,
  #      #ID %in% Level2$ID &
  #      #`InChIKey planar` %in% Level2$`InChIKey planar`,
  #      "Yes", "No"
  #    ),
  #    `In Level3` = ifelse(
  #      ID %in% Level3$ID,
  #      #ID %in% Level3$ID &
  #      #`InChIKey planar` %in% Level3$`InChIKey planar`,
  #      "Yes", "No"
  #    )
  #  )
  
  wb <- openxlsx::createWorkbook()
  
  sheets <- list(
    "Summary-Table" = annotation_table,
    "Level1-Standard-Samples" = Level1,
    "Level2-MS2-Library" = Level2,
    "Level3-In-Silico" = Level3
  )
  
  sheet_intro <- list(
    "Summary-Table" = "This sheet summarizes the annotation table for all samples.",
    "Level1-Standard-Samples" = "Level 1: Standard samples information.",
    "Level2-MS2-Library" = "Level 2: MS2 library information.",
    "Level3-In-Silico" = "Level 3: In-silico predicted compounds."
  )
  
  for (nm in names(sheets)) {
    
    df <- sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    
    intro_style <- openxlsx::createStyle(
      fontSize = 14,
      textDecoration = "Bold"
    )
    
    openxlsx::writeData(wb, nm, sheet_intro[[nm]], startCol = 1, startRow = 1)
    openxlsx::addStyle(wb, nm, style = intro_style, rows = 1, cols = 1, gridExpand = TRUE)
    
    openxlsx::writeData(
      wb, nm, df,
      startCol = 1,
      startRow = 2,
      headerStyle = openxlsx::createStyle(
        fontSize = 12,
        fontColour = "white",
        halign = "center",
        fgFill = "#4F81BD",
        border = "TopBottom",
        textDecoration = "Bold"
      )
    )
    
    openxlsx::setColWidths(wb, nm, cols = 1:ncol(df), widths = "auto")
    
    if (nm == "Level1-Standard-Samples") {
      
      yes_style <- openxlsx::createStyle(
        fontColour = "#006100",
        fgFill = "#C6EFCE"
      )
      
      no_style <- openxlsx::createStyle(
        fontColour = "#9C0006",
        fgFill = "#FFC7CE"
      )
      
      start_row <- 3
      end_row <- nrow(df) + 2
      
      cols_to_format <- (ncol(df)-1):ncol(df)
      
      for (col in cols_to_format) {
        
        openxlsx::conditionalFormatting(
          wb, nm,
          cols = col,
          rows = start_row:end_row,
          rule = '=="Yes"',
          style = yes_style
        )
        
        openxlsx::conditionalFormatting(
          wb, nm,
          cols = col,
          rows = start_row:end_row,
          rule = '=="No"',
          style = no_style
        )
      }
      openxlsx::addStyle(
        wb, nm,
        style = openxlsx::createStyle(border = "TopBottomLeftRight"),
        rows = 2:(nrow(df) + 2),
        cols = (ncol(df)-1):ncol(df),
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }
  
  openxlsx::saveWorkbook(wb, "Supplementary_Table_S1.xlsx", overwrite = TRUE)
}

save_mirror_pdf <- function(plots, file){

  plots <- Filter(Negate(is.null), plots)

  per_page <- 9
  n <- length(plots)

  pdf(file, width = 10, height = 10)

  for(i in seq(1, n, by = per_page)){

    idx <- i:min(i + per_page - 1, n)

    p <- patchwork::wrap_plots(
      plots[idx],
      ncol = 3
    )

    print(p)

  }

  dev.off()
}

create_ST2 <- function(annotation_table, Level1, Level2, Level3,
                      file_name = "Supplementary_Table_S1.xlsx") {
  
  Level1 <- Level1 |>
    dplyr::left_join(
      Level2 |> 
        dplyr::select(ID, `InChIKey planar`) |> 
        dplyr::rename(`InChIKey (Level2)` = `InChIKey planar`),
      by = "ID"
    ) |>
    dplyr::left_join(
      Level3 |> 
        dplyr::select(ID, `InChIKey planar`) |> 
        dplyr::rename(`InChIKey (Level3)` = `InChIKey planar`),
      by = "ID"
    ) |>
    dplyr::mutate(
      `Matched in Level2` = ifelse(!is.na(`InChIKey (Level2)`), "Yes", "No"),
      `Matched in Level3` = ifelse(!is.na(`InChIKey (Level3)`), "Yes", "No")
    )
  
  wb <- openxlsx::createWorkbook()
  
  header_style <- openxlsx::createStyle(
    fontSize = 11,
    fontColour = "black",
    textDecoration = "Bold",
    halign = "center",
    border = "Bottom",
    fgFill = "#D9E1F2"
  )
  
  body_style <- openxlsx::createStyle(
    fontSize = 10,
    border = "TopBottom"
  )
  
  intro_style <- openxlsx::createStyle(
    fontSize = 12,
    textDecoration = "Bold"
  )
  
  note_style <- openxlsx::createStyle(
    fontSize = 9,
    fontColour = "#555555",
    textDecoration = "Italic"
  )
  
  yes_style <- openxlsx::createStyle(
    fontColour = "#006100",
    fgFill = "#C6EFCE"
  )
  
  no_style <- openxlsx::createStyle(
    fontColour = "#9C0006",
    fgFill = "#FFC7CE"
  )
  
  sheets <- list(
    "Summary" = annotation_table,
    "Level1_Standard" = Level1,
    "Level2_MS2" = Level2,
    "Level3_InSilico" = Level3
  )
  
  sheet_intro <- list(
    "Summary" = "Comprehensive annotation summary across all samples.",
    "Level1_Standard" = "Level 1 annotation: confirmed using authentic standards.",
    "Level2_MS2" = "Level 2 annotation: MS/MS spectral library matching.",
    "Level3_InSilico" = "Level 3 annotation: in silico predicted compounds."
  )
  
  sheet_note <- list(
    "Summary" = "All annotations follow Metabolomics Standards Initiative (MSI).",
    "Level1_Standard" = "Level 1 = confirmed structure (RT + MS/MS + standard).",
    "Level2_MS2" = "Level 2 = putatively annotated compounds.",
    "Level3_InSilico" = "Level 3 = tentative candidates based on computational prediction."
  )
  
  for (nm in names(sheets)) {
    
    df <- sheets[[nm]]
    
    openxlsx::addWorksheet(wb, nm)
    
    openxlsx::writeData(wb, nm, sheet_intro[[nm]], startRow = 1)
    openxlsx::addStyle(wb, nm, intro_style, rows = 1, cols = 1)
    
    openxlsx::writeData(
      wb, nm, df,
      startRow = 3,
      headerStyle = header_style,
      borders = "rows"
    )
    
    openxlsx::setColWidths(wb, nm, cols = 1:ncol(df), widths = "auto")
    
    openxlsx::freezePane(wb, nm, firstRow = TRUE, firstCol = TRUE)
    
    if (nm == "Level1_Standard") {
      
      cols_to_format <- (ncol(df)-1):ncol(df)
      start_row <- 4
      end_row <- nrow(df) + 3
      
      for (col in cols_to_format) {
        openxlsx::conditionalFormatting(
          wb, nm,
          cols = col,
          rows = start_row:end_row,
          rule = '=="Yes"',
          style = yes_style
        )
        
        openxlsx::conditionalFormatting(
          wb, nm,
          cols = col,
          rows = start_row:end_row,
          rule = '=="No"',
          style = no_style
        )
      }
    }
    
    openxlsx::writeData(
      wb, nm,
      sheet_note[[nm]],
      startRow = nrow(df) + 5
    )
    
    openxlsx::addStyle(
      wb, nm,
      note_style,
      rows = nrow(df) + 5,
      cols = 1
    )
  }
  
  openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)
}

export_annotation_table <- function(
    annotation_table,
    file = "annotation_table.xlsx"
) {

    wb <- openxlsx::createWorkbook()
    sheet_name <- "Annotation_Table"

    openxlsx::addWorksheet(wb, sheet_name)

    intro_text <- "This table summarizes all annotated metabolites with confidence levels (Level 1–3)."

    intro_style <- openxlsx::createStyle(
        fontSize = 14,
        textDecoration = "Bold"
    )

    openxlsx::writeData(
        wb,
        sheet_name,
        intro_text,
        startRow = 1,
        startCol = 1
    )

    openxlsx::addStyle(
        wb,
        sheet_name,
        intro_style,
        rows = 1,
        cols = 1
    )

    header_style <- openxlsx::createStyle(
        fontSize = 12,
        fontColour = "white",
        halign = "center",
        fgFill = "#4F81BD",
        textDecoration = "Bold",
        border = "TopBottomLeftRight"
    )

    openxlsx::writeData(
        wb,
        sheet_name,
        annotation_table,
        startRow = 2,
        headerStyle = header_style
    )

    openxlsx::setColWidths(
        wb,
        sheet_name,
        cols = 1:ncol(annotation_table),
        widths = "auto"
    )

    openxlsx::freezePane(
        wb,
        sheet_name,
        firstActiveRow = 3
    )

    openxlsx::addFilter(
        wb,
        sheet_name,
        row = 2,
        cols = 1:ncol(annotation_table)
    )

    border_style <- openxlsx::createStyle(
        border = "TopBottomLeftRight"
    )

    openxlsx::addStyle(
        wb,
        sheet_name,
        style = border_style,
        rows = 2:(nrow(annotation_table) + 2),
        cols = 1:ncol(annotation_table),
        gridExpand = TRUE,
        stack = TRUE
    )

    level1_style <- openxlsx::createStyle(
        fgFill = "#DFF0D8",
        fontColour = "#006100"
    )

    level2_style <- openxlsx::createStyle(
        fgFill = "#D9EDF7",
        fontColour = "#1F4E79"
    )

    level3_style <- openxlsx::createStyle(
        fgFill = "#F2F2F2",
        fontColour = "#595959"
    )

    col_idx <- which(colnames(annotation_table) == "Annotation_Level")

    rows_l1 <- which(annotation_table$Annotation_Level == "Level 1") + 2
    rows_l2 <- which(annotation_table$Annotation_Level == "Level 2") + 2
    rows_l3 <- which(annotation_table$Annotation_Level == "Level 3") + 2

    if (length(rows_l1) > 0) {
        openxlsx::addStyle(
            wb,
            sheet_name,
            level1_style,
            rows = rows_l1,
            cols = col_idx,
            gridExpand = TRUE,
            stack = TRUE
        )
    }

    if (length(rows_l2) > 0) {
        openxlsx::addStyle(
            wb,
            sheet_name,
            level2_style,
            rows = rows_l2,
            cols = col_idx,
            gridExpand = TRUE,
            stack = TRUE
        )
    }

    if (length(rows_l3) > 0) {
        openxlsx::addStyle(
            wb,
            sheet_name,
            level3_style,
            rows = rows_l3,
            cols = col_idx,
            gridExpand = TRUE,
            stack = TRUE
        )
    }

    openxlsx::saveWorkbook(
        wb,
        file,
        overwrite = TRUE
    )
}
