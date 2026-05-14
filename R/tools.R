.check_feature_table <- function(x) {

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(
      "'x' must be a matrix or data.frame, got: ",
      class(x),
      call. = FALSE
    )
  }

  rn <- rownames(x)

  if (is.null(rn)) {
    stop("'x' must have rownames (feature IDs)", call. = FALSE)
  }
  if (anyNA(rn)) {
    stop("NA values found in rownames(x)", call. = FALSE)
  }
  if (anyDuplicated(rn)) {
    stop("Duplicated feature IDs in rownames(x)", call. = FALSE)
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  rownames(x) <- rn

  if (nrow(x) == 0L) {
    stop("'x' has 0 rows (no features)", call. = FALSE)
  }
  if (ncol(x) == 0L) {
    stop("'x' has 0 columns (no samples)", call. = FALSE)
  }

  if (is.null(colnames(x))) {
    stop("'x' must have colnames (sample IDs)", call. = FALSE)
  }
  if (anyNA(colnames(x))) {
    stop("NA values found in colnames(x)", call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) {
    stop("Duplicated sample IDs in colnames(x)", call. = FALSE)
  }

  if (!is.numeric(x)) {
    stop("'x' must be numeric", call. = FALSE)
  }

  if (any(x < 0, na.rm = TRUE)) {
    warning(
      "Negative intensities detected; check preprocessing",
      call. = FALSE
    )
  }

  na_count <- sum(is.na(x))
  if (na_count > 0) {
    warning(
      sprintf(
        "Detected %d NA values (%.2f%% of total)",
        na_count,
        100 * na_count / length(x)
      ),
      call. = FALSE
    )
  }

  zero_count <- sum(x == 0, na.rm = TRUE)
  if (zero_count > 0) {
    warning(
      sprintf(
        "Detected %d zero values (%.2f%% of total). In LC-MS data (e.g. MZmine export), zeros often represent missing values.",
        zero_count,
        100 * zero_count / length(x)
      ),
      call. = FALSE
    )
  }

  invisible(x)
}

#' Convert zero to NA
#'
#' @param x numeric matrix
#' @param zero_as_na logical
#'
#' @return matrix
#' @export
zero_to_na <- function(x, zero_as_na = TRUE) {

  .check_feature_table(x)

  if (!zero_as_na) {
    return(x)
  }

  zero_count <- sum(x == 0, na.rm = TRUE)

  if (zero_count > 0) {
    message("Converting ", zero_count, " zeros to NA")
    x[x == 0] <- NA
  }

  return(x)
}
