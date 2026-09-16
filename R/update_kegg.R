.kegg_rest_url <- function(operation, argument) {
  paste0(
    "https://rest.kegg.jp/",
    operation,
    "/",
    argument
  )
}


.kegg_request <- function(
    operation,
    argument,
    quiet = FALSE
) {

  url <- .kegg_rest_url(
    operation = operation,
    argument = argument
  )

  if (!quiet) {
    message("GET ", url)
  }

  response <- httr2::request(url) |>
    httr2::req_perform()

  httr2::resp_check_status(response)

  httr2::resp_body_string(response)
}


.kegg_request_batch <- function(
    ids,
    database,
    batch_size = 10,
    sleep = 0.35,
    quiet = FALSE
) {

  ids <- unique(as.character(ids))
  ids <- ids[nzchar(ids)]

  batches <- split(
    ids,
    ceiling(seq_along(ids) / batch_size)
  )

  result <- vector(
    "list",
    length(batches)
  )

  for (i in seq_along(batches)) {

    if (!quiet) {
      message(
        "Downloading ",
        database,
        ": batch ",
        i,
        "/",
        length(batches)
      )
    }

    batch_id <- paste(
      batches[[i]],
      collapse = "+"
    )

    result[[i]] <- .kegg_request(
      operation = "get",
      argument = batch_id,
      quiet = TRUE
    )

    if (i < length(batches)) {
      Sys.sleep(sleep)
    }
  }

  paste(
    result,
    collapse = "\n"
  )
}


.kegg_parse_list <- function(x) {

  lines <- strsplit(
    x,
    "\n",
    fixed = TRUE
  )[[1]]

  lines <- lines[nzchar(lines)]

  if (!length(lines)) {
    return(
      data.frame(
        id = character(),
        name = character(),
        stringsAsFactors = FALSE
      )
    )
  }

  res <- strsplit(
    lines,
    "\t",
    fixed = TRUE
  )

  data.frame(
    id = vapply(
      res,
      `[`,
      character(1),
      1
    ),
    name = vapply(
      res,
      function(x) {
        if (length(x) >= 2) {
          x[[2]]
        } else {
          NA_character_
        }
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}
