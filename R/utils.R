#' @importFrom rlang as_label
.check_data <- 
  function(object, lst, tip = "(...)"){
    target <- rlang::as_label(substitute(object))
    mapply(lst, names(lst), FUN = function(value, name){
             obj <- match.fun(name)(object)
             if (is.null(obj)) {
               stop(paste0("is.null(", name, "(", target, ")) == T. ",
                           "use `", value, tip, "` previously."))
             }
             if (is.list(obj)) {
               if (length(obj) == 0) {
                 stop(paste0("length(", name, "(", target, ")) == 0. ",
                             "use `", value, tip, "` previously."))
               }
             }
           })
  }

.check_names <- 
  function(param, formal, tip1, tip2){
    if (!is.null(names(param))) {
      if ( any(!names(formal) %in% names(param)) ) {
        stop(paste0("the names of `", tip1, "` must contain all names of ",
                    tip2, "; or without names."
                    ))
      }
    }
  }

#' @importFrom rlang as_label
.check_class <- 
  function(object, class = "layout", tip = "grid::grid.layout"){
    if (!is(object, class)) {
      stop(paste0("`", rlang::as_label(substitute(object)),
                  "` should be a '", class, "' object created by ",
                  "`", tip, "`." ))
    }
  }

.check_columns <- 
  function(obj, lst, tip){
    if (!is.data.frame(obj))
      stop(paste0("'", tip, "' must be a 'data.frame'."))
    lapply(lst, function(col){
             if (is.null(obj[[ col ]]))
               stop(paste0("'", tip, "' must contains a column of '", col, "'."))
           })
  }

.check_type <- 
  function(obj, type, tip){
    fun <- match.fun(paste0("is.", type))
    apply(obj, 2, function(col){
            if (!fun(col))
              stop(paste0("data columns in '", tip, "' must all be '", type, "'."))
           })
  }

.check_path <- 
  function(path){
    if (!file.exists(path)) {
      dir.create(path, recursive = T)
    }
  }

.check_file <- 
  function(file){
    if (!file.exists(file)) {
      stop("file.exists(file) == F, `file` not exists.")
    }
  }
