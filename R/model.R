#' @export
print.sleuth_model <- function(obj) {
  cat('formula: ', deparse(obj$formula), '\n')
  if (!is.null(obj$wald)) {
    cat('tests:\n')
    cat(paste0('\t', names(obj$wald)))
  } else {
    cat('no tests found.\n')
  }

  invisible(obj)
}

#' @export
models <- function(obj) {
  UseMethod('models')
}

#' @export
models.sleuth <- function(obj) {
  for (x in names(obj$fits)) {
    cat('[ ', x,' ]\n')
    models(obj$fits[[x]])
  }

  invisible(obj)
}

models.sleuth_model <- function(obj) {
  print(obj)
}

#' @export
wald_results <- function(obj, which_beta, which_model = 'full') {
  stopifnot( is(obj, 'sleuth') )

  if ( !(which_beta %in% names(obj$fits[[which_model]]$wald)) ) {
    stop("'", which_beta, "' is available in '", which_model,
      "'. Check models(", substitute(obj),
      ") to see list of tests that have been run or run wald_test().")
  }

  obj$fits[[which_model]]$wald[[which_beta]]
}

# TODO: get betas
# TODO: get covars from betas
