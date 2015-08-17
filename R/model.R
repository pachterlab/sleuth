#' Print sleuth model
#'
#' Print a model that has been fit by sleuth
#'
#' @param obj a \code{sleuth_model} object
#' @return obj (invisible)
#' @export
print.sleuth_model <- function(obj) {
  cat('formula: ', deparse(obj$formula), '\n')
  cat('coefficients:\n')
  cat(paste0('\t', colnames(obj$design_matrix), '\n'))
  if (!is.null(obj$wald)) {
    cat('tests:\n')
    cat(paste0('\t', names(obj$wald)), '\n')
  } else {
    cat('no tests found.\n')
  }

  invisible(obj)
}

#' View which models have been fit
#'
#' View which models have been fit
#'
#' @param obj
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

  invisible(obj$fits)
}

#' @export
models.sleuth_model <- function(obj) {
  print(obj)
}


#' @export
tests <- function(obj) {
  UseMethod('tests')
}

#' @export
tests.sleuth_model <- function(obj) {
  names(obj$wald)
}

#' Extract Wald test results from a sleuth object
#'
#' This function extracts Wald test results from a sleuth object.
#'
#' @param obj a \code{sleuth} object
#' @param which_beta a character string denoting which coefficient test to
#' extract
#' @param which_model a character string denoting which model to extract
#' @param rename_cols if \code{TRUE} will rename some columns to be shorter and
#' consistent with vignette
#' @param show_all if \code{TRUE} will show all transcripts (not only the ones
#' passing filters). The transcripts that do not pass filters will have
#' \code{NA} values in most columns.
#' @return a \code{data.frame}
#' @seealso \code{\link{wald_test}} to compute tests, \code{\link{models}} to
#' view which models and betas have been tested
#' @export
sleuth_results <- function(obj, which_beta, which_model = 'full', rename_cols = TRUE, show_all = TRUE) {
  stopifnot( is(obj, 'sleuth') )

  if ( !model_exists(obj, which_model) ) {
    stop("'", which_model, "' does not exist in ", substitute(obj),
      ". Please check  models(", substitute(obj), ") for fitted models.")
  }

  if ( !is(which_beta, 'character') ) {
    stop("'", substitute(which_beta), "' is not a valid character.")
  }

  if ( length(which_beta) != 1) {
    stop("'", substitute(which_beta),
      "' is not a valid length. which_beta must be of length one.")
  }

  if ( !(which_beta %in% names(obj$fits[[which_model]]$wald)) ) {
    stop("'", which_beta, "' is not available in '", which_model,
      "'. Check models(", substitute(obj),
      ") to see list of tests that have been run or run wald_test().")
  }

  # obj$fits[[which_model]]$wald[[which_beta]]
  res <- NULL
  if (rename_cols) {
    res <- dplyr::select(obj$fits[[which_model]]$wald[[which_beta]],
      target_id,
      mean_obs,
      var_obs,
      tech_var = sigma_q_sq,
      sigma_sq,
      smooth_sigma_sq,
      final_sigma_sq = smooth_sigma_sq_pmax,
      b,
      se_b,
      pval,
      qval
      )
  } else {
    res <- dplyr::select(obj$fits[[which_model]]$wald[[which_beta]],
      target_id,
      mean_obs,
      var_obs,
      sigma_q_sq,
      sigma_sq,
      smooth_sigma_sq,
      smooth_sigma_sq_pmax,
      b,
      se_b,
      pval,
      qval
      )
  }

  if ( !is.null(obj$target_mapping) ) {
    res <- dplyr::left_join(
      data.table::as.data.table(res),
      data.table::as.data.table(obj$target_mapping),
      by = 'target_id')
  }

  if (show_all) {
    tids <- adf(target_id = obj$kal[[1]]$abundance$target_id)
    res <- dplyr::left_join(
      data.table::as.data.table(tids),
      data.table::as.data.table(res),
      by = 'target_id'
      )
  }
  res <- as_df(res)

  res
}

# TODO: get betas
# TODO: get covars from betas
