#' kallisto coefficient of variation
#'
#' cv
#'
#' @param kal a kallisto object
#' @return a data.frame with the coefficient of variation
coef_var <- function(kal) {
    stopifnot(is(kal, "kallisto"))

}

#' @export
print.kallisto <- function(obj) {
  cat("\tkallisto object\n")
  cat("\n")
  cat("transcripts: ", length(obj$abundance$target_id), "\n")
  n_bootstrap <- length(obj$bootstrap)

  cat("bootstraps: ", n_bootstrap, "\n")

  invisible(obj)
}

#' @export
head.kallisto <- function(obj) {
  print(obj)
}
