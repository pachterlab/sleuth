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
