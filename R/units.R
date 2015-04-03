#' @export
counts_to_tpm <- function(est_counts, eff_len) {
  stopifnot( length(eff_len) == length(est_counts) )

  which_valid <- which(eff_len > 0)

  num <- (est_counts / eff_len)
  num[-which_valid] <- 0
  denom <- sum(num)

  (1e6 * num) / denom
}

#' @export
tpm_to_alpha <- function(tpm, eff_len) {
  stopifnot( length(tpm) == length(eff_len) )

  num <- tpm * eff_len
  denom <- sum(num)

  num / denom
}
