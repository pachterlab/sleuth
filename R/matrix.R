
#' Convert a sleuth object to matrix
#'
#' Convert a sleuth object to a matrix with the condition names.
#'
#' @param obj a \code{sleuth} object
#' @param which_df character vector of length one. Which type of data to use
#' ("obs_norm" or "obs_raw")
#' @param which_units character vector of length one. Which units to use ("tpm"
#' or "est_counts")
#' @export
sleuth_to_matrix <- function(obj, which_df, which_units) {
  if ( !(which_df %in% c("obs_norm", "obs_raw")) ) {
    stop("Invalid object")
  }
  if ( !(which_units %in% c("tpm", "est_counts")) ) {
    stop("Invalid units")
  }

  data <- as.data.frame(obj[[which_df]])

  res <- list()

  s_data <- data %>%
    select_("target_id", "sample", which_units) %>%
    tidyr::spread_("sample", which_units)
  rownames(s_data) <- s_data$target_id
  s_data$target_id <- NULL
  s_data <- as.matrix(s_data)
  s_data <- s_data[, sample(1:ncol(s_data))]
  res[["data"]] <- s_data

  condition_order <- match(colnames(s_data),
    as.character(obj$sample_to_condition$sample))
  res[["condition"]] <- obj$sample_to_condition$condition[condition_order]

  res
}
