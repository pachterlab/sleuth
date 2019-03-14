#
#    sleuth: inspect your RNA-Seq with a pack of kallistos
#
#    Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @export
print.kallisto <- function(obj) {
  cat("\tkallisto object\n")
  cat("\n")
  cat("transcripts: ", attr(obj, "num_targets"), "\n")
  cat("original number of transcripts: ", attr(obj, "original_num_targets"), "\n")

  subset <- ifelse(is_kallisto_subset(obj), "Original", "Subset")
  cat("Original or Subset: ", subset)

  n_bootstrap <- length(obj$bootstrap)

  cat("bootstraps: ", n_bootstrap, "\n")

  invisible(obj)
}

#' @export
head.kallisto <- function(obj) {
  print(obj)
}

# ' extract the bias table
# '
# ' Extract bias table from either a sleuth or kallisto object
# ' @param obj either a sleuth or kallisto object
# ' @return a \code{data.frame} of bias weights
#' @export
bias_table <- function(obj, ...) {
  UseMethod('bias_table')
}

#' @export
bias_table.sleuth <- function(obj, sample) {
  stopifnot( length(sample) == 1 )

  if ( is(sample, 'numeric') || is(sample, 'integer') ) {
    sample <- as.integer(sample)
  } else {
    sample <- which( obj$sample_to_covariates$sample == sample )
    if ( length(sample) == 0 ) {
      stop('Could not find: "', sample, '"')
    }
  }

  bias_table(obj$kal[[sample]])
}

#' @export
bias_table.kallisto <- function(obj) {
  if ( length(obj$fld) == 1 && all(is.na(obj$fld)) ) {
    stop("kallisto object does not contain the fragment length distribution.",
      "Please rerun with a new version of kallisto.")
  }

  adf(
    hexamer = hexamers,
    expected_counts = obj$bias_normalized,
    observed_counts = obj$bias_observed,
    bias_weights = obj$bias_observed / obj$bias_normalized
    )
}

#' Subset a kallisto object
#'
#' Use a vector of target_ids to subset a single kallisto
#' object or all of kallisto objects in a sleuth object.
#'
#' @param obj either a sleuth or kallisto object
#' @return a new object with only the subset of target_ids included
#' @export
subset_kallisto <- function(obj, ...) {
  UseMethod('subset_kallisto')
}

#' @export
subset_kallisto.sleuth <- function(obj, target_ids) {
  stopifnot(is(target_ids, "character"))
  stopifnot(all(target_ids %in% names(obj$filter_bool)))

  if(length(target_ids) != length(unique(target_ids))) {
    target_tab <- table(target_ids)
    stop("There is at least one target_id that is duplicated. ",
         "Please provide a vector with unique target_ids.\n",
         "Here are the target_ids that are duplicated: ",
         paste(names(target_tab)[target_tab > 1], collapse = ", "))
  }

  obj$kal <- lapply(obj$kal, function(kal) {
    subset_kallisto(kal)
  })

  obj
}

#' @export
subset_kallisto.kallisto <- function(obj, target_ids) {
  stopifnot(is(target_ids, "character"))
  stopifnot(all(target_ids %in% obj$abundance$target_id))

  subset_num <- length(unique(target_ids))
  new_obj <- obj
  new_obj$abundance <- new_obj$abundance[which(new_obj$abundance$target_id %in% target_ids), ]
  new_obj$bootstrap <- lapply(new_obj$bootstrap, function(bs) {
    bs[which(bs$target_id %in% target_ids), ]
  })
  attr(new_obj, "num_targets") <- subset_num
  excluded_ids <- obj$abundance$target_id[which(!(obj$abundance$target_id %in% target_ids))]
  if(length(new_obj$excluded_ids) == 0) {
    new_obj$excluded_ids <- excluded_ids
  } else {
    new_obj$excluded_ids <- c(new_obj$excluded_ids, excluded_ids)
  }

  new_obj
}

#' Is Kallisto Object Subsetted?
#'
#' This function returns a boolean for whether the kallisto object has been
#' subsetted or if it contains the original set of target_ids analyzed by kallisto.
#' If a sleuth object is given, then a vector of booleans for each kallisto object
#' is returned.
#'
#' @param obj a sleuth or kallisto object
#' @return a boolean vector with \code{TRUE} if the kallisto object has been subsetted
#'   (i.e. not all originally analyzed target_ids are represented), or \code{FALSE} if
#'   the kallisto object contains the full original set of target_ids. For a sleuth object,
#'   a boolean vector for all of the contained kallisto objects is returned.
#' @export
is_kallisto_subset <- function(obj) {
  UseMethod("is_kallisto_subset")
}

#' @export
is_kallisto_subset.sleuth <- function(obj) {
  sapply(obj$kal, function(kal) is_kallisto_subset(kal))
}

#' @export
is_kallisto_subset.kallisto <- function(obj) {
  # If the num_targets is not equal to the original num_targets
  # then the kallisto object is a subsetted object
  attr(obj, "num_targets") != attr(obj, "original_num_targets")
}

#' Excluded IDs in Kallisto object
#'
#' Returns the excluded IDs if a kallisto object has been subsetted.
#' Will return an empty character vector if it is the original object.
#'
#' @param obj a kallisto object
#' @return a character vector containing all of the excluded IDs. This will be
#'   an empty character vector if the kallisto is the original object.
#' @export
excluded_ids <- function(obj) {
  stopifnot(is(obj, "kallisto"))

  obj$excluded_ids
}
