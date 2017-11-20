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
  cat("transcripts: ", length(obj$abundance$target_id), "\n")
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
  stopifnot(all(target_ids %in% obj$abundance$target_id))

  subset_num <- length(unique(target_ids))
  new_obj <- obj
  new_obj$abundance <- new_obj$abundance[which(new_obj$abundance %in% target_ids), ]
  new_obj$bootstrap <- lapply(new_obj$bootstrap, function(bs) {
    bs[which(bs$target_id %in% target_ids), ]
  })
  attr(new_obj, "num_targets") <- subset_num
  new_obj
}
