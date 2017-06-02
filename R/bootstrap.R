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

# Bootstrap to matrix
#
# Takes a \code{kallisto} object and converts the bootstrap results into
# a proper \code{matrix}
#
# @param kal a kallisto object with non-null member \code{bootstrap}
# @return a matrix with rownames equal to target_id
bootstrap2mat <- function(kal, column = "tpm") {
    stopifnot(is(kal, "kallisto"))
    # TODO: check if "column" is a valid kallisto column

    # assumes that all bootstrap samples are in same order (from read_kallisto)

    all_boot <- kal$bootstrap
    mat <- matrix(unlist(lapply(all_boot, function(bs) {
            bs[column]
        })), nrow = nrow(all_boot[[1]]))

    rownames(mat) <- all_boot[[1]]$target_id

    mat
}

#' Extract bootstrap for a specific transcript
#'
#' Extract bootstrap for a specific transcript.
#' CURRENTLY DEPRECATED: Will probably be re implemented in the next version.
#' Currently not working because of a complete rewrite of the bootstrap code.
#' If you are interested in getting the bootstraps, you can manually write some code
#' using `read_kallisto()`. Please make sure to comment in the user group if you are using this function.
#'
#' @param obj an object
#' @param ... arguments passed to other functions
#' @return a \code{data.frame} with bootstrap samples
#' @export
get_bootstraps <- function(obj, ...) {
  UseMethod('get_bootstraps')
}

#' @export
get_bootstraps.sleuth <- function(obj, transcript, max_bs = 30) {
  res <- lapply(seq_along(obj$kal),
    function(i) {
      cur <- get_bootstraps(obj$kal[[i]], transcript, max_bs)
      if (nrow(cur) == 0) {
        return(cur)
      }
      samp <- obj$sample_to_covariates$sample[i]
      cur$sample <- samp
      cur
    })
  res <- dplyr::bind_rows(res)
  if (nrow(res) > 0) {
    return(
      dplyr::left_join(res, obj$sample_to_covariates, by = c('sample'))
      )
  }

  res
}

#' @export
get_bootstraps.kallisto <- function(kal, transcript, max_bs = 30) {
  max_bs <- min(max_bs, length(kal$bootstrap))
  t_idx <- which(kal$bootstrap[[1]]$target_id == transcript)
  if (length(t_idx) == 0) {
    return( data.frame() )
  }

  bs <- lapply(1:max_bs,
    function(i) {
      dplyr::slice(kal$bootstrap[[i]], t_idx)
    })

  dplyr::bind_rows(bs)
}


# Convert kallisto bootstraps into a molten data.frame
#
# Melt it!
#
# @param kal a kallisto object
# @param column the column to pull out of the kallisto results (default = "tpm")
# @return a molten data.frame with columns "target_id", "sample" and the selected variable
# @export
melt_bootstrap <- function(kal, column = "tpm", transform = identity) {
    stopifnot(is(kal, "kallisto"))
  stopifnot(length(kal$bootstrap) > 0)

    all_boot <- kal$bootstrap
    boot <- data.frame(lapply(all_boot, select_, .dots = list(column)))
    boot <- transform(boot)
    bs_names <- paste0("bs", 1:ncol(boot))
    data.table::setnames(boot, colnames(boot), bs_names)
    boot <- boot %>%
        mutate(target_id = all_boot[[1]]$target_id)

    tidyr::gather_(boot, "sample", column, bs_names) %>%
      mutate(sample = as.factor(sample))
}

# Aggregate bootstrap samples
#
# A faster way to aggregate bootstrap samples based off of some mapping.
#
# @param kal a \code{kallisto} object
# @param mapping a data.frame containing the mapping with columns
# \code{target_id} and a column specified in 'split_by'
# @param split_by a character string of length one denoting the column in \code{mapping} to split by (such as \code{gene_id})
# @param aggregate_fun a function to aggregate
# @return a data.frame nrow(mapping) rows that has been aggregated
# groupwise using \code{aggregate_fun}
# @export
aggregate_bootstrap <- function(kal, mapping, split_by = "gene_id",
  column = "tpm", aggregate_fun = sum) {

  stopifnot( is(kal, "kallisto") )

  if ( !(column %in% c("tpm", "est_counts")) ) {
    stop("Unit must be 'tpm' or 'est_counts'")
  }

  if ( !("target_id" %in% colnames(mapping)) ) {
    stop("Column 'target_id' not found in 'mapping'")
  }

  if ( !(split_by %in% colnames(mapping)) ) {
    stop("Column 'mapping' not found in 'mapping'")
  }

  if ( any(!complete.cases(mapping)) ) {
    warning("Found some NAs in mapping. Removing them.")
    mapping <- mapping[complete.cases(mapping), ]
  }

  m_bs <- melt_bootstrap(kal, column)

  m_bs <- inner_join(
    data.table::data.table(m_bs, key = "target_id"),
    data.table::data.table(mapping, key = "target_id"),
    by = c("target_id")
    )

  m_bs <- m_bs %>%
    dplyr::group_by_(split_by, "sample") %>%
    dplyr::summarise_(
      value = lazyeval::interp( ~aggregate_fun(x), x = as.name(column) )
      )
  data.table::setnames(m_bs, "value", column)

  as.data.frame(m_bs)
}

# Summarize bootstrap values
#
# Compute the mean, sd, var, and coefficient of variation from a kallisto
# bootstrap
# @param kal a kallisto object with a non-null bootstrap list
# @param column the column to select (rho, tpm, est_counts
# @return a summarized data.frame
# @export
summarize_bootstrap <- function(kal, column = "tpm", transform = identity) {
    stopifnot(is(kal, "kallisto"))
    bs <- melt_bootstrap(kal, column, transform)

    mean_col <- paste0("bs_mean_", column)
    sd_col <- paste0("bs_sd_", column)
    var_col <- paste0("bs_var_", column)
    cv_col <- paste0("bs_cv_", column)

    bs <- bs %>%
        group_by(target_id) %>%
        summarise_(.dots = setNames(list(
                    interp(quote(mean(x)), x = as.name(column)),
                    interp(quote(sd(x)), x = as.name(column)),
                    interp(quote(var(x)), x = as.name(column))
                    ),
                c(mean_col, sd_col, var_col)))

    bs <- bs %>%
        mutate_(.dots = setNames(list(
                    interp(quote(x / y),
                        x = as.name(sd_col), y = as.name(mean_col))),
                c(cv_col)
                ))
    bs
}

# Normalize bootstrap samples
#
# Normalize by dividing by the "size factor" for each TPM and estimated counts
#
# @param kal a kallisto object
# @param tpm_size_factor the size factor (numeric length 1)
# @param est_counts_size_factor the size factor (numeric length 1)
# @export
normalize_bootstrap <- function(kal, tpm_size_factor, est_counts_size_factor) {
  stopifnot(is(kal, "kallisto"))

  calc_norm_tpm <- !missing(tpm_size_factor)
  calc_norm_counts <- !missing(est_counts_size_factor)

  if (calc_norm_tpm) {
    stopifnot(length(tpm_size_factor) == 1)
  }

  if (calc_norm_counts) {
    stopifnot(length(est_counts_size_factor) == 1)
  }

  bs <- lapply(kal$bootstrap, function(bs_tbl) {
      if (calc_norm_tpm)
        bs_tbl$tpm <- bs_tbl$tpm / tpm_size_factor
      if (calc_norm_counts)
        bs_tbl$est_counts <- bs_tbl$est_counts / est_counts_size_factor

      bs_tbl
    })
  kal$bootstrap <- bs

  kal
}

#' bootstrap summary
#'
#' Extract the bootstrap summary from a sleuth object that has been initialized in sleuth_prep.
#'
#' @param obj a \code{sleuth} object such that \code{extra_bootstrap_summary = TRUE} inside of \code{\link{sleuth_prep}}.
#' @param target_id a character vector of length 1 indicating the target_id (transcript or gene name depending on aggregation mode)
#' @param units a character vector of either 'est_counts' or 'tpm' (also requires \code{extra_bootstrap_summary = TRUE} in \code{\link{sleuth_prep}})
#' @return a \code{data.frame} with the summary statistics across all samples for that particular target
#' @export
get_bootstrap_summary <- function(obj, target_id, units = 'est_counts') {
  stopifnot( is(obj, 'sleuth') )

  if (units != 'est_counts' && units != 'tpm' && units != 'scaled_reads_per_base') {
    stop(paste0("'", units, "' is invalid for 'units'. please see documentation"))
  }

  if (is.null(obj$bs_quants)) {
    if (units == 'est_counts') {
      stop("bootstrap summary missing. rerun sleuth_prep() with argument 'extra_bootstrap_summary = TRUE'")
    } else {
      stop("bootstrap summary missing. rerun sleuth_prep() with argument 'extra_bootstrap_summary = TRUE' and 'read_bootstrap_tpm = TRUE'")
    }
  }

  if (!(target_id %in% rownames(obj$bs_quants[[1]][[units]]))) {
    stop(paste0("couldn't find target_id '", target_id, "'"))
  }

  df <- as_df(
    do.call(rbind,
      lapply(obj$bs_quants,
      function(sample_bs) {
        sample_bs[[units]][target_id, ]
      })
      )
    )
  df <- dplyr::bind_cols(df, obj$sample_to_covariates)

  df
}


# Sample bootstraps
#
# From a sleuth object, create experiments by randomly sampling bootstraps from each kallisto object
# @param obj a \code{kallisto} object
# @param n_samples the number of samples to genenerate
# @export
sample_bootstrap <- function(obj, n_samples = 100L) {
  stopifnot( is(obj, "sleuth") )

  n_kal <- length(obj$kal)
  n_bs_per_samp <- unlist(lapply(obj$kal, function(x) length(x$bootstrap)))
  if (any(n_bs_per_samp < n_samples)) {
    warning("You've asked to sample more samples than you have bootstraps.",
      " We recommend you generate more bootstrap samples in kallisto...")
  }

  which_samp <- lapply(seq_along(n_bs_per_samp),
    function(i) {
      cur_n <- n_bs_per_samp[i]
      sample.int(cur_n, n_samples, replace = TRUE)
    })
  # each column contains which bootstrap sample we want from each kallisto
  which_samp <- t(simplify2array(which_samp))

  # allocate the matrices
  sample_mat <- lapply(1:n_samples,
    function(discard) {
      mat <- matrix(NA_real_, nrow = nrow(obj$kal[[1]]$abundance),
        ncol = nrow(which_samp))
      rownames(mat) <- obj$kal[[1]]$abundance$target_id
      colnames(mat) <- obj$sample_to_condition$sample
      mat
    })

  # go through the generated sample ids and place them into the corresponding
  # matrix sample
  for (s in 1:n_samples) {
    for (idx in 1:nrow(which_samp)) {
      b <- which_samp[idx, s]
      sample_mat[[s]][, idx] <- obj$kal[[idx]]$bootstrap[[b]]$est_counts
    }
  }

  sample_mat
}

# @export
dcast_bootstrap <- function(obj, ...) {
  UseMethod("dcast_bootstrap")
}

# @export
dcast_bootstrap.sleuth <- function(obj, units, nsamples = NULL) {
  bs <- lapply(obj[["kal"]], dcast_bootstrap, units, nsamples)

  do.call(cbind, bs)
}

# @export
dcast_bootstrap.kallisto <- function(obj, units, nsamples = NULL) {
  if ( !(units %in% c("est_counts", "tpm")) ) {
    stop(paste0("'", substitute(units),
        "' is not valid for 'units'. Please see documentation"))
  }

  if ( length(obj$bootstrap) < 1 ) {
    stop("No bootstrap samples found.")
  }

  n_bs <- length(obj$bootstrap)
  n_features <- nrow(obj$bootstrap[[1]])

  which_bs <- seq_along(obj$bootstrap)
  if (!is.null(nsamples) && nsamples < n_bs) {
    print(nsamples)
    which_bs <- sample.int(n_bs, nsamples)
  }

  mat <- matrix(NA_real_, nrow = n_features, ncol = length(which_bs))

  for (j in seq_along(which_bs)) {
    mat[, j] <- obj[[ "bootstrap" ]][[which_bs[j]]][[ units ]]
  }
  rownames(mat) <- obj[["bootstrap"]][[1]][["target_id"]]

  mat
}

# Function to process bootstraps for parallelization
process_bootstrap <- function(i, samp_name, kal_path,
                              num_transcripts, est_count_sf,
                              read_bootstrap_tpm, gene_mode,
                              extra_bootstrap_summary,
                              target_id, mappings, which_ids,
                              aggregation_column, transform_fun)
{
  dot(i)
  bs_quants <- list()

  num_bootstrap <- as.integer(rhdf5::h5read(kal_path$path,
                                            "aux/num_bootstrap"))
  if (num_bootstrap == 0) {
    stop(paste0("File ", kal_path, " has no bootstraps.",
                "Please generate bootstraps using \"kallisto quant -b\"."))
  }

  # TODO: only perform operations on filtered transcripts
  eff_len <- rhdf5::h5read(kal_path$path, "aux/eff_lengths")
  bs_mat <- read_bootstrap_mat(fname = kal_path$path,
                               num_bootstraps = num_bootstrap,
                               num_transcripts = num_transcripts,
                               est_count_sf = est_count_sf)

  if (read_bootstrap_tpm) {
    bs_quant_tpm <- aperm(apply(bs_mat, 1, counts_to_tpm,
                                eff_len))
    colnames(bs_quant_tpm) <- colnames(bs_mat)

    # gene level code is analogous here to below code
    if (gene_mode) {
      colnames(bs_quant_tpm) <- target_id
      # Make bootstrap_num an explicit column; each is treated as a "sample"
      bs_tpm_df <- data.frame(bootstrap_num = c(1:num_bootstrap),
                              bs_quant_tpm, check.names = F)
      rm(bs_quant_tpm)
      # Make long tidy table; this step is much faster
      # using data.table melt rather than tidyr gather
      tidy_tpm <- data.table::melt(bs_tpm_df, id.vars = "bootstrap_num",
                                   variable.name = "target_id",
                                   value.name = "tpm")
      tidy_tpm <- data.table::as.data.table(tidy_tpm)
      rm(bs_tpm_df)
      tidy_tpm$target_id <- as.character(tidy_tpm$target_id)
      tidy_tpm <- merge(tidy_tpm, mappings, by = "target_id",
                        all.x = T)
      # Data.table dcast uses non-standard evaluation
      # So quote the full casting formula to make sure
      # "aggregation_column" is interpreted as a variable
      # see: http://stackoverflow.com/a/31295592
      quant_tpm_formula <- paste("bootstrap_num ~",
                                 aggregation_column)
      bs_quant_tpm <- data.table::dcast(tidy_tpm,
                                        quant_tpm_formula, value.var = "tpm",
                                        fun.aggregate = sum)
      bs_quant_tpm <- as.matrix(bs_quant_tpm[, -1])
      rm(tidy_tpm) # these tables are very large
    }
    bs_quant_tpm <- aperm(apply(bs_quant_tpm, 2,
                                quantile))
    colnames(bs_quant_tpm) <- c("min", "lower", "mid",
                                "upper", "max")
    bs_quants$tpm <- bs_quant_tpm
  }

  if (gene_mode) {
    # I can combine target_id and eff_len
    # I assume the order is the same, since it's read from the same kallisto
    # file and each kallisto file has the same order
    eff_len_df <- data.frame(target_id, eff_len,
                             stringsAsFactors = F)
    # make bootstrap number an explicit column to facilitate melting
    bs_df <- data.frame(bootstrap_num = c(1:num_bootstrap),
                        bs_mat, check.names = F)
    rm(bs_mat)
    # data.table melt function is much faster than tidyr's gather function
    # output is a long table with each bootstrap's value for each target_id
    tidy_bs <- data.table::melt(bs_df, id.vars = "bootstrap_num",
                                variable.name = "target_id",
                                value.name = "est_counts")
    rm(bs_df)
    # not sure why, but the melt function always returns a factor,
    # even when setting variable.factor = F, so I coerce target_id
    tidy_bs$target_id <- as.character(tidy_bs$target_id)
    # combine the long tidy table with eff_len and aggregation mappings
    # note that bootstrap number is treated as "sample" here
    # for backwards compatibility
    tidy_bs <- dplyr::select(tidy_bs, target_id,
                             est_counts, sample = bootstrap_num)
    tidy_bs <- merge(data.table::as.data.table(tidy_bs),
                     data.table::as.data.table(eff_len_df), by = "target_id",
                     all.x = TRUE)
    tidy_bs <- merge(tidy_bs, mappings, by = "target_id",
                     all.x = TRUE)
    # create the median effective length scaling factor for each gene
    scale_factor <- tidy_bs[, scale_factor := median(eff_len),
                            by = eval(parse(text=aggregation_column))]
    # use the old reads_per_base_transform method to get gene scaled counts
    scaled_bs <- reads_per_base_transform(tidy_bs,
                                          scale_factor$scale_factor,
                                          aggregation_column,
                                          mappings)
    # this step undoes the tidying to get back a matrix format
    # target_ids here are now the aggregation column ids
    bs_mat <- data.table::dcast(scaled_bs, sample ~ target_id,
                                value.var = "scaled_reads_per_base")
    # this now has the same format as the transcript matrix
    # but it uses gene ids
    bs_mat <- as.matrix(bs_mat[, -1])
    rm(tidy_bs, scaled_bs)
  }

  if (extra_bootstrap_summary) {
    bs_quant_est_counts <- aperm(apply(bs_mat, 2,
                                       quantile))
    colnames(bs_quant_est_counts) <- c("min", "lower",
                                       "mid", "upper", "max")
    bs_quants$est_counts <- bs_quant_est_counts
  }

  bs_mat <- transform_fun(bs_mat)
  # If bs_mat was made at gene-level, already has column names
  # If at transcript-level, need to add target_ids
  if(!gene_mode) {
    colnames(bs_mat) <- target_id
  } else {
    # rename est_counts to scaled_reads_per_base
    bs_quants$scaled_reads_per_base <- bs_quants$est_counts
    bs_quants$est_counts <- NULL
  }
  # all_sample_bootstrap[, i] bootstrap point estimate of the inferential
  # variability in sample i
  # NOTE: we are only keeping the ones that pass the filter
  bootstrap_result <- matrixStats::colVars(bs_mat[, which_ids])

  list(index = i, bs_quants = bs_quants, bootstrap_result = bootstrap_result)
}
