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

#' Basic row filter
#'
#' A basic filter to be used.
#'
#' @param row this is a vector of numerics that will be passedin
#' @param min_reads the minimum mean number of reads
#' @param min_prop the minimum proportion of reads to pass this filter
#' @return a logical of length 1
#' @export
basic_filter <- function(row, min_reads = 5, min_prop = 0.47) {
  mean(row >= min_reads) >= min_prop
}

#' Natural log and offset transformation
#'
#' The default transformation function for converting the normalized counts.
#'
#' @param x numeric that must be >=0. represents an individual observed count (one transcript in one sample).
#' @param offset numeric offset to prevent taking the log of 0.
#' @return log(x + offset)
#' @export
log_transform <- function(x, offset=0.5) {
  log(x + offset)
}

# currently defunct
filter_df_all_groups <- function(df, fun, group_df, ...) {
  grps <- setdiff(colnames(group_df), 'sample')
  res <- sapply(unique(grps),
    function(g) {
      apply(filter_df_by_groups(df, fun, group_df[, c('sample', g)], ...), 1, any)
    })

  vals <- apply(res, 1, any)
  names(vals) <- rownames(df)

  vals
}

# currently defunct
filter_df_by_groups <- function(df, fun, group_df, ...) {
  stopifnot(ncol(group_df) == 2)
  grps <- as.character(group_df[[setdiff(colnames(group_df), 'sample')]])
  sapply(unique(grps),
    function(g) {
      valid_samps <- grps %in% g
      valid_samps <- group_df$sample[valid_samps]
      apply(df[, valid_samps], 1, fun, ...)
    })
}

#' Constructor for a 'sleuth' object
#'
#' A sleuth is a group of kallistos. Borrowing this terminology, a 'sleuth' object stores
#' a group of kallisto results, and can then operate on them while
#' accounting for covariates, sequencing depth, technical and biological
#' variance.
#'
#' @param sample_to_covariates a \code{data.frame} which contains a mapping
#' from \code{sample} (a column) to some set of experimental conditions or
#' covariates. The column \code{path} is also required, which is a character
#' vector where each element points to the corresponding kallisto output directory. The column
#' \code{sample} should be in the same order as the corresponding entry in
#' \code{path}.
#' @param full_model an R \code{formula} which explains the full model (design)
#' of the experiment OR a design matrix. It must be consistent with the data.frame supplied in
#' \code{sample_to_covariates}. You can fit multiple covariates by joining them with '+' (see example)
#' @param filter_fun the function to use when filtering.
#' @param target_mapping a \code{data.frame} that has at least one column
#' 'target_id' and others that denote the mapping for each target. if it is not
#' \code{NULL}, \code{target_mapping} is joined with many outputs where it
#' might be useful. For example, you might have columns 'target_id',
#' 'ensembl_gene' and 'entrez_gene' to denote different transcript to gene
#' mappings.
#' @param max_bootstrap maximum number of bootstrap values to read for each
#' transcript.
#' @param ... additional arguments passed to the filter function
#' @param norm_fun_counts a function to perform between sample normalization on the estimated counts.
#' @param norm_fun_tpm a function to perform between sample normalization on the TPM
#' @param aggregation_column a string of the column name in \code{\link{target_mapping}} to aggregate targets
#' @param read_bootstrap_tpm read and compute summary statistics on bootstraps on the TPM.
#' NOTE: Unnecessary for typical analyses
#' @param extra_bootstrap_summary if \code{TRUE}, compute extra summary
#' statistics needed for some plots (e.g. \code{\link{plot_bootstrap}}).
#' NOTE: Unnecessary for typical analyses
#' @param transformation_function the transformation that should be applied
#' to the normalized counts. Default is \code{'log(x+0.5)'} (i.e. natural log with 0.5 offset)
#' NOTE: be sure you know what you're doing before you change this.
#' @param num_cores an integer of the number of computer cores mclapply should use
#' to speed up sleuth preparation
#' @return a \code{sleuth} object containing all kallisto samples, metadata,
#' and summary statistics
#' @examples # Assume we have run kallisto on a set of samples, and have two treatments,
#' genotype and drug.
#' colnames(s2c)
#' # [1] "sample"  "genotype"  "drug"  "path"
#' so <- sleuth_prep(s2c, ~genotype + drug)
#' @seealso \code{\link{sleuth_fit}} to fit a model, \code{\link{sleuth_wt}} or
#' \code{\link{sleuth_lrt}} to perform hypothesis testing
#' @export
sleuth_prep <- function(
  sample_to_covariates,
  full_model,
  filter_fun = basic_filter,
  target_mapping = NULL,
  max_bootstrap = NULL,
  norm_fun_counts = norm_factors,
  norm_fun_tpm = norm_factors,
  aggregation_column = NULL,
  read_bootstrap_tpm = FALSE,
  extra_bootstrap_summary = FALSE,
  transformation_function = log_transform,
  num_cores = max(1L, parallel::detectCores() - 1L),
  ...) {

  ##############################
  # check inputs

  # data types

  if (!is(sample_to_covariates, "data.frame")) {
    stop(paste0("'", substitute(sample_to_covariates), "' (sample_to_covariates) must be a data.frame"))
  }

  if (!is(full_model, "formula") && !is(full_model, "matrix")) {
    stop(paste0("'", substitute(full_model), "' (full_model) must be a formula or a matrix"))
  }

  if (!("sample" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates),
        "' (sample_to_covariates) must contain a column named 'sample'"))
  }

  if (!("path" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates)),
      "' (sample_to_covariates) must contain a column named 'path'")
  }

  if (!is.null(target_mapping) && !is(target_mapping, 'data.frame')) {
    stop(paste0("'", substitute(target_mapping),
        "' (target_mapping) must be a data.frame or NULL"))
  } else if (is(target_mapping, 'data.frame')){
    if (!("target_id" %in% colnames(target_mapping))) {
      stop(paste0("'", substitute(target_mapping),
          "' (target_mapping) must contain a column named 'target_id'"))
    }
  }

  if (!is.null(max_bootstrap) && max_bootstrap <= 0 ) {
    stop("max_bootstrap must be > 0")
  }

  if (any(is.na(sample_to_covariates))) {
    warning("Your 'sample_to_covariance' data.frame contains NA values. This will likely cause issues later.")
  }

  if (is(full_model, "matrix") &&
      nrow(full_model) != nrow(sample_to_covariates)) {
    stop("The design matrix number of rows are not equal to the number of rows in the sample_to_covariates argument.")
  }

  if (!is(norm_fun_counts, 'function')) {
    stop("norm_fun_counts must be a function")
  }

  if (!is(norm_fun_tpm, 'function')) {
    stop("norm_fun_tpm must be a function")
  }

  if (!is.null(aggregation_column) && is.null(target_mapping)) {
    stop(paste("You provided a 'aggregation_column' to aggregate by,",
               "but not a 'target_mapping'. Please provided a 'target_mapping'."))
  }

  num_cores <- check_num_cores(num_cores)

  # TODO: ensure transcripts are in same order -- if not, report warning that
  # kallisto index might be incorrect

  # done
  ##############################

  msg('reading in kallisto results')
  sample_to_covariates$sample <- as.character(sample_to_covariates$sample)

  kal_dirs <- sample_to_covariates$path
  sample_to_covariates$path <- NULL

  msg('dropping unused factor levels')
  samples_to_covariates <- droplevels(sample_to_covariates)

  nsamp <- 0
  # append sample column to data
  kal_list <- lapply(seq_along(kal_dirs),
    function(i) {
      nsamp <- dot(nsamp)
      path <- kal_dirs[i]
      suppressMessages({
        kal <- read_kallisto(path, read_bootstrap = FALSE,
          max_bootstrap = max_bootstrap)
        })
      kal$abundance <- dplyr::mutate(kal$abundance,
        sample = sample_to_covariates$sample[i])

      kal
    })
  msg('')

  check_result <- check_kal_pack(kal_list)
  kal_versions <- check_result$versions

  obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))

  design_matrix <- NULL
  if (is(full_model, 'formula')) {
    design_matrix <- model.matrix(full_model, sample_to_covariates)
  } else {
    if (is.null(colnames(full_model))) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    design_matrix <- full_model
  }
  rownames(design_matrix) <- sample_to_covariates$sample

  obs_raw <- dplyr::arrange(obs_raw, target_id, sample)

  ###
  # try to deal with weird ensemble names
  ###
  if (!is.null(target_mapping)) {
    tmp_names <- data.frame(target_id = kal_list[[1]]$abundance$target_id,
      stringsAsFactors = FALSE)
    target_mapping <- check_target_mapping(tmp_names, target_mapping)
    rm(tmp_names)
  }

  ret <- list(
      kal = kal_list,
      kal_versions = kal_versions,
      obs_raw = obs_raw,
      sample_to_covariates = sample_to_covariates,
      bootstrap_summary = NA,
      full_formula = full_model,
      design_matrix = design_matrix,
      target_mapping = target_mapping,
      gene_mode = !is.null(aggregation_column),
      gene_column = aggregation_column,
      transform_fun = transformation_function
    )

  # TODO: eventually factor this out
  normalize <- TRUE
  if (normalize ) {

    msg("normalizing est_counts")
    est_counts_spread <- spread_abundance_by(obs_raw, "est_counts",
      sample_to_covariates$sample)
    filter_bool <- apply(est_counts_spread, 1, filter_fun, ...)
    filter_true <- filter_bool[filter_bool]

    msg(paste0(sum(filter_bool), ' targets passed the filter'))
    est_counts_sf <- norm_fun_counts(est_counts_spread[filter_bool, ])

    filter_df <- adf(target_id = names(filter_true))

    est_counts_norm <- as_df(t(t(est_counts_spread) / est_counts_sf))

    est_counts_norm$target_id <- rownames(est_counts_norm)
    est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

    obs_norm <- est_counts_norm
    obs_norm$target_id <- as.character(obs_norm$target_id)
    obs_norm$sample <- as.character(obs_norm$sample)
    rm(est_counts_norm)

    # deal w/ TPM
    msg("normalizing tpm")
    tpm_spread <- spread_abundance_by(obs_raw, "tpm",
      sample_to_covariates$sample)
    tpm_sf <- norm_fun_tpm(tpm_spread[filter_bool, ])
    tpm_norm <- as_df(t(t(tpm_spread) / tpm_sf))
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
    tpm_norm$sample <- as.character(tpm_norm$sample)

    msg('merging in metadata')
    # put everyone in the same order to avoid a slow join
    obs_norm <- dplyr::arrange(obs_norm, target_id, sample)
    tpm_norm <- dplyr::arrange(tpm_norm, target_id, sample)

    stopifnot(all.equal(obs_raw$target_id, obs_norm$target_id) &&
      all.equal(obs_raw$sample, obs_norm$sample))

    suppressWarnings({
      if (!all.equal(dplyr::select(obs_norm, target_id, sample),
          dplyr::select(tpm_norm, target_id, sample))) {
        stop('Invalid column rows. In principle, can simply join. Please report error.')
      }

      # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
      #   by = c('target_id', 'sample'))
      obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, tpm))
    })

    # add in eff_len and len
    obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(obs_raw, eff_len, len))


    obs_norm <- as_df(obs_norm)
    ret$obs_norm <- obs_norm
    ret$est_counts_sf <- est_counts_sf
    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
    ret$obs_norm_filt <- dplyr::semi_join(obs_norm, filter_df, by = 'target_id')
    ret$tpm_sf <- tpm_sf

    #### This code through the for loop is a candidate for moving to another function
    path <- kal_dirs[1]
    kal_path <- get_kallisto_path(path)
    target_id <- as.character(rhdf5::h5read(kal_path$path, "aux/ids"))
    num_transcripts <- length(target_id)
    ret$bs_quants <- list()

    which_target_id <- ret$filter_df$target_id

    if (ret$gene_mode) {
      msg(paste0("aggregating by column: ", aggregation_column))
      # Get list of IDs to aggregate on (usually genes)
      # Also get the filtered list and update the "filter_df" and "filter_bool"
      # variables for the sleuth object
      target_mapping[target_mapping[[aggregation_column]] == "",
                     aggregation_column] <- NA
      agg_id <- unique(target_mapping[, aggregation_column, with = FALSE])
      agg_id <- agg_id[[1]]
      agg_id <- agg_id[!is.na(agg_id)]
      mappings <- dplyr::select_(target_mapping, "target_id", aggregation_column)
      mappings <- data.table::as.data.table(mappings)
      which_tms <- which(mappings$target_id %in% which_target_id)
      which_agg_id <- unique(mappings[which_tms, aggregation_column, with = FALSE])
      which_agg_id <- which_agg_id[[1]]
      which_agg_id <- which_agg_id[!is.na(which_agg_id)]
      filter_df <- adf(target_id = which_agg_id)
      filter_bool <- agg_id %in% which_agg_id

      msg(paste0(length(which_agg_id), " genes passed the filter"))

      # Taken from gene_summary; scale normalized observed counts to "reads/base"
      norm_by_length <- TRUE
      tmp <- data.table::as.data.table(ret$obs_raw)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = TRUE)
      scale_factor <- tmp[, scale_factor := median(eff_len),
                          by=list(sample,eval(parse(text=aggregation_column)))]
      obs_norm_gene <- reads_per_base_transform(ret$obs_norm,
          scale_factor, aggregation_column, mappings, norm_by_length)
      # New code: get gene-level TPM (simple sum of normalized transcript TPM)
      tmp <- data.table::as.data.table(tpm_norm)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = T)
      if (any(is.na(tmp[[aggregation_column]]))) {
        rows_to_remove <- is.na(tmp[[aggregation_column]])
        num_missing <- length(unique(tmp[rows_to_remove, target_id]))
        warning(num_missing, " target_ids are missing annotations for the aggregation_column: ",
                aggregation_column, ".\nThese target_ids will be dropped from the gene-level analysis.",
                "\nIf you did not expect this, check your 'target_mapping' table for missing values.")
        tmp <- tmp[!rows_to_remove]
      }
      tpm_norm_gene <- tmp[, j = list(tpm = sum(tpm)),
                           by = list(sample, eval(parse(text = aggregation_column)))]
      data.table::setnames(tpm_norm_gene, 'parse', 'target_id')
      tpm_norm_gene <- as_df(tpm_norm_gene)

      # Same steps as above to add TPM column to "obs_norm" table
      obs_norm_gene <- dplyr::arrange(obs_norm_gene, target_id, sample)
      tpm_norm_gene <- dplyr::arrange(tpm_norm_gene, target_id, sample)

      stopifnot(all.equal(dplyr::select(obs_norm_gene, target_id, sample),
            dplyr::select(tpm_norm_gene, target_id, sample)))
      suppressWarnings({
        if ( !all.equal(dplyr::select(obs_norm_gene, target_id, sample),
            dplyr::select(tpm_norm_gene, target_id, sample), check.attributes = FALSE) ) {
              stop('Invalid column rows. In principle, can simply join. Please report error.')
            }

        # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
        #   by = c('target_id', 'sample'))
        obs_norm_gene <- dplyr::bind_cols(obs_norm_gene, dplyr::select(tpm_norm_gene, tpm))
      })

      # These are the updated gene-level variables
      ret$filter_df <- adf(target_id = which_agg_id)
      ret$filter_bool <- agg_id %in% which_agg_id
      ret$obs_norm <- obs_norm_gene
      ret$obs_norm_filt <- dplyr::semi_join(obs_norm_gene, filter_df, by = 'target_id')

      rm(obs_norm, tpm_norm, obs_norm_gene, tpm_norm_gene)

      # This is the gene-level version of the matrix
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_agg_id),
                                     ncol = length(ret$kal))
      which_ids <- which_agg_id
    } else {
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_target_id),
                                     ncol = length(ret$kal))
      which_ids <- which_target_id
    }

    msg('summarizing bootstraps')
    bs_results <- parallel::mclapply(seq_along(kal_dirs), function(i) {
      samp_name <- sample_to_covariates$sample[i]
      kal_path <- get_kallisto_path(kal_dirs[i])
      process_bootstrap(i, samp_name, kal_path,
                        num_transcripts, est_counts_sf[[i]],
                        read_bootstrap_tpm, ret$gene_mode,
                        extra_bootstrap_summary,
                        target_id, mappings, which_ids, ret$gene_column,
                        ret$transform_fun)
    }, mc.cores = num_cores)

    # if mclapply results in an error (a warning is shown), then print error and stop
    if (is(bs_results[[1]], "try-error")) {
      print(attributes(bs_results[[1]])$condition)
      stop("mclapply had an error. See the above error message for more details.")
    }

    # mclapply is expected to retun the bootstraps in order; this is a sanity check of that
    indices <- sapply(bs_results, function(result) result$index)
    stopifnot(identical(indices, order(indices)))

    if(read_bootstrap_tpm | extra_bootstrap_summary) {
      ret$bs_quants <- lapply(bs_results, function(result) result$bs_quants)
      names(ret$bs_quants) <- sample_to_covariates$sample
    }

    all_sample_bootstrap <- sapply(bs_results, function(result) result$bootstrap_result)
    rownames(all_sample_bootstrap) <- which_ids

    # end summarize bootstraps
    msg('')

    sigma_q_sq <- rowMeans(all_sample_bootstrap)

    # This is the rest of the gene_summary code
    if (ret$gene_mode) {
      names(sigma_q_sq) <- which_agg_id
      obs_counts <- obs_to_matrix(ret, "scaled_reads_per_base")[which_agg_id, ]
    } else {
      names(sigma_q_sq) <- which_target_id
      obs_counts <- obs_to_matrix(ret, "est_counts")[which_target_id, ]
    }

    sigma_q_sq <- sigma_q_sq[order(names(sigma_q_sq))]
    obs_counts <- ret$transform_fun(obs_counts)
    obs_counts <- obs_counts[order(rownames(obs_counts)),]

    ret$bs_summary <- list(obs_counts = obs_counts, sigma_q_sq = sigma_q_sq)
  }

  class(ret) <- 'sleuth'

  ret
}

# check versions of kallistos and num bootstraps, etc
check_kal_pack <- function(kal_list) {

  versions <- sapply(kal_list, function(x) attr(x, 'kallisto_version'))
  u_versions <- unique(versions)
  if (length(u_versions) > 1) {
    warning('More than one version of kallisto was used: ',
      paste(u_versions, collapse = "\n"))
  }

  ntargs <- sapply(kal_list, function(x) attr(x, 'num_targets'))
  u_ntargs <- unique(ntargs)
  if (length(u_ntargs) > 1 ) {
    stop('Inconsistent number of transcripts. Please make sure you used the same index everywhere.')
  }

  list(versions = u_versions)
}

# this function is mostly to deal with annoying ENSEMBL transcript names that
# have a training .N to keep track of version number
#
# @return the target_mapping if an intersection is found. a target_mapping that
# matches \code{t_id} if no matching is found
check_target_mapping <- function(t_id, target_mapping) {
  t_id <- data.table::as.data.table(t_id)
  target_mapping <- data.table::as.data.table(target_mapping)

  tmp_join <- dplyr::inner_join(t_id, target_mapping, by = 'target_id')

  if (!nrow(tmp_join)) {
    t_id <- dplyr::mutate(t_id, target_id_old = sub('\\.[0-9]+', '', target_id))
    target_mapping_upd <- dplyr::rename(target_mapping, target_id_old = target_id)

    tmp_join <- dplyr::inner_join(t_id, target_mapping_upd,
      c('target_id_old'))

    if (!nrow(tmp_join)) {
      stop("couldn't solve nonzero intersection")
    }

    target_mapping <- dplyr::left_join(target_mapping,
      unique(dplyr::select(tmp_join, target_id_new = target_id,
        target_id = target_id_old),
        by = NULL),
      by = 'target_id')
    # if we couldn't recover, use the old one
    target_mapping <- dplyr::mutate(target_mapping,
      target_id_new = ifelse(is.na(target_id_new), target_id, target_id_new))
    target_mapping <- dplyr::select(target_mapping, -target_id)
    target_mapping <- dplyr::rename(target_mapping, target_id = target_id_new)

    warning(paste0('intersection between target_id from kallisto runs and ',
      'the target_mapping is empty. attempted to fix problem by removing .N ',
      'from target_id, then merging back into target_mapping.',
      ' please check obj$target_mapping to ensure this new mapping is correct.'))
  }

  target_mapping
}

#' Normalization factors
#'
#' More or less the DESeq size factors
#'
#' @param mat a matrix of estimated counts
#' @return a vector of length \code{ncol(mat)} with a sequencing depth factor
#' for each sample
#' @export
norm_factors <- function(mat) {
  nz <- apply(mat, 1, function(row) !any(round(row) == 0))
  mat_nz <- mat[nz, ]
  p <- ncol(mat)
  geo_means <- exp(apply(mat_nz, 1, function(row) mean(log(row))))
  s <- sweep(mat_nz, 1, geo_means, `/`)

  sf <- apply(s, 2, median)
  scaling <- exp( (-1 / p) * sum(log(sf)))

  sf * scaling
}

# Summarize many bootstrap objects
#
# Summarize all the bootstrap samples from a kallisto run. The summarized
# values are then all put into a data.frame and stored in the \code{sleuth}
# object.
#
# @param obj a \code{sleuth} object
# @param force if \code{FALSE}, then will only compute the summary if it has
# not yet been set.
# @param verbose if \code{TRUE}, print verbosely
# @return a \code{kallisto} object with member \code{bootstrap_summary}
# updated and containing a data frame.
sleuth_summarize_bootstrap <- function(obj, force = FALSE, verbose = FALSE) {
  # TODO: make this into an S3 function 'summarize_bootstrap'
  stopifnot(is(obj, "sleuth"))

  if (!is.na(obj$bootstrap_summary) && !force) {
    if (verbose) {
      cat("summarize_bootstrap.sleuth: Already computed summary -- will not recompute.")
    }

    return(obj)
  }

  tpm_bs <- sleuth_summarize_bootstrap_col(obj, "tpm")
  counts_bs <- sleuth_summarize_bootstrap_col(obj, "est_counts")

  s_bs <- inner_join(
    data.table::data.table(tpm_bs),
    data.table::data.table(data.table(counts_bs)),
    #by = c("target_id", "sample", "condition")
    by = c("target_id", "sample")
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  obj$bootstrap_summary <- s_bs

  obj
}

sleuth_summarize_bootstrap_col <- function(obj, col, transform = identity) {
  res <- lapply(seq_along(obj$kal), function(i) {
      cur_samp <- obj$sample_to_covariates$sample[i]

      dplyr::mutate(summarize_bootstrap(obj$kal[[i]], col, transform),
        sample = cur_samp)
    })
  dplyr::bind_rows(res)
}

# Spread abundance by a column
#
# Take a data.frame from a sleuth object (e.g. \code{obs_raw}) and cast it
# into a matrix where the rows are the target_ids and the columns are the
# sample ids. The values are the variable you are "spreading" on.
# @param abund the abundance \code{data.frame} from a \code{sleuth} object
# @param var a character array of length one. The variable for which to get "spread" on (e.g. "est_counts").
# @export
spread_abundance_by <- function(abund, var, which_order) {
  # var <- lazyeval::lazy(var)
  var_spread <- abund %>%
    select_("target_id", "sample", var) %>%
    tidyr::spread_("sample", var) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL

  result <- as.matrix(var_spread)

  result[, which_order, drop = FALSE]
}

#' @export
melt_bootstrap_sleuth <- function(obj) {
  # TODO: make this into a S3 function
  lapply(seq_along(obj$kal), function(i) {
      cur_samp <- obj$sample_to_condition$sample[i]
      cur_cond <- obj$sample_to_condition$condition[i]

      melt_bootstrap(obj$kal[[i]]) %>%
        rename(bs_sample = sample) %>%
        mutate(sample = cur_samp, condition = cur_cond)
    }) %>% rbind_all()
}

#' Get kallisto abundance table from a sleuth object
#'
#' Get back the kallisto abundance table from a sleuth object.
#'
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, return the filtered data
#' @param normalized if TRUE, return the normalized data. Otherwise, return
#' the raw data
#' @param include_covariates if TRUE, join the covariates
#' @return a \code{data.frame} with at least columns
#' @export
kallisto_table <- function(obj,
  use_filtered = FALSE,
  normalized = TRUE,
  include_covariates = TRUE) {

  res <- NULL
  if (use_filtered && normalized) {
    res <- obj$obs_norm_filt
  } else if (normalized) {
    res <- obj$obs_norm
  } else {
    res <- obj$obs_raw
  }

  if (use_filtered && !normalized) {
    res <- dplyr::semi_join(
      data.table::as.data.table(res),
      data.table::as.data.table(obj$filter_df),
      by = 'target_id')
  }

  if (include_covariates) {
    res <- dplyr::left_join(
      data.table::as.data.table(res),
      data.table::as.data.table(obj$sample_to_covariates),
      by = 'sample')
  }

  as_df(res)
}

# TODO: deprecate?
# observations to a matrix
#
# observations to a matrix
#
# @param obj is a sleuth object
# @param value_name either "est_counts" or "tpm"
# @return a matrix with the appropriate names
obs_to_matrix <- function(obj, value_name) {

  obs_counts <- reshape2::dcast(obj$obs_norm, target_id ~ sample,
    value.var = value_name)

  obs_counts <- as.data.frame(obs_counts)
  rownames(obs_counts) <- obs_counts$target_id
  obs_counts$target_id <- NULL
  obs_counts <- as.matrix(obs_counts)
  obs_counts <- obs_counts[, obj$sample_to_covariates$sample]

  obs_counts
}


# TODO: check if works -- currently untested

#' Get a data.frame from all kallisto objects
#'
#' Build a data.frame from all kallisto objects given a column name
#'
#' @param obj a sleuth object
#' @param col_name a column to extract
#' @return a data.frame with columns \code{sample} and column \code{col_name}
get_col <- function(obj, ...) {
  n_samples <- nrow(obj$design_matrix)
  n_trans <- nrow(obj$kal[[1]]$abundance)

  #which_cols <- as.character(...)
  lapply(seq_along(obj$kal),
    function(i) {
      which_sample <- obj$sample_to_condition$sample[i]
      dplyr::select_(obj$kal[[i]]$abundance, "target_id",
        .dots = lazyeval::lazy_dots(...)) %>%
          mutate(sample = which_sample)
    }) %>%
      data.table::rbindlist() %>%
      as.data.frame()
}

#' @export
summary.sleuth <- function(obj, covariates = TRUE) {
  mapped_reads <- sapply(obj$kal, function(k) attr(k, 'num_mapped'))
  n_bs <- sapply(obj$kal, function(k) length(k$bootstrap))
  n_processed <- sapply(obj$kal, function(k) attr(k, 'num_processed'))

  res <- adf(sample = obj$sample_to_covariates[['sample']],
    reads_mapped = mapped_reads,
    reads_proc = n_processed,
    frac_mapped = round(mapped_reads / n_processed, 4),
    bootstraps = n_bs
    )
  if (covariates) {
    res <- dplyr::left_join(res, obj$sample_to_covariates, by = 'sample')
  }

  res
}

#' Create a gene table from a sleuth object
#'
#' Take in a sleuth object with added genes, and return a data table in which genes
#' list the most significant transcript mapping to themselves
#'
#' @param obj a \code{sleuth} object
#' @param test a character string denoting which beta to use
#' @param test_type either 'wt' for Wald test or 'lrt' for likelihood ratio test
#' @param which_model a character string denoting which model to use
#' @param which_group a character string denoting which gene group to use
#' @return a \code{data.frame} with the following columns
#'
#' @return gene name; if ext_gene name specified, it will be legible gene name.
#'         If ens_gene name, it will be an Ensemble gene (assuming user followed the vignette)
#' @return most_sig_trancript: Most significant transcript for the given gene
#' @return pval: p-value for the test chosen
#' @return qval: False discovery rate normalized p-value (Benjamini-Hochberg, see: \code{\link{p.adjust}})
#' @return num_transcripts: Total number of transcripts for the given gene
#' @return list_of_transcripts: All transcripts associated with this gene
#' @examples sleuth_genes <- sleuth_gene_table(so, 'conditionIP', test_type ='wt',
#'                                   which_group = 'ext_gene')
#' head(sleuth_genes) # show info for first 5 genes
#' sleuth_genes[1:5, 6] # show transcripts for first 5 genes
#' @export
sleuth_gene_table <- function(obj, test, test_type = 'lrt', which_model = 'full', which_group = 'ens_gene') {

  if (is.null(obj$target_mapping)) {
    stop("This sleuth object doesn't have added gene names.")
  }
  popped_gene_table <- sleuth_results(obj, test, test_type, which_model)

  popped_gene_table <- dplyr::arrange_(popped_gene_table, which_group, ~qval)
  popped_gene_table <- dplyr::group_by_(popped_gene_table, which_group)

  popped_gene_table <- dplyr::summarise_(popped_gene_table,
    target_id = ~target_id[1],
    pval = ~min(pval, na.rm  = TRUE),
    qval = ~min(qval, na.rm = TRUE),
    num_transcripts = ~n(),
    # all_target_ids = ~toString(target_id[1:length(target_id)])
    all_target_ids = ~paste0(target_id[1:length(target_id)], collapse = ',')
    )

  filter_empty <- nchar(popped_gene_table[[which_group]]) == 0 | # empty transcript name
    is.na(popped_gene_table[[which_group]]) | # empty gene name
    is.na(popped_gene_table$qval) # missing q-value
  filter_empty <- !filter_empty
  popped_gene_table <- popped_gene_table[filter_empty, ]

  popped_gene_table
}


#' Get the names of the transcripts associated to a gene
#'
#' Get the names of the transcripts associated to a gene, assuming genes are added to the input \code{sleuth} object.
#'
#' @param obj a \code{sleuth} object
#' @param test a character string denoting which beta to use
#' @param test_type either 'wt' for wald test or 'lrt' for likelihood ratio test
#' @param which_model a character string denoting which model to use
#' @param gene_colname the name of the column in which the desired gene apperas gene appears. Once genes have been added to a sleuth
#' object, you can inspect the genes names present in your sleuth object via \code{obj$target_mapping}, assuming 'obj' is the name of your sleuth object.
#' This parameter refers to the name of the column that the gene you are searching for appears in. Checkout the column names using \code{names(obj$target_mapping)}
#' @param gene_name a string containing the name of the gene you are interested in
#' @return a vector of strings containing the names of the transcripts that map to a gene
#' @export
transcripts_from_gene <- function(obj, test, test_type,
  which_model, gene_colname, gene_name) {

  # FIXME: this is a work around
  obj$gene_mode <- FALSE

  table <- sleuth_results(obj, test, test_type, which_model)
  table <- dplyr::select_(table, ~target_id, gene_colname, ~qval)
  table <- dplyr::arrange_(table, gene_colname, ~qval)
  if (!(gene_name %in% table[, 2])) {
      stop("Couldn't find gene ", gene_name)
  }
  table$target_id[table[, 2] == gene_name]
}

#' Change sleuth transform function
#'
#' Replace the transformation function of a sleuth object
#'
#' NOTE: if you change the transformation function after having done a fit,
#' the fit(s) will have to be redone using the new transformation.
#' @examples transform_fun(x) <- function(x) log2(x+0.5)
#' @export
`transform_fun<-` <- function(obj, fxn) {
  stopifnot(is.function(fxn))
  obj$transform_fun <- fxn
  if(!is.null(obj$fits)) {
    warning(paste("Your sleuth object has fits based on the old transform function.",
                  "Please rerun sleuth_prep and sleuth_fit."))
    obj$fits <- lapply(obj$fits, function(x) {
                         x$transform_synced <- FALSE
                         x
                       })
  }
  obj
}


#' Extend internal '$<-' for sleuth object
#'
#' This extension is mainly to address case where
#' transform_fun is changed by user.
#' This function informs user that the fits need to be redone
#' and updates those fits.
#' Otherwise it acts normally.
#' @examples obj$transform_fun <- function(x) log2(x+0.5)
#' @export
`$<-.sleuth` <- function(obj, name, value) {
  obj[[name]] <- value
  if(name=="transform_fun") {
    if(!is.null(obj$fits)) {
      warning(paste("Your sleuth object has fits based on the old transform function.",
                    "Please rerun sleuth_prep and sleuth_fit."))
      obj$fits <- lapply(obj$fits, function(x) {
                           x$transform_synced <- FALSE
                           x
                         })
    }
  }
  obj
}
