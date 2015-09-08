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

# currently defunct
filter_df_all_groups <- function(df, fun, group_df, ...) {
  grps <- setdiff(colnames(group_df), 'sample')
  res <- sapply(unique(grps),
    function(g) {
      apply(filter_df_by_groups(df, fun, group_df[,c('sample', g)], ...), 1, any)
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
      apply(df[,valid_samps], 1, fun, ...)
    })
}

#' Constructor for a 'sleuth' object
#'
#' A sleuth is a group of kallistos. Borrowing this terminology, a 'sleuth' object stores
#' a group of kallisto results, and can then operate on them while
#' accounting for covariates, sequencing depth, technical and biological
#' variance.
#'
#' @param kal_dirs a character vector of length greater than one where each
#' string points to a kallisto output directory
#' @param sample_to_covariates is a \code{data.frame} which contains a mapping
#' from \code{sample} (a column) to some set of experimental conditions or
#' covariates. The column \code{sample} should be in the same order as the
#' corresponding entry in \code{kal_dirs}
#' @param full_model is a \code{formula} which explains the full model (design)
#' of the experiment. It must be consistent with the data.frame supplied in
#' \code{sample_to_covariates}
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
#' @return a \code{sleuth} object containing all kallisto samples, metadata,
#' and summary statistics
#' @seealso \code{\link{sleuth_fit}} to fit a model, \code{\link{sleuth_test}} to
#' test whether a coeffient is zero
#' @export
sleuth_prep <- function(
  kal_dirs,
  sample_to_covariates,
  full_model,
  filter_fun = basic_filter,
  target_mapping = NULL,
  max_bootstrap = NULL,
  ...) {

  ##############################
  # check inputs

  # data types
  if( !is(kal_dirs, 'character') ) {
    stop(paste0('"', substitute(kal_dirs),
        '" (kal_dirs) is must be a character vector.'))
  }

  if ( !is(sample_to_covariates, "data.frame") ) {
    stop(paste0("'", substitute(sample_to_covariates), "' must be a data.frame"))
  }

  if (!is(full_model, "formula")) {
    stop(paste0("'",substitute(full_model), "' (full_model) must be a formula"))
  }

  if (length(kal_dirs) != nrow(sample_to_covariates)) {
    stop(paste0("'", substitute(kal_dirs),
        "' must have the same length as the number of rows in '",
        substitute(sample_to_covariates), "'"))
  }

  if (!("sample" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates),
        "' (sample_to_covariates) must contain a column named 'sample'"))
  }

  if ( !is.null(target_mapping) && !is(target_mapping, 'data.frame')) {
    stop(paste0("'", substitute(target_mapping),
        "' (target_mapping) must be a data.frame or NULL"))
  } else if (is(target_mapping, 'data.frame')){
    if (!("target_id" %in% colnames(target_mapping))) {
      stop(paste0("'", substitute(target_mapping),
          "' (target_mapping) must contain a column named 'target_id'"))
    }
  }

  if ( !is.null(max_bootstrap) && max_bootstrap <= 0 ) {
    stop("max_bootstrap must be > 0")
  }

  # TODO: ensure all kallisto have same number of transcripts
  # TODO: ensure transcripts are in same order -- if not, report warning that
  # kallisto index might be correct

  # done
  ##############################

  msg('reading in kallisto results')
  sample_to_covariates$sample <- as.character(sample_to_covariates$sample)

  nsamp <- 0
  # append sample column to data
  kal_list <- lapply(seq_along(kal_dirs),
    function(i) {
      fname <- file.path(kal_dirs[i], 'abundance.h5')

      nsamp <- dot(nsamp)

      if ( !file.exists(fname) ) {
        stop(paste0('Could not find HDF5 file: ', fname))
      }

      suppressMessages({kal <- read_kallisto_h5(fname, read_bootstrap = TRUE, max_bootstrap = max_bootstrap)})
      kal$abundance <- dplyr::mutate(kal$abundance,
        sample = sample_to_covariates$sample[i])

      kal
    })
  msg('')

  kal_versions <- check_kal_pack(kal_list)

  obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))

  obs_raw <- dplyr::arrange(obs_raw, target_id, sample)
  ret <- list(
      kal = kal_list,
      kal_versions = kal_versions,
      obs_raw = obs_raw,
      sample_to_covariates = sample_to_covariates,
      bootstrap_summary = NA,
      full_formula = full_model,
      design_matrix = model.matrix(full_model, sample_to_covariates),
      target_mapping = target_mapping
    )

  # TODO: eventually factor this out
  normalize <- TRUE
  if ( normalize ) {
    msg("normalizing est_counts")
    est_counts_spread <- spread_abundance_by(obs_raw, "est_counts")
    filter_bool <- apply(est_counts_spread, 1, filter_fun, ...)
    filter_true <- filter_bool[filter_bool]

    msg(paste0(sum(filter_bool), ' targets passed the filter'))
    est_counts_sf <- norm_factors(est_counts_spread[filter_bool,])

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
    tpm_spread <- spread_abundance_by(obs_raw, "tpm")
    tpm_sf <- norm_factors(tpm_spread[filter_bool,])
    tpm_norm <- as_df(t(t(tpm_spread) / tpm_sf))
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
    tpm_norm$sample <- as.character(tpm_norm$sample)

    msg('merging in metadata')
    # put everyone in the same order to avoid a slow join
    obs_norm <- dplyr::arrange(obs_norm, target_id, sample)
    tpm_norm <- dplyr::arrange(tpm_norm, target_id, sample)

    stopifnot( all.equal(obs_raw$target_id, obs_norm$target_id) &&
      all.equal(obs_raw$sample, obs_norm$sample) )

    suppressWarnings({
      if ( !all.equal(dplyr::select(obs_norm, target_id, sample),
          dplyr::select(tpm_norm, target_id, sample)) ) {
        stop('Invalid column rows. In principle, can simply join. Please report error.')
      }

      # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
      #   by = c('target_id', 'sample'))
      obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, tpm))
    })



    # add in eff_len and len
    obs_norm <- dplyr::mutate(obs_norm,
      eff_len = obs_raw$eff_len,
      len = obs_raw$len)

    # suppressWarnings({
    #   obs_norm <- dplyr::inner_join(
    #     data.table::as.data.table(obs_norm),
    #     data.table::as.data.table(
    #       dplyr::select(obs_raw, target_id, sample, len, eff_len)
    #       ),
    #     by = c('target_id', 'sample')
    #     )
    # })

    msg("normalizing bootstrap samples")
    ret$kal <- lapply(seq_along(ret$kal), function(i) {
      normalize_bootstrap(ret$kal[[i]],
        tpm_size_factor = tpm_sf[i],
        est_counts_size_factor = est_counts_sf[i])
      })

    obs_norm <- as_df(obs_norm)
    ret$obs_norm <- obs_norm
    ret$est_counts_sf <- est_counts_sf
    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
    ret$obs_norm_filt <- dplyr::semi_join(obs_norm, filter_df, by = 'target_id')
    ret$tpm_sf <- tpm_sf
  }

  class(ret) <- 'sleuth'

  ret
}


# check versions of kallistos and num bootstraps, etc
check_kal_pack <- function(kal_list) {

  versions <- sapply(kal_list, function(x) attr(x, 'kallisto_version'))
  u_versions <- unique(versions)
  if ( length(u_versions) > 1) {
    warning('More than one version of kallisto was used: ', u_versions)
  }

  ntargs <- sapply(kal_list, function(x) attr(x, 'num_targets'))
  u_ntargs <- unique(ntargs)
  if ( length(u_ntargs) > 1 ) {
    stop('Inconsistent number of transcripts. Please make sure you used the same index everywhere.')
  }

  u_versions
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
  mat_nz <- mat[nz,]
  p <- ncol(mat)
  geo_means <- exp(apply(mat_nz, 1, function(row) (1/p) * sum(log(row)) ))
  s <- sweep(mat_nz, 1, geo_means, `/`)

  apply(s, 2, median)
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
  res <- lapply(seq_along(obj$kal), function(i)
    {
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
spread_abundance_by <- function(abund, var) {
  # var <- lazyeval::lazy(var)
  var_spread <- abund %>%
    select_("target_id", "sample", var) %>%
    tidyr::spread_("sample", var) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL

  as.matrix(var_spread)
}

#' @export
melt_bootstrap_sleuth <- function(obj) {
  # TODO: make this into a S3 function
  lapply(seq_along(obj$kal), function(i)
    {
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
    function(i)
    {
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
  n_processed = sapply(obj$kal, function(k) attr(k, 'num_processed'))

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
#' @param which_beta a character string denoting which beta to use
#' @param which_model a character string denoting which model to use
#' @param which_group a character string denoting which gene group to use
#' @return a \code{data.frame} containing gene names, transcript names, and significance
#' @export


sleuth_gene_table <- function(obj, which_beta, which_model = 'full', which_group = 'ens_gene') {
  
    if(is.null(obj$target_mapping))
    {
        stop("This sleuth object doesn't have added gene names.")
    }
    popped_gene_table = sleuth_results(obj, which_beta, which_model)

    
    popped_gene_table = dplyr::arrange_(popped_gene_table, which_group, ~qval)
    popped_gene_table = dplyr::group_by_(popped_gene_table, which_group)

    popped_gene_table = dplyr::summarise_(popped_gene_table, most_sig_transcript = ~target_id[1], pval = ~min(pval, na.rm  = TRUE), qval = ~min(qval, na.rm = TRUE), num_transcripts = ~n(), list_of_transcripts = ~toString(target_id[1:length(target_id)]))
    
    popped_gene_table = popped_gene_table[!popped_gene_table[,1] == "",]
    popped_gene_table = popped_gene_table[!is.na(popped_gene_table[,1]),] #gene_id
    popped_gene_table = popped_gene_table[!is.na(popped_gene_table$qval),]
    popped_gene_table
}


#' Get the names of the transcripts associated to a gene
#'
#' Get the names of the transcripts associated to a gene, assuming genes are added to the input \code{sleuth} object.
#'
#' @param obj a \code{sleuth} object
#' @param which_beta a character string denoting which beta to use
#' @param which_model a character string denoting which model to use
#' @param gene_colname the name of the column in which the desired gene apperas gene appears. Once genes have been added to a sleuth
#' object, you can inspect the genes names present in your sleuth object via \code{obj$target_mapping}, assuming 'obj' is the name of your sleuth object.
#' This parameter refers to the name of the column that the gene you are searching for appears in. Checkout the column names using \code{names(obj$target_mapping)}
#' @param gene_name a string containing the name of the gene you are interested in
#' @return a vector of strings containing the names of the transcripts that map to a gene
#' @export


transcripts_from_gene <- function(obj, which_beta, which_model, gene_colname, gene_name)
{
    table = sleuth_results(obj, which_beta, which_model)
    table = dplyr::select_(table, ~target_id, gene_colname, ~qval)
    table = dplyr::arrange_(table, gene_colname, ~qval)
    if(!(gene_name %in% table[,2]))
    {
        stop("Couldn't find gene ", gene_name)
    }
    table$target_id[table[,2] == gene_name]
}