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

      suppressMessages({kal <- read_kallisto_h5(fname, read_bootstrap = TRUE)})
      kal$abundance <- dplyr::mutate(kal$abundance,
        sample = sample_to_covariates$sample[i])

      kal
    })
  msg('')

  obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))

  ret <- list(
      kal = kal_list,
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
    # filter_bool <- filter_df_all_groups(est_counts_spread, filter_fun,
    #   sample_to_covariates)
    filter_true <- filter_bool[filter_bool]

    msg(paste0(sum(filter_bool), ' targets passed the filter'))
    est_counts_sf <- norm_factors(est_counts_spread[filter_bool,])

    filter_df <- adf(target_id = names(filter_true))

    est_counts_norm <- as_df(t(t(est_counts_spread) / est_counts_sf))

    est_counts_norm$target_id <- rownames(est_counts_norm)
    est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

    obs_norm <- est_counts_norm
    obs_norm$target_id <- as.character(obs_norm$target_id)
    rm(est_counts_norm)

    # deal w/ TPM
    msg("normalizing tpm")
    tpm_spread <- spread_abundance_by(obs_raw, "tpm")
    tpm_sf <- norm_factors(tpm_spread[filter_bool,])
    tpm_norm <- as_df(t(t(tpm_spread) / tpm_sf))
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)

    suppressWarnings({
      if ( !all.equal(dplyr::select(obs_norm, target_id, sample),
          dplyr::select(tpm_norm, target_id, sample)) ) {
        stop('Invalid column rows. In principle, can simply join. Please report error.')
      }

      # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
      #   by = c('target_id', 'sample'))
      obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, tpm))
      obs_norm <- dplyr::left_join(obs_norm,
        data.table::as.data.table(sample_to_covariates),
        by = c("sample"))
    })

    msg("normalizing bootstrap samples")
    kal_list <- lapply(seq_along(kal_list), function(i) {
      normalize_bootstrap(kal_list[[i]],
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
  mapped_reads <- sapply(obj$kal, function(k) sum(k$abundance$est_counts))
  n_bs <- sapply(obj$kal, function(k) length(k$bootstrap))

  res <- adf(sample = obj$sample_to_covariates[['sample']],
    mapped_reads = mapped_reads,
    n_bootstraps = n_bs)
  if (covariates) {
    res <- dplyr::left_join(res, obj$sample_to_covariates, by = 'sample')
  }

  res
}

