#' Basic row filter
#'
#' A basic filter to be used.
#'
#' @param row this is a vector of numerics that will be passedin
#' @param mean_reads the minimum mean number of reads
#' @param min_prop the minimum proportion of reads to pass this filter
#' @return a logical of length 1
#' @export
basic_filter <- function(row, mean_reads = 5, min_prop = 0.8) {
  mean(row > mean_reads) > min_prop
}

#' Constructor for a 'sleuth' object
#'
#' Conceptually, a sleuth is a pack of kallistos. A 'sleuth' object stores
#' a pack of kallisto results, and can then intelligently operate them
#' depending on conditions and sequencing depth.
#'
#' @param kal_dirs a character vector of length greater than one where each
#' string points to a kallisto output directory
#' @param sample_to_covariates is a \code{data.frame} which contains a mapping
#' from \code{sample} (a column) to some set of experimental conditions or
#' covariates
#' @param full_model is a \code{formula} which explains the full model (design)
#' of the experiment. It must be consistent with the data.frame supplied in
#' \code{sample_to_covariates}
#' @param trans is a character string of either 'log', or 'cpm'. Currently
#' ignored.
#' @param filter_fun the function to use when filtering.
#' @param ... additional arguments passed to other functions
#' @return a \code{sleuth} object containing all kallisto samples, metadata,
#' and summary statistics
#' @seealso \code{\link{sleuth_fit}} to fit a model, \code{\link{wald_test}} to
#' test whether a coeffient is zero
#' @export
sleuth_prep <- function(
  kal_dirs,
  sample_names,
  sample_to_covariates,
  full_model,
  normalize = TRUE,
  trans = 'log',
  filter_fun = basic_filter,
  annotation = NULL,
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
        "' must contain a column names 'sample'"))
  }

  # TODO: check if paths exists and contain a HDF5 file
  # TODO: ensure all kallisto have same number of transcripts
  # TODO: ensure transcripts are in same order -- if not, sort them

  # done
  ##############################

  msg('Reading in kallisto results')

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

  if (FALSE) {
    # ignore TPM stuff for now
    msg("Normalizing 'tpm'")

    tpm_spread <- spread_abundance_by(obs_raw, "tpm")
    tpm_sf <- DESeq2::estimateSizeFactorsForMatrix(tpm_spread)
    tpm_norm <- t(t(tpm_spread) / tpm_sf) %>%
      as.data.frame(stringsAsFactors = FALSE)
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
  }

  ret <- list(
      kal = kal_list,
      obs_raw = obs_raw,
      sample_to_covariates = sample_to_covariates,
      bootstrap_summary = NA,
      full_formula = full_model,
      design_matrix = model.matrix(full_model, sample_to_covariates)
    )

  if ( normalize ) {
    msg("Normalizing 'est_counts'")
    est_counts_spread <- spread_abundance_by(obs_raw, "est_counts")
    filter_bool <- apply(est_counts_spread, 1, filter_fun)
    filter_true <- filter_bool[filter_bool]

    msg(paste0(sum(filter_bool), ' targets passed the filter.'))
    est_counts_sf <- DESeq2::estimateSizeFactorsForMatrix(est_counts_spread[filter_bool,])

    filter_df <- adf(target_id = names(filter_true))

    est_counts_norm <- as_df(t(t(est_counts_spread[filter_bool,]) / est_counts_sf))
    est_counts_norm$target_id <- rownames(est_counts_norm)
    est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

    obs_norm <- est_counts_norm

    obs_norm <- dplyr::left_join(obs_norm,
      data.table::as.data.table(sample_to_covariates),
      by = c("sample"))
    obs_norm <- as_df(obs_norm)

    msg("Normalizing bootstrap samples")
    kal_list <- lapply(seq_along(kal_list), function(i) {
        normalize_bootstrap(kal_list[[i]],
          est_counts_size_factor = est_counts_sf[i])
      })

    # add relevant objs back into sleuth obj
    ret$obs_norm <- obs_norm
    ret$est_counts_sf <- est_counts_sf
    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
  }

  class(ret) <- 'sleuth'

  ret
}

#' Summarize many bootstrap objects
#'

#' Summarize all the bootstrap samples from a kallisto run. The summarized
#' values are then all put into a data.frame and stored in the \code{sleuth}
#' object.
#'
#' @param obj a \code{sleuth} object
#' @param force if \code{FALSE}, then will only compute the summary if it has
#' not yet been set.
#' @param verbose if \code{TRUE}, print verbosely
#' @return a \code{kallisto} object with member \code{bootstrap_summary}
#' updated and containing a data frame.
#' @export
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
    by = c("target_id", "sample", "condition")
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  obj$bootstrap_summary <- s_bs

  obj
}

var_fit <- function(obj) {
  stopifnot(is(obj, "sleuth"))
  all_data <- inner_join(
    data.table::data.table(obj$bootstrap_summary),
    data.table::data.table(obj$obs_norm),
    by = c("target_id", "sample", "condition")
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  all_data
}

sleuth_summarize_bootstrap_col <- function(obj, col, transform = identity) {
  res <- lapply(seq_along(obj$kal), function(i)
    {
      cur_samp <- obj$sample_to_covariates$sample[i]
      cur_cond <- obj$sample_to_covariates$condition[i]

      dplyr::mutate(summarize_bootstrap(obj$kal[[i]], col, transform),
        sample = cur_samp, condition = cur_cond)
    })

  dplyr::bind_rows(res)
}

#' Spread abundance by a column
#'
#' Take a data.frame from a sleuth object (e.g. \code{obs_raw}) and cast it
#' into a matrix where the rows are the target_ids and the columns are the
#' sample ids. The values are the variable you are "spreading" on.
#' @param abund the abundance \code{data.frame} from a \code{sleuth} object
#' @param var a character array of length one. The variable for which to get "spread" on (e.g. "est_counts").
#' @export
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
obs_transcript_summary <- function(data, pool = TRUE)
{
  stopifnot(is(data, "sleuth"))

  obs_abundance <- data$obs_abundance

  if (pool) {
    obs_abundance %>%
      group_by(target_id) %>%
      summarize(
        mean_tpm = mean(tpm),
        var_tpm = var(tpm),
        mean_est_counts = mean(est_counts),
        var_est_counts = var(est_counts)
        )
  }

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

#' @export
null_mean_var <- function(obj, transform = identity, min_reads = 1) {
  stopifnot( is(obj, "sleuth") )

  which_pass <- obj$obs_norm %>%
    group_by(target_id) %>%
    summarise(pass = all(est_counts > min_reads)) %>%
    filter(pass)

  obj$obs_norm %>%
    data.table::data.table() %>%
    inner_join(data.table::data.table(which_pass), by = c("target_id")) %>%
    mutate(trans_counts = transform(est_counts)) %>%
    group_by(target_id) %>%
    summarise(
      counts_mean = mean(trans_counts),
      counts_var = var(trans_counts)
      )
}

#' observations to a matrix
#'
#' observations to a matrix
#'
#' @param obj is a sleuth object
#' @param value_name either "est_counts" or "tpm"
#' @return a matrix with the appropriate names
obs_to_matrix <- function(obj, value_name) {

  obs_counts <- reshape2::dcast(obj$obs_norm, target_id ~ sample,
    value.var = value_name)

  obs_counts <- as.data.frame(obs_counts)
  rownames(obs_counts) <- obs_counts$target_id
  obs_counts$target_id <- NULL
  obs_counts <- as.matrix(obs_counts)

  obs_counts
}


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
