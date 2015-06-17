
#' Constructor for a 'sleuth' object
#'
#' Conceptually, a sleuth is a pack of kallistos. A 'sleuth' object stores
#' a pack of kallisto results, and can then intelligently operate them
#' depending on conditions and sequencing depth.
#'
#' @param kal_list a list of \code{kallisto} objects
#' @param sample_to_condition is a \code{data.frame} which contains a mapping from \code{sample} (a column) to some set of experimental conditions
#' @param design is a \code{formula} which explains the full model (design) of the experiment
#' @param norm_boostrap if \code{TRUE} use the size factors calculated from the
#' raw observations to normalize bootstrap values
#' @return a \code{sleuth} object
#' @export
new_sleuth <- function(
  kal_list, sample_to_condition,
  design,
  norm_bootstraps = TRUE, verbose = TRUE) {

  ##############################
  # check inputs

  # data types
  if (!all(unlist(lapply(kal_list, is, "kallisto")))) {
    stop(paste0("One or more objects in '", substitute(kal_list),
        "' (kal_list) is NOT a 'kallisto' object"))
  }

  if (!is(sample_to_condition, "data.frame")) {
    stop(paste0("'", substitute(sample_to_condition), "' must be a data.frame"))
  }

  if (!is(design, "formula")) {
    stop(paste0("'",substitute(design), "' (design) must be a formula"))
  }

  if (length(kal_list) != nrow(sample_to_condition)) {
    stop(paste0("'", substitute(kal_list), "' must have the same length as the number of rows in '",
        substitute(sample_to_condition), "'"))
  }


  if (!("sample" %in% colnames(sample_to_condition))) {
    stop(paste0("'", substitute(sample_to_condition),
        "' must contain a column names 'sample'"))
  }

  # TODO: ensure all kallisto have same number of transcripts
  # TODO: ensure transcripts are in same order -- if not, sort them

  # done
  ##############################

  if (verbose) cat("Appending sample names\n")

  # append sample ane condition columns to data
  kal_list <- lapply(seq_along(kal_list), function(it)
    {
      kal_list[[it]]$abundance <- kal_list[[it]]$abundance %>%
        mutate(sample = sample_to_condition$sample[it])
      kal_list[[it]]
    })


  obs_abundance_raw <- rbind_all(lapply(kal_list, function(k) k$abundance))

  if (verbose) cat("Normalizing 'tpm'\n")
  # normalize TPM using TMM
  tpm_spread <- spread_abundance_by(obs_abundance_raw, "tpm")
  # tpm_sf <- edgeR::calcNormFactors(tpm_spread)
  tpm_sf <- DESeq2::estimateSizeFactorsForMatrix(tpm_spread)
  tpm_norm <- t(t(tpm_spread) / tpm_sf) %>%
    as.data.frame(stringsAsFactors = FALSE)
  tpm_norm$target_id <- rownames(tpm_norm)
  tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)

  if (verbose) cat("Normalizing 'est_counts'\n")
  # normalize est_counts using TMM
  est_counts_spread <- spread_abundance_by(obs_abundance_raw, "est_counts")
  # est_counts_sf <- edgeR::calcNormFactors(est_counts_spread)
  est_counts_sf <- DESeq2::estimateSizeFactorsForMatrix(est_counts_spread)
  est_counts_norm <- t(t(est_counts_spread) / est_counts_sf) %>%
    as.data.frame(stringsAsFactors = FALSE)
  est_counts_norm$target_id <- rownames(est_counts_norm)
  est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

  if (verbose) cat("Joining tpm and est_counts tables\n")
  obs_norm <- inner_join(data.table::as.data.table(est_counts_norm),
    data.table::as.data.table(tpm_norm), by = c("target_id", "sample"))

  # sample_to_condition <- data.frame(sample = sample_names,
  #   condition = condition_names, stringsAsFactors = FALSE)

  if (verbose) cat("Joining tpm and est_counts tables\n")
  obs_norm <- left_join(obs_norm,
    data.table::as.data.table(sample_to_condition),
    by = c("sample")) %>%
    as.data.frame(stringsAsFactors = FALSE)

  # Normalize all the bootstrap samples
  if (norm_bootstraps) {
    if (verbose) cat("Normalizing bootstrap samples\n")
    kal_list <- lapply(seq_along(kal_list), function(i)
      {
        normalize_bootstrap(kal_list[[i]], tpm_sf[i], est_counts_sf[i])
      })
  }

  structure(list(
      kal = kal_list,
      obs_raw = obs_abundance_raw,
      obs_norm = obs_norm,
      sample_to_condition = sample_to_condition,
      bootstrap_summary = NA,
      tpm_sf = tpm_sf,
      est_counts_sf = est_counts_sf,
      design = design,
      design_matrix = model.matrix(design, sample_to_condition)
      ),
    class = "sleuth")
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
      cur_samp <- obj$sample_to_condition$sample[i]
      cur_cond <- obj$sample_to_condition$condition[i]

      summarize_bootstrap(obj$kal[[i]], col, transform) %>%
        mutate(sample = cur_samp,
          condition = cur_cond)
    })

  rbind_all(res)
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
