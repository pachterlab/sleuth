
#' Constructor for a 'sleuth' object
#'
#' Conceptually, a sleuth is a pack of kallistos. A 'sleuth' object stores
#' a pack of kallisto results, and can then intelligently operate them
#' depending on conditions and sequencing depth.
#' @param kal_list a list of \code{kallisto} objects
#' @param sample_names a character vector of \code{length(kal_list)} with
#' identifiers for the sample.
#' @param condition_names a character vector of \code{length(kal_list)} with
#' identifiers for the condition. These typically refer to the cell type (e.g.
#' tumor vs not) and should be the same for samples that are part of the same
#' experimental condition
#' @param norm_boostrap if \code{TRUE} use the size factors calculated from the
#' raw observations to normalize bootstrap values
#' @return a \code{sleuth} object
#' @export
new_sleuth <- function(kal_list, sample_names, condition_names,
  norm_bootstraps = TRUE) {
  if (length(kal_list) != length(sample_names)) {
    stop(paste0("'", substitute(kal_list), "' must be the same length as '",
        substitute(sample_names), "'"))
  }

  if (length(kal_list) != length(condition_names)) {
    stop(paste0("'", substitute(kal_list), "' must be the same length as '",
        substitute(condition_names), "'"))
  }

  # append sample ane condition columns to data
  kal_list <- lapply(seq_along(kal_list), function(it)
    {
      kal_list[[it]]$abundance <- kal_list[[it]]$abundance %>%
        mutate(sample = sample_names[it], condition = condition_names[it])
      kal_list[[it]]
    })


  obs_abundance_raw <- rbind_all(lapply(kal_list, function(k) k$abundance))

  # normalize TPM using TMM
  tpm_spread <- spread_abundance_by(obs_abundance_raw, "tpm")
  tpm_sf <- edgeR::calcNormFactors(tpm_spread)
  tpm_norm <- t(t(tpm_spread) / tpm_sf) %>%
    as.data.frame(stringsAsFactors = FALSE)
  tpm_norm$target_id <- rownames(tpm_norm)
  tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)

  # normalize est_counts using TMM
  est_counts_spread <- spread_abundance_by(obs_abundance_raw, "est_counts")
  est_counts_sf <- edgeR::calcNormFactors(est_counts_spread)
  est_counts_norm <- t(t(est_counts_spread) / est_counts_sf) %>%
    as.data.frame(stringsAsFactors = FALSE)
  est_counts_norm$target_id <- rownames(est_counts_norm)
  est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

  obs_norm <- inner_join(data.table::as.data.table(est_counts_norm),
    data.table::as.data.table(tpm_norm), by = c("target_id", "sample"))

  sample_to_condition <- data.frame(sample = sample_names,
    condition = condition_names, stringsAsFactors = FALSE)

  obs_norm <- left_join(obs_norm, data.table::as.data.table(sample_to_condition),
    by = c("sample")) %>%
    as.data.frame(stringsAsFactors = FALSE)

  print(tpm_sf)
  print(est_counts_sf)

  # Normalize all the bootstrap samples
  kal_list <- lapply(seq_along(kal_list), function(i)
    {
      normalize_bootstrap(kal_list[[i]], tpm_sf[i], est_counts_sf[i])
    })

  structure(list(
      kal = kal_list,
      obs_raw = obs_abundance_raw,
      obs_norm = obs_norm,
      sample_to_condition = sample_to_condition,
      bootstrap_summary = NA,
      tpm_sf = tpm_sf,
      est_counts_sf = est_counts_sf
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

sleuth_summarize_bootstrap_col <- function(obj, col) {
  res <- lapply(seq_along(obj$kal), function(i)
    {
      cur_samp <- obj$sample_to_condition$sample[i]
      cur_cond <- obj$sample_to_condition$condition[i]

      summarize_bootstrap(obj$kal[[i]], col) %>%
        mutate(sample = cur_samp,
          condition = cur_cond)
    })

  rbind_all(res)
}

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
