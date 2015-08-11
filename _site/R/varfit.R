# Constants
.MIN_VARIANCE <- 0.001

#' Fit locfit and predict
#'
#' Wrapper to fit a locfit model and predict in one call
#'
#' @param x predictor variable
#' @param y response variable
#' @param maxk See \code{\link{locfit::locfit.raw}}
#' @param ... Addition arguments sent to \code{\link{locfit::locfit.raw}}
#' @return predicted values from fitting y ~ x
lf_predict <- function(x, y, maxk = 1000, ...) {
  locfit.raw(x, y, maxk = maxk, ...) %>%
    predict(x)
}

#' Fit the variance as a function of the mean
#'
#' Fit the variance as a function of the mean, by conditioning on the boostrap
#' samples.
#'
#' @param sleu a \code{sleuth} object
#' @param pool if \code{TRUE} pool all conditions together, otherwise only pool
#' by condition. NOTE: \code{pool = FALSE} is not yet implemented
#' @return a data frame with the fitted mean and variance
#' @export
varfit_smooth_bio <- function(sleu, pool = TRUE) {
  stopifnot( is(sleu, "sleuth") )

  sleu <- sleuth_summarize_bootstrap(sleu)

  # summary of bootstrap
  data <- inner_join(
    data.table::data.table(sleu$bootstrap_summary),
    data.table::data.table(sleu$obs_norm),
    by = c("target_id", "sample", "condition")
    )

  # compute the summary statistics over samples
  mv <- data %>%
    group_by(target_id) %>%
    summarise(
      mean_tpm = mean(tpm),
      mean_tpm_bs = mean(bs_mean_tpm),
      var_tpm_raw = var(tpm),
      var_tpm_bs = mean(bs_var_tpm)
      ) %>%
    mutate(
      var_tpm_bio_raw = var_tpm_raw - var_tpm_bs,
      var_tpm_bio_raw = ifelse(var_tpm_bio_raw < .MIN_VARIANCE,
        .MIN_VARIANCE,
        var_tpm_bio_raw)
      )

  # compute the smooth estimates
  mv <- mv %>%
    mutate(
      var_tpm_bio_smooth = lf_predict(mean_tpm, var_tpm_bio_raw),
      var_tpm_bio_smooth_total = var_tpm_bio_smooth + var_tpm_bs,
      var_tpm_smooth = lf_predict(mean_tpm, var_tpm_raw),
      raw_gt_bs = var_tpm_raw > var_tpm_bs)

  mv
}

#' Fit estimated variance from a bunch of bootstraps
#'
#' Take a bunch of bootstrap samples and compute the Monte Carlo estimate of
#' the variance using bootstrap samples.
#'
#' @param obj a \code{sleuth} object
#' @return a data frame with a summary statistics
#' @export
varfit_bootstrap_est <- function(obj, pool = TRUE) {
  stopifnot(is(obj, "sleuth"))

  cat("Melting sleuth object\n")
  bs <- data.table(melt_bootstrap_sleuth(obj), key = c("target_id", "bs_sample"))

  cat("Summarizing bootstraps by sample\n")
  bs <- bs[, list(mean_tpm_bs = mean(tpm), var_tpm_bs = var(tpm)),
    by = c("target_id", "bs_sample")]

  cat("Summarizing across bootstraps\n")
  bs <- bs[, list(mean_tpm_bs_pooled = mean(mean_tpm_bs),
      var_tpm_bs_pooled = mean(var_tpm_bs)), by = c("target_id")]
  # bs <- bs %>%
  #   group_by(target_id, bs_sample) %>%
  #   summarise(
  #     mean_tpm_bs = mean(tpm),
  #     var_tpm_bs = var(tpm)
  #     )

  # bs %>%
  #   group_by(target_id) %>%
  #   summarise(mean_tpm_bs_pooled = mean(mean_tpm_bs),
  #     var_tpm_bs_pooled = mean(var_tpm_bs))
  bs
}
