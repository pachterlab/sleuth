# Constants
.MIN_VARIANCE <- 0.001

#' Fit locfit and predict
#'
#' Wrapper to fit a locfit model and predict in one call
#'
#' @param x predictor variable
#' @param y response variable
#' @param maxk See \link{code{locfit::locfit.raw}}
#' @param ... Addition arguments sent to \link{code{locfit::locfit.raw}}
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
#' @export
fit_bootstrap_bio <- function(sleu) {
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
      var_tpm_smooth = lf_predict(mean_tpm, var_tpm_raw))

  mv
}

