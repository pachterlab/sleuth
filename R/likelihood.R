# compute the likelihood of which_model, then return the original object with
# new element: obj$fits[[which_model]]$likelihood
compute_likelihood <- function(obj, which_model) {
  stopifnot(is(obj, "sleuth"))
  model_exists(obj, which_model)

  # we basically do lapply on all of the models
  #
  # the fitted values are here:
  #   obj$fits[[which_model]]$models$ols_fit$fitted.values
  #
  # the observations can be recovered by:
  #   obj$fits$full$models[[1]]$ols_fit$residuals + fitted values

  # TODO: move this elsewhere
  obj$fits[[which_model]]$summary <- obj$fits[[which_model]]$summary[
    match(names(obj$fits[[which_model]]$models),
      obj$fits[[which_model]]$summary$target_id), ]

  all_likelihood <- sapply(seq_along(obj$fits[[which_model]]$models),
    function( i ) {
      cur_model <- obj$fits[[which_model]]$models[[ i ]]
      cur_mu <- cur_model$ols_fit$fitted.values
      obs <- cur_model$ols_fit$residuals + cur_mu

      cur_summary <- obj$fits[[which_model]]$summary

      cur_var <- cur_summary[i, "smooth_sigma_sq_pmax"] +
        cur_summary[i, "sigma_q_sq"]

      sum(dnorm(obs, mean = cur_mu, sd = sqrt(cur_var), log = TRUE))
    })

  names(all_likelihood) <- names(obj$fits[[which_model]]$models)

  obj$fits[[which_model]]$likelihood <- all_likelihood

  obj
}

get_likelihood <- function(obj, which_model) {
  stopifnot( is(obj, "sleuth") )

   obj$fits[[which_model]]$likelihood
}

likelihood_exists <- function(obj, which_model) {
  !is.null( get_likelihood(obj, which_model) )
}

#'  sleuth likelihood ratio test
#'
#' compute the likelihood ratio test for 2 models. this requires that the null
#' model be nested in the alternate model
#'
#' @param obj a sleuth object
#' @param null_model the null (or "reduced") model
#' @param alt_model the alternate (or "full") model
#' @return an updated sleuth object with a likelihood ratio test computation
#' @seealso \code{\link{models}} to view which models have been fit and which
#' coefficients can be tested, \code{\link{sleuth_results}} to get back
#' a data.frame of the results
#' @export
sleuth_lrt <- function(obj, null_model, alt_model) {
  stopifnot( is(obj, "sleuth") )
  model_exists(obj, null_model)
  model_exists(obj, alt_model)

  if(!obj$fits[[alt_model]]$transform_synced) {
    stop("Model '", alt_model, "' was not computed using the sleuth object's",
         " current transform function. Please rerun sleuth_fit for this model.")
  }

  if(!obj$fits[[null_model]]$transform_synced) {
    stop("Model '", null_model, "' was not computed using the sleuth object's",
         " current transform function. Please rerun sleuth_fit for this model.")
  }

  if ( !likelihood_exists(obj, null_model) ) {
    obj <- compute_likelihood(obj, null_model)
  }
  if ( !likelihood_exists(obj, alt_model) ) {
    obj <- compute_likelihood(obj, alt_model)
  }

  n_ll <- get_likelihood(obj, null_model)
  a_ll <- get_likelihood(obj, alt_model)

  obj$fits[[null_model]]$likelihood <- n_ll
  obj$fits[[alt_model]]$likelihood <- a_ll

  test_statistic <- 2 * (a_ll - n_ll)

  degrees_free <- obj$fits[[null_model]]$models[[1]]$ols_fit$df.residual -
    obj$fits[[alt_model]]$models[[1]]$ols_fit$df.residual

  # P(chisq > test_statistic)
  p_value <- pchisq(test_statistic, degrees_free, lower.tail = FALSE)
  result <- adf(target_id = names(obj$fits[[alt_model]]$likelihood),
    test_stat = test_statistic, pval = p_value)
  result <- dplyr::mutate(result, qval = p.adjust(pval, method = "BH"))
  model_info <- data.table::data.table(obj$fits[[null_model]]$summary)
  model_info <- dplyr::select(model_info, -c(iqr))
  result <- dplyr::left_join(
    data.table::data.table(result),
    model_info,
    by = 'target_id')
  result <- dplyr::mutate(result, degrees_free = degrees_free)

  test_name <- paste0(null_model, ':', alt_model)
  obj <- add_test(obj, result, test_name, 'lrt')

  obj
}
