# this is an internal function, you know it carries the prefix "sleuth"
sleuth_likelihood <- function(obj, which_model) {
  stopifnot(is(obj, "sleuth"))
  model_exists(obj, which_model)

  # we basically do lapply on all of the models
  #
  # the fitted values are here:
  #   obj$fits[[which_model]]$models$ols_fit$fitted.values
  #
  # TODO: after computing bootstrap summary once, get the observations from
  # there
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

  all_likelihood
}

sleuth_lrt <- function() {

}