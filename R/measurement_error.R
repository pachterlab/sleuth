#' Measurement error model
#'
#' Fit the measurement error model across all samples
#'
#' @param obj a \code{sleuth} object
#' @param design a design matrix
#' @param bs_summary a list from \code{bs_sigma_summary}
#' @return a list with a bunch of objects that are useful for shrinking
#' @export
me_model_by_row <- function(obj, design, bs_summary) {
  stopifnot( is(obj, "sleuth") )

  stopifnot( all.equal(names(bs_summary$sigma_q_sq), rownames(bs_summary$obs_counts)) )
  stopifnot( length(bs_summary$sigma_q_sq) == nrow(bs_summary$obs_counts))

  models <- lapply(1:nrow(bs_summary$obs_counts),
    function(i)
    {
      me_model(design, bs_summary$obs_counts[i,], bs_summary$sigma_q_sq[i])
    })
  names(models) <- rownames(bs_summary$obs_counts)

  models
}

#' @export
bs_sigma_summary <- function(obj, transform = identity) {
  obs_counts <- obs_to_matrix(obj, "est_counts")
  obs_counts <- transform( obs_counts )

  bs_summary <- sleuth_summarize_bootstrap_col(obj, "est_counts", transform)
  bs_summary <- bs_summary %>%
    group_by(target_id) %>%
    summarise(sigma_q_sq = mean(bs_var_est_counts))

  bs_summary <- as.data.frame(bs_summary)

  bs_sigma <- bs_summary$sigma_q_sq
  names(bs_sigma) <- bs_summary$target_id
  bs_sigma <- bs_sigma[rownames(obs_counts)]

  list(obs_counts = obs_counts, sigma_q_sq = bs_sigma)
}

me_model <- function(X, y, sigma_q_sq)
{
  n <- nrow(X)
  degrees_free <- n - ncol(X)

  ols_fit <- lm.fit(X, y)
  rss <- sum(ols_fit$residuals ^ 2)
  sigma_sq <- rss / (degrees_free) - sigma_q_sq

  mean_obs <- mean(y)
  var_obs <- var(y)

  data.frame(
    b1 = ols_fit$coefficients[2],
    rss = rss,
    sigma_sq = sigma_sq,
    sigma_q_sq = sigma_q_sq,
    mean_obs = mean_obs,
    var_obs = var_obs,
    degrees_free = degrees_free
    )
}

#' @export
compute_t_me <- function(data, which_sigma, Sxx, n_data) {

  var_b <- (data[,which_sigma] + data[,"sigma_q_sq"]) / (n_data * Sxx)
  se_b <- sqrt( var_b )

  data$t_value <- data$b1 / se_b
  data$se_b <- se_b
  data$pval <- 2 * pt(abs(data$t_value), data$degrees_free[1], lower.tail = FALSE)

  data
}
