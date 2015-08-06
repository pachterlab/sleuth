#' Measurement error model with equal variances
#'
#' Fit the measurement error model with equal variances across samples and
#' noise
#'
#' @param obj a \code{sleuth} object
#' @param pass_filt a data.frame with columns \code{target_id} and
#' \code{count_filt}. \code{count_filt} is a logical where a target labeled as
#' FALSE does not pass the filter and thus will not be used in testing
#' @param xform a function to transform a matrix
#' @return a \code{data.frame}
#' @export
me_equal_var <- function(obj, pass_filt, xform = function(x) log(x + 0.5)) {
  cat("Summarizing bootstraps\n")
  bs_summary <- bs_sigma_summary(obj, xform)

  # filter things first
  filt_names <- dplyr::filter(pass_filt, count_filt)[['target_id']]

  bs_summary$obs_counts <- bs_summary$obs_counts[filt_names,]
  bs_summary$sigma_q_sq <- bs_summary$sigma_q_sq[filt_names]

  cat("Fitting ME models\n")
  mes <- me_model_by_row(obj, obj$design_matrix, bs_summary)
  tid <- names(mes)

  mes_df <- dplyr::bind_rows(lapply(mes,
    function(x) {
      data.frame(rss = x$rss, sigma_sq = x$sigma_sq, sigma_q_sq = x$sigma_q_sq,
        mean_obs = x$mean_obs, var_obs = x$var_obs)
    }))

  mes_df$target_id <- tid
  rm(tid)

  mes_df <- dplyr::mutate(mes_df, sigma_sq_pmax = pmax(sigma_sq, 0))

  # mes <- semi_join(mes, dplyr::filter(pass_filt, count_filt),
  #   by = "target_id")
  # mes <- dplyr::mutate(mes, sigma_sq_pmax = pmax(sigma_sq, 0))

  cat("Grouping by sliding window\n")
  swg <- sliding_window_grouping(mes_df, "mean_obs", "sigma_sq_pmax",
    ignore_zeroes = TRUE)

  cat("Shrinkage estimation\n")
  l_smooth <- shrink_df(swg, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, "iqr")
  l_smooth <- select(mutate(l_smooth, smooth_sigma_sq = shrink ^ 4), -shrink)

  l_smooth <- mutate(l_smooth,
    smooth_sigma_sq_pmax = pmax(smooth_sigma_sq, sigma_sq))

  X <- obj$design_matrix
  A <- solve( t(X) %*% X )

  cat("Computing var(beta)\n")
  beta_covars <- lapply(1:nrow(l_smooth),
    function(i) {
      row <- l_smooth[i,]
      with(row,
        covar_beta(smooth_sigma_sq_pmax + sigma_q_sq, X, A)
        )
    })
  names(beta_covars) <- l_smooth$target_id

  list(mes = mes, summary = l_smooth, beta_covars = beta_covars)
}


#' Compute the covariance on beta under OLS
#'
#' Compute the covariance on beta under OLS
#' @param sigma a numeric of either length 1 or nrow(X) defining the variance
#' on D_i
#' @param X the design matrix
#' @param A inv(t(X) X) (for speedup)
#' @return a covariance matrix on beta
covar_beta <- function(sigma, X, A) {
  if (length(sigma) == 1) {
    return( sigma * A )
  }

  # sammich!
  A %*% (t(X) %*% diag(sigma) %*% X) %*% A
}

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

#' non-equal var
#'
#' word
#'
#' @param obj a sleuth object
#' @param design a design matrix
#' @param samp_bs_summary the sample boostrap summary computed by sleuth_summarize_bootstrap_col
#' @return a list with a bunch of objects used for shrinkage :)
#' @export
me_heteroscedastic_by_row <- function(obj, design, samp_bs_summary, obs_counts) {
  stopifnot( is(obj, "sleuth") )

  cat("dcasting...\n")
  sigma_q_sq <- dcast(
    select(samp_bs_summary, target_id, bs_var_est_counts, sample),
    target_id ~ sample,
    value.var  = "bs_var_est_counts")
  sigma_q_sq <- as.data.frame(sigma_q_sq)
  rownames(sigma_q_sq) <- sigma_q_sq$target_id
  sigma_q_sq$target_id <- NULL
  sigma_q_sq <- as.matrix(sigma_q_sq)

  stopifnot( all.equal(rownames(sigma_q_sq), rownames(obs_counts)) )
  stopifnot( dim(sigma_q_sq) == dim(obs_counts))

  X <- design
  A <- solve(t(X) %*% X) %*% t(X)

  models <- lapply(1:nrow(bs_summary$obs_counts),
    function(i) {
      res <- me_white_model(design, obs_counts[i,], sigma_q_sq[i,], A)
      res$df$target_id = rownames(obs_counts)[i]
      res
    })
  names(models) <- rownames(obs_counts)

  models
}


me_white_model <- function(X, y, bs_sigma_sq, A) {
  n <- nrow(X)
  degrees_free <- n - ncol(X)

  ols_fit <- lm.fit(X, y)

  # estimate of sigma_i^2 + sigma_{qi}^2
  r_sq <- ols_fit$residuals ^ 2
  sigma_sq <- r_sq - bs_sigma_sq

  mean_obs <- mean(y)
  var_obs <- var(y)

  df <- data.frame(mean_obs = mean_obs, var_obs = var_obs,
    sigma_q_sq = bs_sigma_sq, sigma_sq = sigma_sq, r_sq = r_sq,
    sample = names(bs_sigma_sq))

  list(
    ols = ols_fit,
    r_sq = r_sq,
    sigma_sq = sigma_sq,
    bs_sigma_sq = bs_sigma_sq,
    mean_obs = mean_obs,
    var_obs = var_obs,
    df = df
    )
}

me_white_var <- function(df, sigma_col, sigma_q_col, X, tXX_inv) {
  # TODO: ensure X is in the same order as df
  sigma <- df[[sigma_col]] + df[[sigma_q_col]]
  df <- mutate(df, sigma = sigma)
  beta_var <- tXX_inv %*% (t(X) %*% diag(df$sigma) %*% X) %*% tXX_inv

  res <- as.data.frame(t(diag(beta_var)))
  res$target_id <- df$target_id[1]

  res
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

  list(
    ols_fit = ols_fit,
    b1 = ols_fit$coefficients[2],
    rss = rss,
    sigma_sq = sigma_sq,
    sigma_q_sq = sigma_q_sq,
    mean_obs = mean_obs,
    var_obs = var_obs
    )
}

#' @export
compute_t_me <- function(data, which_sigma, Sxx, n_data, adjust_se = NULL) {
  #var_b <- (data[,which_sigma] + data[,"sigma_q_sq"]) / (n_data * Sxx)
  var_b <- (data[,which_sigma] + data[,"sigma_q_sq"]) / (Sxx)
  se_b <- sqrt( var_b )

  s0 <- 0
  if (!is.null(adjust_se)) {
    s0 <- quantile(se_b, adjust_se, na.rm = TRUE)
    #s0 <- s0 * 10
  }
  cat(s0, "\n")
  cat(quantile(se_b, probs = seq(0, 1, length.out = 10), na.rm = TRUE), "\n")

  data$t_value <- data$b1 / (se_b + s0)
  data$se_b <- se_b
  data$pval <- 2 * pt(abs(data$t_value), data$degrees_free[1], lower.tail = FALSE)

  data
}
