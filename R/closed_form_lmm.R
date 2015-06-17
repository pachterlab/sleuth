#' Simple special case LMM
#'
#' This is an implementation of the special linear mixed effects model. The
#' contraints are the design is balanced (within group sizes are equal) and
#' there is one common random intercept.
#'
#' Notation is borrowed from Ch2: "Mixed Models: Theory and Applications with R,
#' Second Edition" by Eugene Demidenko
#'
#' @param X the design matrix for the fixed effect
#' @param y the response
#' @param nsamples the number of biological samples (in LMM terms, the "number
#' of clusters")
#' @param adjust_sigma if TRUE, adjust sigma when computing test
#' @export
lmm <- function(X, y, nsamples, adjust_sigma = TRUE) {
  # let's assume that the groups all have equal size

  # solve betas by QR decomposition (OLS-style)
  # qrX <- qr(X)
  # Q <- qr.Q(qrX)
  # R <- qr.R(qrX)
  ols <- lm.fit(X, y)

  b <- ols$coefficients

  n <- nrow(X)
  n_i <- as.integer(n / nsamples)

  # aka, SSE
  S <- sum(ols$residuals ^ 2)

  # split X and y into their "within-cluster" groupings
  start <- seq(1, n, by = n_i)
  end <- start + n_i - 1

  X_i <- Map(function(s, e) X[s:e,], start, end)
  y_i <- Map(function(s, e) y[s:e], start, end)
  resid_i <- Map(function(s, e) ols$residuals[s:e], start, end)

  h_i <- sapply(seq_along(X_i),
    function(i)
    {
      mean(y_i[[i]]) - t(ols$coefficients) %*% apply(X_i[[i]], 2, mean)
    })

  A <- sum(h_i ^ 2)

  n_i2A <- (n_i^2) * A
  d <- (n_i2A - S) / ( n_i*S - n_i2A )

  # this is the (I + Z_i D Z_i')^-1 expression
  # note that in this case, it reduces to:
  # I - d/(1 + n_i *d) 1 1'
  covar_inv <- diag(n_i) - (d / (1 + n_i * d)) *
    matrix(1, nrow = n_i, ncol = n_i)

  sigma_sq <- sapply(seq_along(X_i),
    function(i)
    {
      (t(resid_i[[i]]) %*% covar_inv) %*% resid_i[[i]]
    })

  sigma_sq <- sum(sigma_sq) / (n - nsamples)

  # can simplify this when we know the design ahead of time
  cov_b <- lapply(seq_along(X_i),
    function(i)
    {
      (t(X_i[[i]]) %*% covar_inv) %*% X_i[[i]]
    })
  cov_b <- Reduce(function(x, y) x + y, cov_b) %>%
    solve()
  # XXX: cov_b is unadjusted
  # cov_b <- cov_b * sigma_sq
  fixed_sd <- sqrt( diag(cov_b) * sigma_sq )

  adjustment <- NA_real_
  if (adjust_sigma) {
    adjustment <- n / (n - length(fixed_sd))
    fixed_sd <- fixed_sd * sqrt( adjustment )
  }

  t_value <- b / fixed_sd

  log_lik <- -1/2 * (nsamples * log(1 + n_i * d) +
    nsamples * n_i * log(S - (n_i^(2) * d * A)/(1 + n_i * d)))

  res <- list(
    sigma_sq = sigma_sq,
    coef = b,
    log_lik = log_lik,
    cov_b = cov_b,
    var_b  = diag(cov_b),
    t_value = t_value,
    fixed_sd = fixed_sd,
    adjustment = adjustment
    )
  class(res) <- "sleuth_lmm"

  res
}

#' @export
compute_lmms <- function(obj, transform_fn, filter_df) {
  stopifnot( is(obj, "sleuth") )
  stopifnot( is(transform_fn, "function") )

  filter_df <- as.data.frame(filter_df)
  rownames(filter_df) <- filter_df$target_id
  filter_df$target_id <- NULL

  bs <- dcast_bootstrap(obj, "est_counts")

  # design <- lmm_design(obj)

  targs <- rownames(bs)
  filt <- filter_df[targs,1]
  stopifnot(length(targs) == length(filt))

  design <- model.matrix(~ condition, lmm_design(obj))

  bs_lmm <- lapply(1:nrow(bs),
    function(i)
    {
      tryCatch(
        {
          if (filt[i]) {
            lmm(design, transform_fn(bs[i,]), length(obj$kal))
          } else {
            stop("didn't pass filter")
          }
        },
        error = function(e) {e},
        finally = function() {}
        )
    }
    )
  names(bs_lmm) <- obj$kal[[1]]$bootstrap[[1]]$target_id

  filt_lmm <- bs_lmm %>%
    Filter(function(x) is(x, "sleuth_lmm"), .)

  obj$lmms <- filt_lmm
  obj$bs_means <- apply(bs, 1, function(row) mean(transform_fn(row)))
  obj$bs_means <- obj$bs_means[names(filt_lmm)]
  obj$sigma$bs_mean <- obj$bs_means

  reshape_lmms(obj)
}

reshape_lmms <- function(obj) {
  stopifnot( is(obj, "sleuth") )

  obj$fixed_sd <- rbind_matrix(obj$lmms, "fixed_sd")
  obj$b <- rbind_matrix(obj$lmms, "coef")
  obj$sigma <- data.frame(raw_sigma = sapply(obj$lmms,
      function(x) sqrt(x$sigma_sq)))
  obj$adjustment <- rbind_matrix(obj$lmms, "adjustment")[1L,1L]
  obj$var_b <- rbind_matrix(obj$lmms, "var_b")

  obj
}

rbind_matrix <- function(list_obj, which_slot) {
  res <- do.call(rbind,
    lapply(list_obj,
    function(x) matrix(x[[which_slot]], nrow = 1)))
  dimnames(res) <- list(names(list_obj), names(list_obj[[1]][[which_slot]]))

  res
}

#' Naive shrinkage using locfit
#'
#' The most naive shrinkage you could perform using locfit and the raw
#' variances from the LMM
#'
#' @param obj a \code{sleuth} object
#' @return a \code{sleuth} object with an extra column in the slot \code{sigma}
#' @export
naive_shrink_lmm <- function(obj) {
  stopifnot( is(obj, "sleuth") )

  # model the sqrt like limma
  smooth_model <- locfit::locfit(sqrt(raw_sigma) ~ bs_mean, data = obj$sigma)
  smooth <- locfit:::predict.locfit(smooth_model, obj$sigma$bs_mean)

  obj$sigma$naive_locfit <- smooth^2
  obj$sigma$naive_shrink_sigma <- pmax(smooth^2, obj$sigma$raw_sigma)

  obj
}


lf_pred <- function(data) {
  model <- locfit::locfit(sqrt(raw_sigma) ~ bs_mean, data = data)
  smooth_pred <- locfit:::predict.locfit(model, data$bs_mean) ^ 2
  max_pred <- pmax(smooth_pred, data$raw_sigma)

  result <- data.frame(smooth_pred = smooth_pred, max_pred = max_pred)
  rownames(result) <- rownames(data)
  result$target_id <- data$target_id

  result
}

#' Shrink lmm by grouping
#'
#' Given a data.frame with columns 'target_id' and 'shrinkage_group', fit
#' different shrinkage lines
#' @param obj a \code{sleuth} object
#' @param grouping a data.frame containing the columns 'target_id' and
#' 'shrinkage_group'
#' @return a \code{sleuth} object with an updated \code{sigma} slot with
#' updated column \code{group_shrink_sigma}
#' @export
shrink_by_group_lmm <- function(obj, grouping) {
  # TODO: check whether 'target_id' and 'shrinkage_group' are in the DF
  sigma_by_grouping <- inner_join(obj$sigma, select(grouping, target_id, shrinkage_group),
    by = "target_id")

  orig_order <- rownames(obj$sigma)

  pred <- sigma_by_grouping %>%
    group_by(shrinkage_group) %>%
    do({
      lf_pred(.)
    })

  pred <- pred %>%
    rename(group_shrink_smooth = smooth_pred, group_shrink_sigma = max_pred)

  # obj$sigma[pred$target_id, "group_shrink_smooth"] <- pred$group_shrink_smooth
  # obj$sigma[pred$target_id, "group_shrink_sigma"] <- pred$group_shrink_sigma

  obj$sigma <- inner_join(pred, obj$sigma, by = "target_id")
  obj$sigma <- as.data.frame(obj$sigma)
  rownames(obj$sigma) <- obj$sigma$target_id
  obj$sigma <- obj$sigma[orig_order,]

  obj
}

#' Compute t-statistics
#'
#' Compute t-statistics from a "sleuth" object
#'
#' @param obj a \code{sleuth} object
#' @param which_var which variance to use
#' @param adjust_sigma if TRUE, adjust sigma
#' @return a data.frame with t-statistics
#' @export
compute_t <- function(obj, which_var, adjust_sigma = TRUE) {
  stopifnot( which_var %in% colnames(obj$sigma) )

  #fixed_sd <- sqrt(obj$var_b * obj$adjustment) * matrix(rep(as.data.frame(obj$sigma)[,which_var], 2), ncol = 2)
  fixed_sd <- sqrt(obj$var_b * obj$adjustment) * obj$sigma[,which_var]
  #fixed_sd$target_id <- obj$sigma$target_id
  #stopifnot( rownames(fixed_sd) == rownames(obj$b) )

  # FIXME: figure out what DF is in general
  degrees_free <- 6
  t_stats <- obj$b / fixed_sd
  p_vals <- apply(t_stats, 2, function(col) 2 * pt(-abs(col), df = degrees_free))

  list(t_stats = t_stats, p_vals = p_vals)
}

#' OLS on observations
#'
#' @param obj a \code{sleuth} object
#' @param design a matrix
#' @return a \code{sleuth} object with an updated slot containing slot "obs_beta"
#' @export
ols_by_row <- function(obj, design, transform = identity) {
  stopifnot( is(obj, "sleuth") )

  obs_counts <- reshape2::dcast(obj$obs_norm, target_id ~ sample,
    value.var = "est_counts")
  obs_counts <- as.data.frame(obs_counts)
  rownames(obs_counts) <- obs_counts$target_id
  obs_counts$target_id <- NULL
  obs_counts <- as.matrix(obs_counts)

  obs_counts <- transform(obs_counts)

  all_fits <- lapply(1:nrow(obs_counts),
    function(it)
    {
      matrix(lm.fit(design, obs_counts[it,])$coefficients, nrow = 1)
    })

  coef_matrix <- do.call("rbind", all_fits)
  dimnames(coef_matrix) <- list(rownames(obs_counts), colnames(design))

  coef_matrix
}
