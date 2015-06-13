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
#' @export
lmm <- function(X, y, nsamples) {
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

  h_i <- lapply(seq_along(X_i),
    function(i)
    {
      mean(y_i[[i]]) - t(ols$coefficients) %*% apply(X_i[[i]], 2, mean)
    }) %>%
      unlist()

  A <- sum(h_i ^ 2)

  n_i2A <- (n_i^2) * A
  d <- (n_i2A - S) / ( n_i*S - n_i2A )

  # this is the (I + Z_i D Z_i')^-1 expression
  # note that in this case, it reduces to:
  # I - d/(1 + n_i *d) 1 1'
  covar_inv <- diag(n_i) - (d / (1 + n_i * d)) *
    matrix(1, nrow = n_i, ncol = n_i)

  sigma_sq <- lapply(seq_along(X_i),
    function(i)
    {
      (t(resid_i[[i]]) %*% covar_inv) %*% resid_i[[i]]
    })

  sigma_sq <- sum(unlist(sigma_sq)) / (n - nsamples)

  # can simplify this when we know the design ahead of time
  cov_b <- lapply(seq_along(X_i),
    function(i)
    {
      (t(X_i[[i]]) %*% covar_inv) %*% X_i[[i]]
    })
  cov_b <- Reduce(function(x, y) x + y, cov_b) %>%
    solve()
  cov_b <- cov_b * sigma_sq

  log_lik <- -1/2 * (nsamples * log(1 + n_i * d) +
    nsamples * n_i * log(S - (n_i^(2) * d * A)/(1 + n_i * d)))

  list(sigma_sq = sigma_sq, coef = b, log_lik = log_lik, cov_b = cov_b)
}
