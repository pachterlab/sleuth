#' @export
get_quantile <- function(data, which_col, lwr, upr, ignore_zeroes = TRUE) {
  data <- as.data.frame(data)

  valid <- rep(TRUE, nrow(data))

  data$iqr <- valid

  if (ignore_zeroes) {
    valid <- data[,which_col] > 0
    if (sum(valid) == 0) {
      # in this case, everything in this bin is zero, so we just return the
      # entire bin
      return(data)
    }

    # replace all the invalid points with false
    data$iqr[!valid] <- FALSE
  }

  data$iqr[valid] <- data[valid,which_col] %>%
    ecdf(.)(.) %>%
    ifelse(lwr <= . & . <= upr)

  data
}

#' @export
sliding_window_grouping <- function(data, x_col, y_col,
  n_bins = 100, lwr = 0.25, upr = 0.75, ignore_zeroes = TRUE) {
  data %>%
    mutate(x_ecdf = ecdf(data[,x_col])(data[,x_col])) %>%
    mutate(x_group = cut(x_ecdf, n_bins)) %>%
    group_by(x_group) %>%
    do(get_quantile(., y_col, lwr, upr, ignore_zeroes)) %>%
    select(-c(x_ecdf, x_group))
}

#' @export
shrink_df <- function(data, shrink_formula, filter_var) {
  data <- as.data.frame(data)
  s_formula <- substitute(shrink_formula)
  fit <- eval(loess(s_formula, data[data[,filter_var],]))
  predict(fit, data)
}
