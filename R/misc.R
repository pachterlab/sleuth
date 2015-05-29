#' @export
get_quantile <- function(data, which_col, lwr, upr) {
  data <- as.data.frame(data)
  data$cdf <- data[,which_col] %>%
    ecdf(.)(.) %>%
    ifelse(lwr <= . & . <= upr)
  data
}
