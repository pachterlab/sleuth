#
#    sleuth: inspect your RNA-Seq with a pack of kallistos
#
#    Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

  data <- as.data.frame(data)

  data <- mutate(data,x_ecdf = ecdf(data[,x_col])(data[,x_col]))
  data <- mutate(data, x_group = cut(x_ecdf, n_bins))
  data <- group_by(data, x_group)

  res <- do(data, get_quantile(., y_col, lwr, upr, ignore_zeroes))
  res <- select(res, -c(x_ecdf, x_group))

  ungroup(res)
}

#' @export
shrink_df <- function(data, shrink_formula, filter_var) {
  data <- as.data.frame(data)
  s_formula <- substitute(shrink_formula)
  fit <- eval(loess(s_formula, data[data[,filter_var],]))
  data.frame(data, shrink = predict(fit, data))
}

msg <- function(..., nl = TRUE) {
  message(..., appendLF = nl)
}

# shortcut for as.data.frame(x, stringsAsFactors = FALSE)
as_df <- function(x, ...) {
  as.data.frame(x, stringsAsFactors = FALSE, ...)
}

adf <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}

dot <- function(count, max_dots = 50) {
  count <- count + 1
  new_line <- FALSE
  if ( (count %% max_dots) == 0 ) {
    new_line <- TRUE
  }

  msg('.', nl = new_line)

  count
}

kld <- function(p, q) {
  stopifnot(length(p) == length(q))

  which_valid <- which(p > 0 & q > 0)

  p <- p[which_valid]
  q <- q[which_valid]

  sum(p * (log(p) - log(q)))
}

jsd <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)

  m <- (p + q)/2
  (kld(p, m) + kld(q, m)) / 2
}

apply_all_pairs <- function(mat, fun) {
  ids <- colnames(mat)

  res <- matrix(NA, nrow = length(ids), ncol = length(ids))
  dimnames(res) <- list(ids, ids)

  all_pairs <- utils::combn(ids, 2)
  for (i in 1:ncol(all_pairs)) {
    j <- all_pairs[1,i]
    k <- all_pairs[2,i]

    cur <- fun(mat[,j], mat[,k])

    res[j,k] <- cur
    res[k,j] <- cur
  }

  for (i in 1:length(ids)) {
    res[i,i] <- fun(mat[,i], mat[,i])
  }


  res
}
