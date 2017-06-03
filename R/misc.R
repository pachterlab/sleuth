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
    valid <- data[, which_col] > 0
    if (sum(valid) == 0) {
      # in this case, everything in this bin is zero, so we just return the
      # entire bin
      return(data)
    }

    # replace all the invalid points with false
    data$iqr[!valid] <- FALSE
  }

  data$iqr[valid] <- data[valid, which_col] %>%
    ecdf(.)(.) %>%
    ifelse(lwr <= . & . <= upr)

  data
}

#' @export
sliding_window_grouping <- function(data, x_col, y_col,
  n_bins = 100, lwr = 0.25, upr = 0.75, ignore_zeroes = TRUE) {

  data <- as.data.frame(data)

  data <- mutate(data, x_ecdf = ecdf(data[, x_col])(data[, x_col]))
  data <- mutate(data, x_group = cut(x_ecdf, n_bins))
  data <- group_by(data, x_group)

  res <- do(data, get_quantile(., y_col, lwr, upr, ignore_zeroes))
  res <- ungroup(res)
  res <- select(res, -c(x_ecdf, x_group))

  res
}

#' @export
shrink_df <- function(data, shrink_formula, filter_var) {
  data <- as.data.frame(data)
  s_formula <- substitute(shrink_formula)
  fit <- eval(loess(s_formula, data[data[, filter_var], ], model = T))
  shrink_data <- data.frame(data, shrink = predict(fit, data))

  # include a column to report any target IDs that failed shrinkage estimation
  # 'failed_ise' is short for 'failed initial shrinkage estimation'
  shrink_data$failed_ise <- is.na(shrink_data$shrink)
  if (any(shrink_data$failed_ise)) {
    na_rows <- which(shrink_data$failed_ise)
    num_na <- length(na_rows)
    which_na <- shrink_data$target_id[na_rows]
    max_mean_obs <- summary(fit$model$mean_obs)["Max."]
    min_mean_obs <- summary(fit$model$mean_obs)["Min."]
    # check for mean_obs values that fall outside of the window that
    # was used to fit the LOESS curve. If any values do fall outside,
    # repeat the LOESS fit using "surface = direct" to get an exact surface fit,
    # which can then be used for extrapolation of these values.
    if (any(shrink_data$mean_obs > max_mean_obs |
            shrink_data$mean_obs < min_mean_obs)) {
      message(num_na, " NA values were found during variance shrinkage estimation",
              " due to mean observation values outside of the range used for the LOESS fit.\n",
              "The LOESS fit will be repeated using exact computation of the fitted ",
              "surface to extrapolate the missing values.\n",
              "These are the target ids with NA values: ", paste(which_na))
      direct_fit <- eval(loess(s_formula, data[data[, filter_var], ], model = T,
                  control = loess.control(surface = "direct")))
      na_data <- data[na_rows, ]
      na_shrink_data <- data.frame(na_data, shrink = predict(direct_fit, na_data))
      shrink_data[na_rows, "shrink"] <- na_shrink_data$shrink
    }
    # repeat the check for NA values; if there were no values outside of the mean_obs window used
    # (first check above) or if there were new NA values from the direct fit, then
    # the cause of the NA values is unknown. The user will be warned that this will cause
    # all downstream testing to fail for these target IDs.
    if (any(is.na(shrink_data$shrink))) {
      stop(num_na, " NA values were found during variance shrinkage estimation",
              " due to an unknown cause. These values will result in NAs ",
              "with any downstream testing for these target ids.\n",
              "Please submit a bug report at the Sleuth Github website.\n",
              "These are the target ids with NA values: ", paste(which_na))
    }
  }
  shrink_data
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

  m <- (p + q) / 2
  (kld(p, m) + kld(q, m)) / 2
}

apply_all_pairs <- function(mat, fun) {
  ids <- colnames(mat)

  res <- matrix(NA, nrow = length(ids), ncol = length(ids))
  dimnames(res) <- list(ids, ids)

  all_pairs <- utils::combn(ids, 2)
  for (i in 1:ncol(all_pairs)) {
    j <- all_pairs[1, i]
    k <- all_pairs[2, i]

    cur <- fun(mat[, j], mat[, k])

    res[j, k] <- cur
    res[k, j] <- cur
  }

  for (i in 1:length(ids)) {
    res[i, i] <- fun(mat[, i], mat[, i])
  }


  res
}

# Check if the number of cores is a sane number
# Also check for situations where cores needs to be set because
# of the user environment
check_num_cores <- function(cores) {
  if (is.null(cores) || is.na(suppressWarnings(as.integer(cores))) ||
       cores < 1 || cores > parallel::detectCores()) {
    warning("'num_cores' must be an integer between 1 and the number ",
            "of cores on your machine. Since the value given was '", num_cores,
            "', 'num_cores' has been set to 1 so that your analysis can run.")
    cores <- 1
  }

  # The technique used here to test if the user is running sleuth from
  # RStudio was taken from the following stackoverflow thread:
  # https://stackoverflow.com/a/17804414
  if(cores > 1 && Sys.getenv("RSTUDIO") == "1") {
    warning("It appears that you are running Sleuth from within Rstudio.\n",
            "Because of concerns with forking processes from a GUI, ",
            "'num_cores' is being set to 1.\nIf you wish to take ",
            "advantage of multiple cores, please consider running ",
            "sleuth from the command line.")
    cores <- 1
  }

  cores
}
