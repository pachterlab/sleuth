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

#' Fit a measurement error model
#'
#' This function is a wrapper for fitting a measurement error model using
#' \code{sleuth}. It performs the technical variance estimation from the boostraps, biological
#' variance estimation, and shrinkage estimation.
#'
#' For most users, simply providing the sleuth object should be sufficient. By
#' default, this behavior will fit the full model initially specified and store
#' it in the sleuth object under 'full'.
#'
#' To see which models have been fit, users will likely find the function
#' \code{\link{models}} helpful.
#'
#' @param obj a \code{sleuth} object
#' @param formula an R formula specifying the design to fit OR a design matrix.
#' If you are interested in only fitting the model that was specified in \code{sleuth_prep}
#' you do not need to specify it again (will be fit as the 'full' model).
#' @param fit_name the name to store the fit in the sleuth object (at so$fits$fit_name).
#' If \code{NULL}, the model will be named 'full'.
#' @param ... additional arguments passed to \code{sliding_window_grouping} and
#' \code{shrink_df}
#' @return a sleuth object with updated attributes.
#' @examples # If you specified the formula in sleuth_prep, you can simply run to run the full model
#' so <- sleuth_fit(so)
#' # The intercept only model can be fit like this
#' so <- sleuth_fit(so, ~1, 'reduced')
#' @seealso \code{\link{models}} for seeing which models have been fit,
#' \code{\link{sleuth_prep}} for creating a sleuth object,
#' \code{\link{sleuth_wt}} to test whether a coefficient is zero,
#' \code{\link{sleuth_lrt}} to test nested models.
#' @export
sleuth_fit <- function(obj, formula = NULL, fit_name = NULL, ...) {
  stopifnot( is(obj, 'sleuth') )

  if ( is.null(formula) ) {
    formula <- obj$full_formula
    if (is.null(formula)) {
      stop("'formula' was not specified and the 'full' model was not specified in `sleuth_prep`.",
        " Please specify a formula and a label.")
    }
  } else if ( !is(formula, 'formula') && !is(formula, 'matrix') ) {
    stop("'", substitute(formula), "' is not a valid 'formula' or 'matrix'")
  }

  if ( is.null(fit_name) ) {
    fit_name <- 'full'
  } else if ( !is(fit_name, 'character') ) {
    stop("'", substitute(fit_name), "' is not a valid 'character'")
  }

  if ( length(fit_name) > 1 ) {
    stop("'", substitute(fit_name), "' is of length greater than one.",
      " Please only supply one string.")
  }

  # TODO: check if model matrix is full rank
  X <- NULL
  if ( is(formula, 'formula') ) {
    X <- model.matrix(formula, obj$sample_to_covariates)
  } else {
    if ( is.null(colnames(formula)) ) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    X <- formula
  }
  rownames(X) <- obj$sample_to_covariates$sample
  A <- solve(t(X) %*% X)

  msg("fitting measurement error models")
  mes <- me_model_by_row(obj, X, obj$bs_summary)
  tid <- names(mes)

  mes_df <- dplyr::bind_rows(lapply(mes,
    function(x) {
      data.frame(rss = x$rss, sigma_sq = x$sigma_sq, sigma_q_sq = x$sigma_q_sq,
        mean_obs = x$mean_obs, var_obs = x$var_obs)
    }))

  mes_df$target_id <- tid
  rm(tid)

  mes_df <- dplyr::mutate(mes_df, sigma_sq_pmax = pmax(sigma_sq, 0))

  # FIXME: sometimes when sigma is negative the shrinkage estimation becomes NA
  # this is for the few set of transcripts, but should be able to just do some
  # simple fix
  msg('shrinkage estimation')
  swg <- sliding_window_grouping(mes_df, 'mean_obs', 'sigma_sq_pmax',
    ignore_zeroes = TRUE, ...)

  l_smooth <- shrink_df(swg, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, 'iqr', ...)
  l_smooth <- dplyr::select(
    dplyr::mutate(l_smooth, smooth_sigma_sq = shrink ^ 4),
    -shrink)

  l_smooth <- dplyr::mutate(l_smooth,
    smooth_sigma_sq_pmax = pmax(smooth_sigma_sq, sigma_sq))

  msg('computing variance of betas')
  beta_covars <- lapply(1:nrow(l_smooth),
    function(i) {
      row <- l_smooth[i, ]
      with(row,
          covar_beta(smooth_sigma_sq_pmax + sigma_q_sq, X, A)
        )
    })
  names(beta_covars) <- l_smooth$target_id

  if ( is.null(obj$fits) ) {
    obj$fits <- list()
  }

  obj$fits[[fit_name]] <- list(
    models = mes,
    summary = l_smooth,
    beta_covars = beta_covars,
    formula = formula,
    design_matrix = X,
    transform_synced = TRUE)

  class(obj$fits[[fit_name]]) <- 'sleuth_model'

  obj
}

# if 'fail' is set to TRUE, and fail if the model is not found and give an
# error message
model_exists <- function(obj, which_model, fail = TRUE) {
  stopifnot( is(obj, 'sleuth') )
  stopifnot( is(which_model, 'character') )
  stopifnot( length(which_model) == 1 )

  result <- which_model %in% names(obj$fits)
  if (fail && !result) {
    stop("model '", which_model, "' not found")
  }

  result
}

#' Wald test for a sleuth model
#'
#' This function computes the Wald test on one specific 'beta' coefficient on
#' every transcript.
#'
#' @param obj a \code{sleuth} object
#' @param which_beta a character string of denoting which grouping to test.
#' For example, if you have a model fit to 'treatment,' with values of neg_ctl, pos_ctl,
#' and drug, you would need to run \code{sleuth_wt} once each for pos_ctl and drug
#' @param which_model a character string of length one denoting which model to
#' use
#' @return an updated sleuth object
#' @examples # Assume we have a sleuth object with a model fit to both genotype and drug,
#' models(so)
#' # formula:  ~genotype + drug
#' # coefficients:
#' #   (Intercept)
#' #   genotypeKO
#' #   drugDMSO
#' so <- sleuth_wt(so, 'genotypeKO')
#' so <- sleuth_wt(so, 'drugDMSO')
#' @seealso \code{\link{models}} to view which models have been fit and which
#' coefficients can be tested, \code{\link{sleuth_results}} to get back
#' a \code{data.frame} of the results
#' @export
sleuth_wt <- function(obj, which_beta, which_model = 'full') {
  stopifnot( is(obj, 'sleuth') )

  if ( !model_exists(obj, which_model) ) {
    stop("'", which_model, "' is not a valid model. Please see models(",
      substitute(obj), ") for a list of fitted models")
  }

  if(!obj$fits[[which_model]]$transform_synced) {
    stop("Model '", which_model, "' was not computed using the sleuth object's",
         " current transform function. Please rerun sleuth_fit for this model.")
  }

  d_matrix <- obj$fits[[which_model]]$design_matrix

  # get the beta index
  beta_i <- which(colnames(d_matrix) %in% which_beta)

  if ( length(beta_i) == 0 ) {
    stop(paste0("'", which_beta,
        "' doesn't appear in your design. Try one of the following:\n",
        paste(colnames(d_matrix), collapse = ' ')))
  } else if ( length(beta_i) > 1 ) {
    stop(paste0("Sorry. '", which_beta, "' is ambiguous for columns: ",
        paste(colnames(d_matrix[beta_i]), collapse = ' ')))
  }

  b <- sapply(obj$fits[[ which_model ]]$models,
    function(x) {
      x$ols_fit$coefficients[ beta_i ]
    })
  names(b) <- names(obj$fits[[ which_model ]]$models)

  res <- obj$fits[[ which_model ]]$summary
  res$target_id <- as.character(res$target_id)
  res <- res[match(names(b), res$target_id), ]

  stopifnot( all.equal(res$target_id, names(b)) )

  se <- sapply(obj$fits[[ which_model ]]$beta_covars,
    function(x) {
      x[beta_i, beta_i]
    })
  se <- sqrt( se )
  se <- se[ names(b) ]

  stopifnot( all.equal(names(b), names(se)) )

  res <- dplyr::mutate(res,
    b = b,
    se_b = se,
    wald_stat = b / se,
    pval = 2 * pnorm(abs(wald_stat), lower.tail = FALSE),
    qval = p.adjust(pval, method = 'BH')
    )

  obj <- add_test(obj, res, which_beta, 'wt', which_model)

  obj
}

# Compute the covariance on beta under OLS
#
# Compute the covariance on beta under OLS
# @param sigma a numeric of either length 1 or nrow(X) defining the variance
# on D_i
# @param X the design matrix
# @param A inv(t(X) X) (for speedup)
# @return a covariance matrix on beta
covar_beta <- function(sigma, X, A) {
  if (length(sigma) == 1) {
    return( sigma * A )
  }

  # sammich!
  A %*% (t(X) %*% diag(sigma) %*% X) %*% A
}

# Measurement error model
#
# Fit the measurement error model across all samples
#
# @param obj a \code{sleuth} object
# @param design a design matrix
# @param bs_summary a list from \code{bs_sigma_summary}
# @return a list with a bunch of objects that are useful for shrinking
me_model_by_row <- function(obj, design, bs_summary) {
  stopifnot( all.equal(names(bs_summary$sigma_q_sq), rownames(bs_summary$obs_counts)) )
  stopifnot( length(bs_summary$sigma_q_sq) == nrow(bs_summary$obs_counts))

  models <- lapply(1:nrow(bs_summary$obs_counts),
    function(i) {
      me_model(design, bs_summary$obs_counts[i, ], bs_summary$sigma_q_sq[i])
    })
  names(models) <- rownames(bs_summary$obs_counts)

  models
}

# non-equal var
#
# word
#
# @param obj a sleuth object
# @param design a design matrix
# @param samp_bs_summary the sample boostrap summary computed by sleuth_summarize_bootstrap_col
# @return a list with a bunch of objects used for shrinkage :)
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
      res <- me_white_model(design, obs_counts[i, ], sigma_q_sq[i, ], A)
      res$df$target_id <- rownames(obs_counts)[i]
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
bs_sigma_summary <- function(obj, transform = identity, norm_by_length = FALSE) {
  # if (norm_by_length) {
  #   scaling_factor <- get_scaling_factors(obj$obs_raw)
  #   reads_per_base_transform()
  # }
  obs_counts <- obs_to_matrix(obj, "est_counts")
  obs_counts <- transform( obs_counts )

  bs_summary <- sleuth_summarize_bootstrap_col(obj, "est_counts", transform)
  bs_summary <- dplyr::group_by(bs_summary, target_id)
  bs_summary <- dplyr::summarise(bs_summary,
    sigma_q_sq = mean(bs_var_est_counts))

  bs_summary <- as_df(bs_summary)

  bs_sigma <- bs_summary$sigma_q_sq
  names(bs_sigma) <- bs_summary$target_id
  bs_sigma <- bs_sigma[rownames(obs_counts)]

  list(obs_counts = obs_counts, sigma_q_sq = bs_sigma)
}

# transform reads into reads per base
#
#
reads_per_base_transform <- function(reads_table, scale_factor_input,
  collapse_column = NULL,
  mapping = NULL,
  norm_by_length = TRUE) {

  reads_table <- data.table::as.data.table(reads_table)

  if (is(scale_factor_input, 'data.frame')) {
    scale_factor_input <- data.table::as.data.table(dplyr::select(scale_factor_input, target_id,
                                                                  sample, scale_factor))
    reads_table <- merge(reads_table, scale_factor_input,
      by = c('sample', 'target_id'), all.x=T)
  } else {
    reads_table[, scale_factor := scale_factor_input]
  }
  # browser()
  reads_table[, reads_per_base := est_counts / eff_len]
  reads_table[, scaled_reads_per_base := scale_factor * reads_per_base]

  if (!is.null(collapse_column)) {
    mapping <- data.table::as.data.table(mapping)
    # old stuff
    if (!(collapse_column %in% colnames(reads_table))) {
      reads_table <- merge(reads_table, mapping, by = 'target_id', all.x=T)
    }
    # browser()
    # reads_table <- dplyr::left_join(reads_table, mapping, by = 'target_id')

    rows_to_remove <- !is.na(reads_table[[collapse_column]])
    reads_table <- reads_table[rows_to_remove]
    if ('sample' %in% colnames(reads_table)) {
      reads_table <- reads_table[, j = list(scaled_reads_per_base = sum(scaled_reads_per_base)),
                  by = list(sample, eval(parse(text=collapse_column)))]
    } else {
      reads_table <- reads_table[, j = list(scaled_reads_per_base = sum(scaled_reads_per_base)),
                  by = eval(parse(text=collapse_column))]
    }

    data.table::setnames(reads_table, 'parse', 'target_id')
  }

  as_df(reads_table)
}

gene_summary <- function(obj, which_column, transform = identity,
                         norm_by_length = TRUE, num_cores=2) {
  # stopifnot(is(obj, 'sleuth'))
  msg(paste0('aggregating by column: ', which_column))
  apply_function <- if (num_cores == 1) {
    lapply
  } else {
    function(x, y) parallel::mclapply(x, y, mc.cores = num_cores)
  }
  obj_mod <- obj
  if (norm_by_length) {
    tmp <- obj$obs_raw
    tmp <- dplyr::left_join(
      data.table::as.data.table(tmp),
      data.table::as.data.table(obj$target_mapping),
      by = 'target_id')
    tmp <- dplyr::group_by_(tmp, 'sample', which_column)
    scale_factor <- dplyr::mutate(tmp, scale_factor = median(eff_len))
  } else {
    scale_factor <- median(obj_mod$obs_norm_filt$eff_len)
  }
  obj_mod$obs_norm_filt <- reads_per_base_transform(obj_mod$obs_norm_filt,
    scale_factor, which_column, obj$target_mapping, norm_by_length)
  obj_mod$obs_norm <- reads_per_base_transform(obj_mod$obs_norm,
    scale_factor, which_column, obj$target_mapping, norm_by_length)

  obs_counts <- obs_to_matrix(obj_mod, "scaled_reads_per_base")
  obs_counts <- transform(obs_counts)

  msg("starting mclapply process now")
  obj_mod$kal <- lapply(seq_along(obj_mod$kal),
    function(i) {
      k <- obj_mod$kal[[i]]
      current_sample <- obj_mod$sample_to_covariates$sample[i]
      msg(paste('aggregating across sample: ', current_sample))
      k$bootstrap <- apply_function(seq_along(k$bootstrap), function(j) {
        b <- k$bootstrap[[j]]
        b <- dplyr::mutate(b, sample = current_sample)
        reads_per_base_transform(b, scale_factor, which_column,
          obj$target_mapping, norm_by_length)
      })

      k
    })

  bs_summary <- sleuth_summarize_bootstrap_col(obj_mod, "scaled_reads_per_base",
    transform)

  bs_summary <- dplyr::group_by(bs_summary, target_id)
  # FIXME: the column name 'bs_var_est_counts' is incorrect. should actually rename it above
  bs_summary <- dplyr::summarise(bs_summary,
    sigma_q_sq = mean(bs_var_scaled_reads_per_base))

  bs_summary <- as_df(bs_summary)

  bs_sigma <- bs_summary$sigma_q_sq
  names(bs_sigma) <- bs_summary$target_id
  bs_sigma <- bs_sigma[rownames(obs_counts)]

  list(obs_counts = obs_counts, sigma_q_sq = bs_sigma)
}

me_model <- function(X, y, sigma_q_sq) {
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
