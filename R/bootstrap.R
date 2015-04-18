#' Bootstrap to matrix
#'
#' Takes a \code{kallisto} object and converts the bootstrap results into
#' a proper \code{matrix}
#'
#' @param kal a kallisto object with non-null member \code{bootstrap}
#' @return a matrix with rownames equal to target_id
#' @export
bootstrap2mat <- function(kal, column = "tpm")
{
    stopifnot(is(kal, "kallisto"))
    # TODO: check if "column" is a valid kallisto column

    # assumes that all bootstrap samples are in same order (from read_kallisto)

    all_boot <- kal$bootstrap
    mat <- matrix(unlist(lapply(all_boot, function(bs)
        {
            bs[column]
        })), nrow = nrow(all_boot[[1]]))

    rownames(mat) <- all_boot[[1]]$target_id

    mat
}


#' Convert kallisto bootstraps into a molten data.frame
#'
#' Melt it!
#'
#' @param kal a kallisto object
#' @param column the column to pull out of the kallisto results (default = "tpm")
#' @return a molten data.frame with columns "target_id", "sample" and the selected variable
#' @export
melt_bootstrap <- function(kal, column = "tpm")
{
    stopifnot(is(kal, "kallisto"))
  stopifnot(length(kal$bootstrap) > 0)

    all_boot <- kal$bootstrap
    boot <- data.frame(lapply(all_boot, select_, .dots = list(column)))
    bs_names <- paste0("bs", 1:ncol(boot))
    data.table::setnames(boot, colnames(boot), bs_names)
    boot <- boot %>%
        mutate(target_id = all_boot[[1]]$target_id)

    tidyr::gather_(boot, "sample", column, bs_names) %>%
      mutate(sample = as.factor(sample))
}

#' Summarize bootstrap values
#'
#' Compute the mean, sd, var, and coefficient of variation from a kallisto
#' bootstrap
#' @param kal a kallisto object with a non-null bootstrap list
#' @param column the column to select (rho, tpm, est_counts
#' @return a summarized data.frame
#' @export
summarize_bootstrap <- function(kal, column = "tpm")
{
    stopifnot(is(kal, "kallisto"))
    bs <- melt_bootstrap(kal, column)

    mean_col <- paste0("bs_mean_", column)
    sd_col <- paste0("bs_sd_", column)
    var_col <- paste0("bs_var_", column)
    cv_col <- paste0("bs_cv_", column)

    print(cv_col)

    bs <- bs %>%
        group_by(target_id) %>%
        summarise_(.dots = setNames(list(
                    interp(quote(mean(x)), x = as.name(column)),
                    interp(quote(sd(x)), x = as.name(column)),
                    interp(quote(var(x)), x = as.name(column))
                    ),
                c(mean_col, sd_col, var_col)))

    bs <- bs %>%
        mutate_(.dots = setNames(list(
                    interp(quote(x / y),
                        x = as.name(sd_col), y = as.name(mean_col))),
                c(cv_col)
                ))
    bs
}

#' Normalize bootstrap samples
#'
#' Normalize by dividing by the "size factor" for each TPM and estimated counts
#'
#' @param kal a kallisto object
#' @param tpm_size_factor the size factor (numeric length 1)
#' @param est_counts_size_factor the size factor (numeric length 1)
#' @export
normalize_bootstrap <- function(kal, tpm_size_factor, est_counts_size_factor) {
  stopifnot(is(kal, "kallisto"))
  stopifnot(length(tpm_size_factor) == 1 && length(est_counts_size_factor) == 1)

  calc_norm_tpm <- !missing(tpm_size_factor)
  calc_norm_counts <- !missing(est_counts_size_factor)
  bs <- lapply(kal$bootstrap, function(bs_tbl)
    {
      if (calc_norm_tpm)
        bs_tbl$tpm <- bs_tbl$tpm / tpm_size_factor
      if (calc_norm_counts)
        bs_tbl$est_counts <- bs_tbl$est_counts / est_counts_size_factor

      bs_tbl
    })
  kal$bootstrap <- bs

  kal
}


#' Sample bootstraps
#'
#' From a sleuth object, create experiments by randomly sampling bootstraps from each kallisto object
#' @param obj a \code{kallisto} object
#' @param n_samples the number of samples to genenerate
#' @export
sample_bootstrap <- function(obj, n_samples = 100L) {
  stopifnot( is(obj, "sleuth") )

  n_kal <- length(obj$kal)
  n_bs_per_samp <- unlist(lapply(obj$kal, function(x) length(x$bootstrap)))
  if (any(n_bs_per_samp < n_samples)) {
    warning("You've asked to sample more samples than you have bootstraps.",
      " We recommend you generate more bootstrap samples in kallisto...")
  }

  which_samp <- lapply(seq_along(n_bs_per_samp),
    function(i)
    {
      cur_n <- n_bs_per_samp[i]
      sample.int(cur_n, n_samples, replace = TRUE)
    })
  # each column contains which bootstrap sample we want from each kallisto
  which_samp <- t(simplify2array(which_samp))

  # allocate the matrices
  sample_mat <- lapply(1:n_samples,
    function(discard)
    {
      mat <- matrix(NA_real_, nrow = nrow(obj$kal[[1]]$abundance),
        ncol = nrow(which_samp))
      rownames(mat) <- obj$kal[[1]]$abundance$target_id
      colnames(mat) <- obj$sample_to_condition$sample
      mat
    })

  # go through the generated sample ids and place them into the corresponding
  # matrix sample
  for (s in 1:n_samples) {
    for (idx in 1:nrow(which_samp)) {
      b <- which_samp[idx,s]
      sample_mat[[s]][,idx] <- obj$kal[[idx]]$bootstrap[[b]]$est_counts
    }
  }

  sample_mat
}

#' @export
dcast_bootstrap <- function(obj, ...) {
  UseMethod("dcast_bootstrap")
}

#' @export
dcast_bootstrap.sleuth <- function(obj, units) {
  bs <- lapply(obj[["kal"]], dcast_bootstrap, units)

  do.call(cbind, bs)
}

#' @export
dcast_bootstrap.kallisto <- function(obj, units) {
  if ( !(units %in% c("est_counts", "tpm")) ) {
    stop(paste0("'", substitute(units),
        "' is not valid for 'units'. Please see documentation"))
  }

  if ( length(obj$bootstrap) < 1 ) {
    stop("No bootstrap samples found.")
  }

  n_bs <- length(obj$bootstrap)
  n_features <- nrow(obj$bootstrap[[1]])
  mat <- matrix(NA_real_, nrow = n_features, ncol = n_bs)

  for (j in seq_along(obj$bootstrap)) {
    mat[ ,j] <- obj[[ "bootstrap" ]][[j]][[ units ]]
  }
  rownames(mat) <- obj[["bootstrap"]][[1]][["target_id"]]

  mat
}
