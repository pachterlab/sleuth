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

    tidyr::gather_(boot, "sample", column, bs_names)
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
