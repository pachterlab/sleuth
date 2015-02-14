#' Read a kallisto data set
#'
#' Read a kallisto data set
#'
#' @param output_dir the directory of the output data
#' @param read_bootstrap if TRUE, then searches for bootstrap data, else doesn't read it.
#' @return a S3 \code{kallisto} object with the following members:
#' @export
read_kallisto <- function(output_dir, read_bootstrap = TRUE)
{
    if (!file.exists(output_dir) || !file.info(output_dir)$isdir) {
        stop(paste0("'", output_dir, "' is not a valid path"))
    }

    cat("Reading main abundance estimates\n")
    exp_fname <- file.path(output_dir, "expression.txt")
    if (!file.exists( exp_fname )) {
        stop("'", exp_fname, "' does not exists. Are you sure you have the right kallisto output directory?")
    }

    trans_abund <- suppressWarnings(fread(exp_fname, data.table = FALSE)) %>%
        arrange(target_id)

    bs_samples <- NULL

    if (read_bootstrap) {
        bs_fnames <- Sys.glob(file.path(output_dir, "bs_expression_*"))
        if (length(bs_fnames) > 0) {
            cat("Found",length(bs_fnames), "bootstrap files. Reading them in.\n")
            bs_samples <- lapply(seq_along(bs_fnames), function(b)
                {
                    cat(".")
                    if (b %% 50 == 0 && b > 0) {
                        cat("\n")
                    }
                    fname <- bs_fnames[b]
                    suppressWarnings(fread(fname, data.table = FALSE)) %>%
                        arrange(target_id)
                })
        } else {
            warning("No bootstrap samples found!")
        }
    }
    # XXX: should we check that all targets are identical?

    invisible(structure(list(abundance = trans_abund, bootstrap = bs_samples),
        class = "kallisto"))
}

