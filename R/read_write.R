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

    exp_fname <- file.path(output_dir, "expression.txt")
    if (!file.exists( exp_fname )) {
        stop("'", exp_fname, "' does not exists. Are you sure you have the right kallisto output directory?")
    }

    trans_abund <- suppressWarnings(fread(exp_fname, data.table = FALSE))

    bs_samples <- NULL

    if (read_bootstrap) {
        bs_fnames <- Sys.glob(file.path(output_dir, "bs_expression_*"))
        if (length(bs_fnames) > 0) {
            bs_samples <- lapply(bs_fnames, function(fname)
                {
                    suppressWarnings(fread(fname, data.table = FALSE))
                })
        } else {
            warning("No bootstrap samples found!")
        }
    }

    invisible(structure(list(abundance = trans_abund, bootstrap = bs_samples),
        class = "kallisto"))
}
