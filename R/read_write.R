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

    bs_samples <- list()

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
                    suppressWarnings(data.table::fread(fname, data.table = FALSE)) %>%
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


#' @export
kv_vec_to_df <- function(x, cols = c("gene_id", "transcript_id")) {
  stopifnot(length(x) %% 2 == 0)

  key_idx <- seq(1, length(x), 2)
  val_idx <- key_idx + 1

  vals <- x[val_idx]
  vals <- vals %>%
    gsub('"', "", .) %>%
    gsub(";", "", .)

  res <- setNames(as.list(vals), x[key_idx])
  res[cols]
}

#' @export
gtf_attributes_to_gene_trans <- function(gtf_attr) {
  stopifnot(is(gtf_attr, "character"))

  lapply(strsplit(gtf_attr, " "), kv_vec_to_df) %>%
    rbind_all()
}

#' @export
gtf_gene_names <- function(gtf_attr) {
  stopifnot(is(gtf_attr, "character"))
  all_attr <- strsplit(gtf_attr, " ")

  gene_id <- vector("character", length(all_attr))
  trans_id <- vector("character", length(all_attr))

  for (i in 1:length(all_attr)) {
    j <- 1
    while ((nchar(gene_id[i]) < 1 || nchar(trans_id[i]) < 1) &&
      j <= length(all_attr[[i]]) ) {
      if (all_attr[[i]][j] == "gene_id") {
        gene_id[i] <- all_attr[[i]][j+1] %>%
          gsub('"', "", .) %>%
          sub(";", "", .)
        j <- j + 2
      } else if (all_attr[[i]][j] == "transcript_id") {
        trans_id[i] <- all_attr[[i]][j+1] %>%
          gsub('"', "", .) %>%
          sub(";", "", .)
        j <- j + 2
      } else {
        j <- j + 1
      }
    }
  }

  data.frame(gene_id = gene_id, transcript_id = trans_id)
}

#' @export
read_gtf <- function(fname) {
  gtf <- data.table::fread(fname, sep = "\t", header = FALSE,
    data.table = FALSE)

  gtf_colnames <- c("seqname", "source", "feature", "start", "end", "score",
    "strand", "frame", "attribute")

  data.table::setnames(gtf, colnames(gtf), gtf_colnames)


  gtf
}
