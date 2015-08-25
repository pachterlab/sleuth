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

#' Read a kallisto object from an HDF5 file
#'
#' Read a kallisto object from an HDF5 file.
#'
#' @param fname the file name for the HDF5 file
#' @param read_bootstrap if \code{TRUE} load bootstraps, otherwise do not
#' @param max_bootstrap an integer denoting the number of bootstraps to read.
#' if \code{NULL} read everything available
#' @return a \code{kallisto} object
#' @export
read_kallisto_h5 <- function(fname, read_bootstrap = TRUE, max_bootstrap = NULL) {
  stopifnot(is(fname, "character"))

  fname <- path.expand(fname)

  if (!file.exists(fname)) {
    stop("Can't find file: '", fname, "'")
  }

  target_id <- as.character(rhdf5::h5read(fname, "aux/ids"))
  if ( length(target_id) != length(unique(target_id))) {
    tid_counts <- table(target_id)
    warning('Some target_ids in your kallisto index are exactly the same. We will make these unique but strongly suggest you change the names of the FASTA and recreate the index.',
      ' These are the repeats: ',
      paste(names(tid_counts[which(tid_counts > 1)]), collapse = ', '))
    rm(tid_counts)

    target_id <- make.unique(target_id, sep = '_')
  }

  abund <- adf(target_id = target_id)
  abund$est_counts <- as.numeric(rhdf5::h5read(fname, "est_counts"))
  abund$eff_len <- as.numeric(rhdf5::h5read(fname, "aux/eff_lengths"))
  abund$len <- as.numeric(rhdf5::h5read(fname, "aux/lengths"))

  num_processed <- if ( h5check(fname, '/aux', 'num_processed') ) {
    as.integer(rhdf5::h5read(fname, 'aux/num_processed'))
  } else {
    NA_integer_
  }

  bs_samples <- list()
  if (read_bootstrap) {
    num_bootstrap <- as.integer(rhdf5::h5read(fname, "aux/num_bootstrap"))
    if (num_bootstrap > 0) {
      msg("Found ", num_bootstrap, " bootstrap samples")
      if (!is.null(max_bootstrap) && max_bootstrap < num_bootstrap) {
        msg("Only reading ", max_bootstrap, " bootstrap samples")
        num_bootstrap <- max_bootstrap
      }
      bs_samples <- lapply(0:(num_bootstrap[1]-1), function(i)
        {
          .read_bootstrap_hdf5(fname, i, abund)
        })
    } else {
      msg("No bootstrap samples found")
    }
  }

  abund$tpm <- counts_to_tpm(abund$est_counts, abund$eff_len)

  res <- list(abundance = abund, bootstrap = bs_samples)
  class(res) <- 'kallisto'

  attr(res, 'index_version') <- rhdf5::h5read(fname, 'aux/index_version')
  attr(res, 'kallisto_version') <- rhdf5::h5read(fname, 'aux/kallisto_version')
  attr(res, 'start_time') <- rhdf5::h5read(fname, 'aux/start_time')
  attr(res, 'num_targets') <- nrow(abund)
  attr(res, 'num_mapped') <- sum(abund$est_counts)
  attr(res, 'num_processed') <- num_processed

  invisible(res)
}

h5check <- function(fname, group, name) {
  objs <- rhdf5::h5ls(fname)
  objs <- dplyr::rename(objs, grp = group, nm = name)
  objs <- dplyr::filter(objs, grp == group, nm == name)

  nrow(objs) == 1
}

# read a bootstrap from an HDF5 file and return a \code{data.frame}
.read_bootstrap_hdf5 <- function(fname, i, main_est) {
  bs <- adf( target_id = main_est$target_id )
  bs$est_counts <- as.numeric(rhdf5::h5read(fname, paste0("bootstrap/bs", i)))
  bs$tpm <- counts_to_tpm(bs$est_counts, main_est$eff_len)

  bs
}

# Read a kallisto data set
#
# Read a kallisto data set
#
# @param output_dir the directory of the output data
# @param read_bootstrap if TRUE, then searches for bootstrap data, else doesn't read it.
# @return a S3 \code{kallisto} object with the following members:
read_kallisto <- function(output_dir, read_bootstrap = TRUE)
{
  # TODO: function needs to be reworked for plaintext case
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
print.sleuth <- function(obj) {
  cat("\tsleuth object\n")
  cat("\n")
  cat("bears:", length(obj$kal), "\n")
  cat("design:", deparse(obj$full_formula), "\n")
}


# @export
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

# @export
gtf_attributes_to_gene_trans <- function(gtf_attr) {
  stopifnot(is(gtf_attr, "character"))

  lapply(strsplit(gtf_attr, " "), kv_vec_to_df) %>%
    rbind_all()
}

# @export
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

# @export
read_gtf <- function(fname) {
  gtf <- data.table::fread(fname, sep = "\t", header = FALSE,
    data.table = FALSE)

  gtf_colnames <- c("seqname", "source", "feature", "start", "end", "score",
    "strand", "frame", "attribute")

  data.table::setnames(gtf, colnames(gtf), gtf_colnames)


  gtf
}

# @export
trans_to_genes_from_gtf <- function(fname) {
  trans <- rtracklayer::import(fname)
  trans <- data.frame(GenomicRanges::mcols(trans), stringsAsFactors = FALSE)
  trans <- dplyr::select(trans, transcript_id, gene_id)
  trans <- dplyr::distinct(trans)

  trans
}

# Write a kallisto object to HDF5
#
# Write a kallisto object to HDF5.
#
# @param kal the kallisto object to write out
# @param fname the file name to write out to
# @return the kallisto object \code{kal} invisibly.
# @export
write_kallisto_hdf5 <- function(kal, fname, overwrite = TRUE, write_bootstrap = TRUE, compression = 6L) {
  stopifnot( is(kal, "kallisto") )
  stopifnot( is(fname, "character") )
  stopifnot( is(compression, "integer") )
  stopifnot( length(compression) == 1 )

  # TODO: ensure that all bootstraps are sorted according to abundance
  if (compression < 0 || compression > 7 ) {
    stop("'compression' must be in [0, 7]")
  }

  fname <- path.expand(fname)

  if (file.exists(fname)) {
    if (overwrite) {
      warning(paste0("'", fname, "' already exists. Overwritting."))
      file.remove(fname)
    } else {
      stop(paste0("'", fname, "' already exists."))
    }
  }

  if ( !rhdf5::h5createFile( fname ) ) {
    stop(paste0("Error: Couldn't open '", fname, "' to write out."))
  }

  dims <- c(nrow(kal$abundance), 1)
  cat("dims: ", class(dims), "\n")

  # write out auxilary info
  rhdf5::h5createGroup(fname, "aux")
  # stopifnot( rhdf5::h5writeDataset(fname, "aux/ids", dims = dims,
  #   storage.mode = "character", size = 100, level = compression) )
  stopifnot( rhdf5::h5writeDataset(fname, "aux/ids", dims = dims,
    storage.mode = "character", size = 100, level = compression) )
  # rhdf5::h5write(kal$abundance$target_id, fname, "aux/ids")

  if (write_bootstrap) {
    rhdf5::h5write(length(kal$bootstrap), fname, "aux/num_bootstrap")
  } else {
    rhdf5::h5write(0L, fname, "aux/num_bootstrap")
  }

  rhdf5::h5write(kal$abundance$eff_len, fname, "aux/eff_lengths")
  rhdf5::h5write(kal$abundance$len, fname, "aux/lengths")

  # TODO: put lengths in aux
  # TODO: put effective lengths in aux
  # TODO: put version in aux

  # rhdf5::h5createDataset(fname, "est_counts",
  #   storage.mode = "double", level = compression)
  rhdf5::h5write(kal$abundance$est_counts, fname, "est_counts")

  if (write_bootstrap && length(kal$bootstrap) > 0) {
    rhdf5::h5createGroup(fname, "bootstrap")
    for (i in seq_along(kal$bootstrap)) {
      bs <- kal$bootstrap[[i]]$est_counts
      bs_name <- paste0("bootstrap/bs", i)
      # rhdf5::h5createDataset(fname, bs_name,
      #   storage.mode = "double", level = compression)
      rhdf5::h5write(bs, fname, bs_name)
    }
  }

  invisible(kal)
}
