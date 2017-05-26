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

#' Read kallisto output
#'
#' This is a general driver function to read kallisto output. It can read either
#' H5 or tsv.
#' @param either the kallisto directory name or the file name of a h5 or tsv output from kallisto
#' @param read_bootstrap if \code{TRUE}, bootstraps will be read (h5 only)
#' @param max_bootstrap an integer denoting the number of bootstraps to read.
#' if \code{NULL} read everything available (h5 only)
#' @return a \code{kallisto} object
#' @export
read_kallisto <- function(path, read_bootstrap = TRUE, max_bootstrap = NULL) {
  stopifnot(is(path, "character"))

  kal_path <- get_kallisto_path(path)

  if ( kal_path$ext == "tsv" && read_bootstrap ) {
    warning("You specified to read bootstraps, but we won't do so for plaintext")
  }

  result <- NULL
  if ( kal_path$ext == "h5" ) {
    result <- read_kallisto_h5(kal_path$path, read_bootstrap = read_bootstrap,
      max_bootstrap = max_bootstrap)
  } else {
    result <- read_kallisto_tsv(kal_path$path)
  }

  result
}

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
  stopifnot( is.null(max_bootstrap) ||
    is(max_bootstrap, "numeric") ||
    is(max_bootstrap, "integer") )

  fname <- path.expand(fname)

  if (!file.exists(fname)) {
    stop("Can't find file: '", fname, "'")
  }

  target_id <- as.character(rhdf5::h5read(fname, "aux/ids"))
  if ( length(target_id) != length(unique(target_id))) {
    tid_counts <- table(target_id)
    warning(
      'Some target_ids in your kallisto index are exactly the same.',
      ' We will make these unique but strongly suggest you change the names',
      ' of the FASTA and recreate the index.',
      ' These are the repeats: ',
      paste(names(tid_counts[which(tid_counts > 1)]), collapse = ', '))
    rm(tid_counts)

    target_id <- make.unique(target_id, sep = '_')
  }

  abund <- adf(target_id = target_id)
  abund$est_counts <- as.numeric(rhdf5::h5read(fname, "est_counts"))
  abund$eff_len <- as.numeric(rhdf5::h5read(fname, "aux/eff_lengths"))
  abund$len <- as.numeric(rhdf5::h5read(fname, "aux/lengths"))

  num_processed <- if ( h5check(fname, '/aux', 'num_processed') ) { # nolint
    as.integer(rhdf5::h5read(fname, 'aux/num_processed'))
  } else {
    NA_integer_
  }

  fld <- if ( h5check(fname, '/aux', 'fld') ) { # nolint
    as.integer(rhdf5::h5read(fname, 'aux/fld'))
  } else {
    NA_integer_
  }

  bias_observed <- NA
  bias_normalized <- NA
  if ( h5check(fname, '/aux', 'bias_observed') ) { # nolint
    bias_observed <- rhdf5::h5read(fname, 'aux/bias_observed')
    bias_normalized <- rhdf5::h5read(fname, 'aux/bias_normalized')
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
      bs_samples <- lapply(0:(num_bootstrap[1] - 1),
        function(i) {
          .read_bootstrap_hdf5(fname, i, abund)
        })
    } else {
      msg("No bootstrap samples found")
    }
  }

  abund$tpm <- counts_to_tpm(abund$est_counts, abund$eff_len)

  res <- list(
    abundance = abund,
    bias_normalized = bias_normalized,
    bias_observed = bias_observed,
    bootstrap = bs_samples,
    fld = fld
    )
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

# @return a matrix with each row being a bootstrap sample
read_bootstrap_mat <- function(fname,
  num_bootstraps,
  num_transcripts,
  est_count_sf) {
  bs_mat <- matrix(ncol = num_bootstraps, nrow = num_transcripts)
  for (i in 1:ncol(bs_mat)) {
    bs_mat[, i] <- rhdf5::h5read(fname, paste0("bootstrap/bs", i - 1)) / est_count_sf
  }
  bs_mat <- t(bs_mat)
  target_id <- as.character(rhdf5::h5read(fname, "aux/ids"))
  colnames(bs_mat) <- target_id

  bs_mat
}

#' Read kallisto plaintext output
#'
#' This function reads kallisto plaintext output. Note, it cannot be used with
#' sleuth. It also does not read bootstraps since reading plaintext bootstraps
#' is quite slow.
#' @param fname the filename for the tsv file
#' @return a \code{kallisto} object (currently missing attributes and
#' bootstraps)
#' @export
read_kallisto_tsv <- function(fname) {
  stopifnot(is(fname, "character"))

  fname <- path.expand(fname)

  if (!file.exists(fname)) {
    stop("Can't find file: '", fname, "'")
  }

  abundance <- suppressWarnings(data.table::fread(fname, data.table = FALSE))
  abundance <- dplyr::rename(abundance,
    len = length,
    eff_len = eff_length)
  abundance <- dplyr::arrange(abundance, target_id)

  result <- list(abundance = abundance, bootstrap = NULL)
  class(result) <- 'kallisto'

  invisible(result)
}

# this function takes a path and tries to infer whether it is a h5, tsv, or neither
get_kallisto_path <- function(path) {
  output <- list()

  if ( dir.exists(path) ) {
    if (file.exists(file.path(path, "abundance.h5"))) {
      # standard case where the user has not changed the filename
      output$ext <- "h5"
      output$path <- file.path(path, "abundance.h5")
    } else if ( file.exists(file.path(path, 'abundance.tsv')) ){
      # HDF5 doesn't exist, but we have plaintext
      output$ext <- "tsv"
      output$path <- file.path(path, "abundance.tsv")
    } else {
      stop(path, 'exists, but does not contain kallisto output (abundance.h5)')
    }
  } else if ( file.exists(path) ){
    # make an assumption that the user has kept the correct extension
    base <- basename(path)
    s <- strsplit(base, '\\.')
    ext <- s[[1]][length(s[[1]])]

    if (ext == 'h5') {
      output$ext <- 'h5'
    } else if (ext == 'tsv') {
      output$ext <- 'tsv'
    } else {
      stop("'", path, "' exists, is not a recognized extension")
    }
    output$path <- path
  } else {
    stop("'", path, "' does not exist.")
  }

  output
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
    while ( (nchar(gene_id[i]) < 1 || nchar(trans_id[i]) < 1) &&
      j <= length(all_attr[[i]]) ) {
      if (all_attr[[i]][j] == "gene_id") {
        gene_id[i] <- all_attr[[i]][j + 1] %>%
          gsub('"', "", .) %>%
          sub(";", "", .)
        j <- j + 2
      } else if (all_attr[[i]][j] == "transcript_id") {
        trans_id[i] <- all_attr[[i]][j + 1] %>%
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

#' save a sleuth object
#'
#' save a sleuth object
#'
#' @param obj a \code{sleuth} object
#' @param the location to save the object to
#' @seealso \code{\link{sleuth_load}}, \code{\link{sleuth_deploy}}
#' @export
sleuth_save <- function(obj, file) {
  if (!is(obj, 'sleuth')) {
    stop('please provide a sleuth object')
  }
  saveRDS(obj, file=file)
}

#' load a sleuth object
#'
#' load a sleuth object previously saved with \code{sleuth_save}
#'
#' @param file the file to load
#' @return a \code{sleuth} object
#' @seealso \code{\link{sleuth_save}}, \code{\link{sleuth_deploy}}
#' @export
sleuth_load <- function(file) {
  readRDS(file)
}
