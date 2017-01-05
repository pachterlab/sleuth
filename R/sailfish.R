#' Constructor for a 'sleuth' object using sailfish data
#'
#' Converts sailfish results to kallisto HDF5
#' format, then return the results of \code{\link[sleuth]{sleuth_prep}} on
#' the converted data.
#'
#' See \code{\link[sleuth]{sleuth_prep}} for parameters other than sf_dirs.
#'
#' @param sf_dirs a character vector of length greater than one where each
#' string points to a sailfish output directory
#' @export
sleuth_prep_sailfish <- function(
  sf_dirs,
  sample_to_covariates,
  full_model,
  filter_fun = basic_filter,
  target_mapping = NULL,
  max_bootstrap = NULL,
  ...) {

  # check sf_dirs
  if( !is(sf_dirs, 'character') ) {
    stop(paste0('"', substitute(sf_dirs),
                '" (sf_dirs) is must be a character vector.'))
  }

  kal_dirs <- package_sf_as_kal(sf_dirs)

  library("sleuth")
  sleuth::sleuth_prep(
    kal_dirs,
    sample_to_covariates,
    full_model,
    filter_fun,
    target_mapping,
    max_bootstrap
  )
}

#' Convert sailfish results for one or more samples to kallisto HDF5
#'
#' @param sf_dirs a character vector of length greater than one where each
#' string points to a sailfish output directory
package_sf_as_kal <- function(sf_dirs) {
  sapply(sf_dirs, sf_to_hdf5)
  sf_dirs
}

#' Convert sailfish results for one sample to kallisto HDF5
#'
#' @param sf_dir path to a sailfish output directory
sf_to_hdf5 <- function(sf_dir) {
  library(data.table)

  h5file <- file.path(sf_dir, 'abundance.h5')
  if (file.exists(h5file)) {
    print(paste("Skipping conversion: abundance.h5 already in ", sf_dir))
    return()
  }

  # load quantification data
  quant <- fread(file.path(sf_dir, 'quant.sf'))
  setnames(quant, c('target_id', 'length', 'tpm', 'est_counts'))
  setkey(quant, 'target_id')

  # load bootstrap data if it exists
  bootspath <- file.path(sf_dir, 'quant_bootstraps.sf')
  numBoot <- 0
  if (file.exists(bootspath)) {
    boots <- fread(bootspath)
    target_ids <- names(boots)
    boots <- data.table(t(boots))
    setnames(boots, sapply(0:(ncol(boots)-1), function(i) paste('bs', i, sep='')))
    numBoot <- ncol(boots)
    boots[, target_id:=target_ids]
    setkey(boots, 'target_id')
    quant <- merge(quant, boots)
  }

  # load stats
  stats_tbl <- fread(file.path(sf_dir, 'stats.tsv'))
  stats <- stats_tbl$V2
  names(stats) <- stats_tbl$V1
  stats_tbl <- stats_tbl[-1]
  setnames(stats_tbl, c('target_id', 'eff_length'))
  setkey(stats_tbl, 'target_id')
  quant <- merge(quant, stats_tbl)

  numProcessed <- stats[['numObservedFragments']]

  # build the hdf5
  library(rhdf5)
  h5createFile(h5file)

  # counts are at root
  h5write(quant$est_counts, h5file, 'est_counts')

  # aux group has metadata about the run and targets
  h5createGroup(h5file, 'aux')
  h5write(numProcessed, h5file, 'aux/num_processed')
  h5write(numBoot, h5file, 'aux/num_bootstrap')
  h5write(quant$length, h5file, 'aux/lengths')
  h5write(quant$eff_length, h5file, 'aux/eff_lengths')
  h5write(quant$target_id, h5file, 'aux/ids')
  h5write('10', h5file, 'aux/index_version')
  h5write('sailfish', h5file, 'aux/kallisto_version')
  h5write(timestamp(prefix="", suffix=""), h5file, "aux/start_time")

  # bootstrap group has (.. wait for it ..) bootstrap data
  if (numBoot > 0) {
    h5createGroup(h5file, 'bootstrap')
    sapply(0:(numBoot-1), function(i) {
      bootid <- paste('bs', i, sep='')
      h5write(unlist(quant[, bootid, with=FALSE]),
              h5file, paste('bootstrap', bootid, sep='/'))
    })
  }

  print(paste("Successfully converted Sailfish results in", sf_dir, "to kallisto HDF5 format"))
}
