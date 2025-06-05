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

#' Basic row filter
#'
#' A basic filter to be used.
#'
#' @param row this is a vector of numerics that will be passedin
#' @param min_reads the minimum mean number of reads
#' @param min_prop the minimum proportion of samples to pass this filter
#' @return a logical of length 1
#' @export
basic_filter <- function(row, min_reads = 5, min_prop = 0.47) {
  mean(row >= min_reads) >= min_prop
}

#' Natural log and offset transformation
#'
#' The default transformation function for converting the normalized counts.
#'
#' @param x numeric that must be >=0. represents an individual observed count (one transcript in one sample).
#' @param offset numeric offset to prevent taking the log of 0.
#' @return log(x + offset)
#' @export
log_transform <- function(x, offset=0.5) {
  log(x + offset)
}

# currently defunct
filter_df_all_groups <- function(df, fun, group_df, ...) {
  grps <- setdiff(colnames(group_df), 'sample')
  res <- sapply(unique(grps),
    function(g) {
      apply(filter_df_by_groups(df, fun, group_df[, c('sample', g)], ...), 1, any)
    })

  vals <- apply(res, 1, any)
  names(vals) <- rownames(df)

  vals
}

# currently defunct
filter_df_by_groups <- function(df, fun, group_df, ...) {
  stopifnot(ncol(group_df) == 2)
  grps <- as.character(group_df[[setdiff(colnames(group_df), 'sample')]])
  sapply(unique(grps),
    function(g) {
      valid_samps <- grps %in% g
      valid_samps <- group_df$sample[valid_samps]
      apply(df[, valid_samps], 1, fun, ...)
    })
}

#' Constructor for a 'sleuth' object
#'
#' A sleuth is a group of kallistos. Borrowing this terminology, a 'sleuth' object stores
#' a group of kallisto results, and can then operate on them while
#' accounting for covariates, sequencing depth, technical and biological
#' variance.
#'
#' @param sample_to_covariates a \code{data.frame} which contains a mapping
#' from \code{sample} (a required column) to some set of experimental conditions or
#' covariates. The column \code{path} is also required, which is a character
#' vector where each element points to the corresponding kallisto output directory. The column
#' \code{sample} should be in the same order as the corresponding entry in
#' \code{path}.
#' @param full_model an R \code{formula} which explains the full model (design)
#' of the experiment OR a design matrix. It must be consistent with the data.frame supplied in
#' \code{sample_to_covariates}. You can fit multiple covariates by joining them with '+' (see example)
#' @param target_mapping a \code{data.frame} that has at least one column
#' 'target_id' and others that denote the mapping for each target. if it is not
#' \code{NULL}, \code{target_mapping} is joined with many outputs where it
#' might be useful. For example, you might have columns 'target_id',
#' 'ensembl_gene' and 'entrez_gene' to denote different transcript to gene
#' mappings. Note that sleuth_prep will treat all columns as having the 'character' data type.
#' @param aggregation_column a string of the column name in \code{\link{target_mapping}} to aggregate targets
#' (typically to summarize the data on the gene level). The aggregation is done using a p-value aggregation
#' method when generating the results table. See \code{\link{sleuth_results}} for more information.
#' @param num_cores an integer of the number of computer cores mclapply should use
#' to speed up sleuth preparation
#' @param ... any of several other arguments that can be used as advanced options for
#' sleuth preparation. See details.
#'
#' @details This method takes a list of samples with kallisto results and returns a sleuth
#'   object with the defined normalization of the data across samples (default is the DESeq method;
#'   See \code{\link{basic_filter}}), and then the defined transformation of the data (default is log(x + 0.5)).
#'   This also collects all of the bootstraps for the modeling done using \code{\link{sleuth_fit}}. This
#'   function also takes several advanced options that can be used to customize your analysis.
#'   Here are the advanced options for \code{sleuth_prep}:
#'
#'   Extra arguments related to Bootstrap Summarizing:
#'   \itemize{
#'     \item \code{extra_bootstrap_summary}: if \code{TRUE}, compute extra summary
#'     statistics for estimated counts. This is not necessary for typical analyses; it is only needed
#'     for certain plots (e.g. \code{\link{plot_bootstrap}}). Default is \code{FALSE}.
#'     \item \code{read_bootstrap_tpm}: read and compute summary statistics on bootstraps on the TPM.
#'     This is not necessary for typical analyses; it is only needed for some plots (e.g. \code{\link{plot_bootstrap}})
#'     and if TPM values are used for \code{\link{sleuth_fit}}. Default is \code{FALSE}.
#'     \item \code{max_bootstrap}: the maximum number of bootstrap values to read for each
#'     transcript. Setting this lower than the total bootstraps available will save some time, but
#'     will likely decrease the accuracy of the estimation of the inferential noise.
#'   }
#'
#'   Advanced Options for Filtering:
#'   \itemize{
#'     \item \code{filter_fun}: the function to use when filtering. This function will be applied to the raw counts
#'     on a row-wise basis, meaning that each feature will be considered individually. The default is to filter out
#'     any features that do not have at least 5 estimated counts in at least 47% of the samples (see \code{\link{basic_filter}}
#'     for more information). If the preferred filtering method requires a matrix-wide transformation or otherwise
#'     needs to consider multiple features simultaneously instead of independently, please consider using
#'     \code{filter_target_id} below.
#'     \item \code{filter_target_id}: character vector of target_ids to filter using methods that
#'     can't be implemented using \code{filter_fun}. If non-NULL, this will override \code{filter_fun}.
#'   }
#'
#'   Advanced Options for the Normalization Step:
#'   (NOTE: Be sure you know what you're doing before you use these options)
#'   \itemize{
#'     \item \code{normalize}: boolean for whether normalization and other steps should be performed.
#'     If this is set to false, bootstraps will not be read and transformation of the data will not be done.
#'     This should only be set to \code{FALSE} if one desires to do a quick check of the raw data.
#'     The default is \code{TRUE}.
#'     \item \code{norm_fun_counts}: a function to perform between sample normalization on the estimated counts.
#'     The default is the DESeq method. See \code{\link{norm_factors}} for details.
#'     \item \code{norm_fun_tpm}: a function to perform between sample normalization on the TPM.
#'     The default is the DESeq method. See \code{\link{norm_factors}} for details.
#'   }
#'
#'   Advanced Options for the Transformation Step:
#'   (NOTE: Be sure you know what you're doing before you use these options)
#'   \itemize{
#'     \item \code{transform_fun_counts}: the transformation that should be applied
#'     to the normalized counts. Default is \code{'log(x+0.5)'} (i.e. natural log with 0.5 offset).
#'     \item \code{transform_fun_tpm}: the transformation that should be applied
#'     to the TPM values. Default is \code{'x'} (i.e. the identity function / no transformation)
#'   }
#'
#'   Advanced Options for Gene Aggregation:
#'   \itemize{
#'     \item \code{gene_mode}: Set this to \code{TRUE} to get the old counts-aggregation method
#'     for doing gene-level analysis. This requires \code{aggregation_column} to be set. If
#'     \code{TRUE}, this will override the p-value aggregation mode, but will allow for gene-centric
#'     modeling, plotting, and results.
#'   }
#'
#' @return a \code{sleuth} object containing all kallisto samples, metadata,
#' and summary statistics
#' @examples # Assume we have run kallisto on a set of samples, and have two treatments,
#' genotype and drug.
#' colnames(s2c)
#' # [1] "sample"  "genotype"  "drug"  "path"
#' so <- sleuth_prep(s2c, ~genotype + drug)
#' @seealso \code{\link{sleuth_fit}} to fit a model, \code{\link{sleuth_wt}} or
#' \code{\link{sleuth_lrt}} to perform hypothesis testing
#' @export
sleuth_prep <- function(
  sample_to_covariates,
  full_model = NULL,
  target_mapping = NULL,
  aggregation_column = NULL,
  num_cores = max(1L, parallel::detectCores() - 1L),
  ...) {

  ##############################
  # check extra options
  extra_opts <- list(...)
  if ("extra_bootstrap_summary" %in% names(extra_opts)) {
    extra_bootstrap_summary <- extra_opts$extra_bootstrap_summary
  } else {
    extra_bootstrap_summary <- FALSE
  }
  if ("read_bootstrap_tpm" %in% names(extra_opts)) {
    read_bootstrap_tpm <- extra_opts$read_bootstrap_tpm
  } else {
    read_bootstrap_tpm <- FALSE
  }
  if ("max_bootstrap" %in% names(extra_opts)) {
    max_bootstrap <- extra_opts$max_bootstrap
  } else {
    max_bootstrap <- NULL
  }
  if ("filter_fun" %in% names(extra_opts)) {
    filter_fun <- extra_opts$filter_fun
  } else {
    filter_fun <- basic_filter
  }
  if ("filter_target_id" %in% names(extra_opts)) {
    filter_target_id <- extra_opts$filter_target_id
  } else {
    filter_target_id <- NULL
  }
  if ("normalize" %in% names(extra_opts)) {
    normalize <- extra_opts$normalize
  } else {
    normalize <- TRUE
  }
  if ("norm_fun_counts" %in% names(extra_opts)) {
    norm_fun_counts <- extra_opts$norm_fun_counts
  } else {
    norm_fun_counts <- norm_factors
  }
  if ("norm_fun_tpm" %in% names(extra_opts)) {
    norm_fun_tpm <- extra_opts$norm_fun_tpm
  } else {
    norm_fun_tpm <- norm_factors
  }
  if ("transform_fun_counts" %in% names(extra_opts)) {
    transform_fun_counts <- extra_opts$transform_fun_counts
  } else {
    transform_fun_counts <- log_transform
  }
  if ("transform_fun_tpm" %in% names(extra_opts)) {
    transform_fun_tpm <- extra_opts$transform_fun_tpm
  } else {
    transform_fun_tpm <- identity
  }
  if ("gene_mode" %in% names(extra_opts)) {
    gene_mode <- extra_opts$gene_mode
  } else {
    gene_mode <- FALSE
  }

  # check inputs

  # data types

  if (!is(sample_to_covariates, "data.frame")) {
    stop(paste0("'", substitute(sample_to_covariates), "' (sample_to_covariates) must be a data.frame"))
  }

  if (!is(full_model, "formula") && !is(full_model, "matrix") && !is.null(full_model)) {
    stop(paste0("'", substitute(full_model), "' (full_model) must be a formula or a matrix"))
  }

  if (!("sample" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates),
        "' (sample_to_covariates) must contain a column named 'sample'"))
  }

  if (!("path" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates)),
      "' (sample_to_covariates) must contain a column named 'path'")
  }

  if (!is.null(target_mapping) && !is(target_mapping, 'data.frame')) {
    stop(paste0("'", substitute(target_mapping),
        "' (target_mapping) must be a data.frame or NULL"))
  } else if (is(target_mapping, 'data.frame')){
    if (!("target_id" %in% colnames(target_mapping))) {
      stop(paste0("'", substitute(target_mapping),
          "' (target_mapping) must contain a column named 'target_id'"))
    }
  }

  if (!is.null(max_bootstrap) && max_bootstrap <= 0 ) {
    stop("max_bootstrap must be > 0")
  }

  if (any(is.na(sample_to_covariates))) {
    warning("Your 'sample_to_covariates' data.frame contains NA values. This will likely cause issues later.")
  }

  if (is(full_model, "matrix") &&
      nrow(full_model) != nrow(sample_to_covariates)) {
    stop("The design matrix number of rows are not equal to the number of rows in the sample_to_covariates argument.")
  }

  if (!is(filter_fun, 'function')) {
    stop("filter_fun must be a function")
  }

  if (!is.null(filter_target_id) & !is.character(filter_target_id)) {
    stop("if filter_target_id is used, it must be a character vector")
  }

  if (!is(norm_fun_counts, 'function')) {
    stop("norm_fun_counts must be a function")
  }

  if (!is(norm_fun_tpm, 'function')) {
    stop("norm_fun_tpm must be a function")
  }

  if (!is(transform_fun_counts, 'function')) {
    stop("transform_fun_counts must be a function")
  }

  if (!is(transform_fun_tpm, 'function')) {
    stop("transform_fun_tpm must be a function")
  }

  if (is.null(aggregation_column) && gene_mode) {
    stop("You set 'gene_mode' to TRUE, but did not provide an 'aggregation_column' ",
         "to aggregate by. Please provide an 'aggregation_column'.")
  } else if (gene_mode) {
    message("'gene_mode' is TRUE. Sleuth will do counts aggregation at the gene level ",
            "for downstream normalization, transformation, and modeling steps, as well as ",
            "for plotting and results.")
  }

  if (!is.null(aggregation_column) && is.null(target_mapping)) {
    stop(paste("You provided an 'aggregation_column' to aggregate by,",
               "but not a 'target_mapping'. Please provide a 'target_mapping'."))
  }

  pval_aggregate <- !is.null(aggregation_column) && !gene_mode
  num_cores <- check_num_cores(num_cores)

  # TODO: ensure transcripts are in same order -- if not, report warning that
  # kallisto index might be incorrect

  # done
  ##############################

  msg('reading in kallisto results')
  sample_to_covariates <- as.data.frame(sample_to_covariates)
  sample_to_covariates$sample <- as.character(sample_to_covariates$sample)

  if(nrow(sample_to_covariates) == 1 && !is.null(full_model)) {
    warning("There is only one sample present, but you also provided a model. ",
            "The model will be set to NULL to prevent downstream errors.\n",
            "The sample can be viewed using sleuth_live after preparation, ",
            "but you need more than one sample to run the other aspects of Sleuth.")
    full_model <- NULL
  }

  kal_dirs <- sample_to_covariates$path
  sample_to_covariates$path <- NULL

  msg('dropping unused factor levels')
  sample_to_covariates <- droplevels(sample_to_covariates)

  nsamp <- 0
  # append sample column to data
  kal_list <- lapply(seq_along(kal_dirs),
    function(i) {
      nsamp <- dot(nsamp)
      path <- kal_dirs[i]
      suppressMessages({
        kal <- read_kallisto(path, read_bootstrap = FALSE,
          max_bootstrap = max_bootstrap)
        })
      kal$abundance <- dplyr::mutate(kal$abundance,
        sample = sample_to_covariates$sample[i])

      kal
    })
  msg('')

  check_result <- check_kal_pack(kal_list)
  kal_versions <- check_result$versions

  obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))

  counts_test <- data.table::as.data.table(obs_raw)
  counts_test <- counts_test[, .(total = sum(est_counts)), by = "sample"]
  if (any(counts_test$total == 0)) {
    zero_names <- counts_test$sample[which(counts_test$total == 0)]
    formatted_names <- paste(zero_names, collapse = ", ")
    warning("At least one sample have no reads aligned. ",
            "Here are the samples with zero counts:\n",
            formatted_names)
  }

  design_matrix <- NULL
  if (is(full_model, 'formula')) {
    design_matrix <- model.matrix(full_model, sample_to_covariates)
  } else if (is(full_model, 'matrix')) {
    if (is.null(colnames(full_model))) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    design_matrix <- full_model
  }

  if (!is.null(full_model)) {
    rownames(design_matrix) <- sample_to_covariates$sample
    # check if the resulting design_matrix is singular (i.e. non-invertible)
    # followed the suggested method found here: https://stackoverflow.com/a/24962470
    M <- t(design_matrix) %*% design_matrix
    det_mod <- determinant(M)$modulus
    if(!is.finite(det_mod)) {
      stop("The full model you provided seems to result in a singular design matrix. ",
           "This frequently happens when one of the covariates is a linear ",
           "combination of one or more other covariates (e.g. one covariate ",
           "yields identical groupings as another covariate). Check your ",
           "sample_to_covariates table and your full model.")
    }
    rm(M, det_mod)
  }

  obs_raw <- dplyr::arrange(obs_raw, target_id, sample)

  ###
  # try to deal with weird ensemble names
  ###
  if (!is.null(target_mapping)) {
    tmp_names <- data.frame(target_id = kal_list[[1]]$abundance$target_id,
      stringsAsFactors = FALSE)
    target_mapping <- check_target_mapping(tmp_names, target_mapping,
                                           !is.null(aggregation_column))
    rm(tmp_names)
  }

  ret <- list(
      kal = kal_list,
      kal_versions = kal_versions,
      obs_raw = obs_raw,
      sample_to_covariates = sample_to_covariates,
      bootstrap_summary = NA,
      full_formula = full_model,
      design_matrix = design_matrix,
      target_mapping = target_mapping,
      gene_mode = gene_mode,
      gene_column = aggregation_column,
      norm_fun_counts = norm_fun_counts,
      norm_fun_tpm = norm_fun_tpm,
      transform_fun_counts = transform_fun_counts,
      transform_fun_tpm = transform_fun_tpm,
      pval_aggregate = pval_aggregate
    )

  if (normalize ) {

    msg("normalizing est_counts")
    est_counts_spread <- spread_abundance_by(obs_raw, "est_counts",
      sample_to_covariates$sample)
    if(!is.null(filter_target_id)) {
       msg("A list of target IDs for filtering was found. Using this for filtering")
       target_ids <- rownames(est_counts_spread)
       filter_bool <- target_ids %in% filter_target_id
       names(filter_bool) <- target_ids
    } else {
      filter_bool <- apply(est_counts_spread, 1, filter_fun)
    }
    filter_true <- filter_bool[filter_bool]

    if (sum(filter_bool) == 0) {
      stop("Zero targets passed the filter you used. Please double check the filter used.")
    }

    msg(paste0(sum(filter_bool), ' targets passed the filter'))
    est_counts_sf <- norm_fun_counts(est_counts_spread[filter_bool, , drop = FALSE])

    filter_df <- adf(target_id = names(filter_true))

    est_counts_norm <- as_df(t(t(est_counts_spread) / est_counts_sf))

    est_counts_norm$target_id <- rownames(est_counts_norm)
    est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)

    obs_norm <- est_counts_norm
    obs_norm$target_id <- as.character(obs_norm$target_id)
    obs_norm$sample <- as.character(obs_norm$sample)
    rm(est_counts_norm)

    # deal w/ TPM
    msg("normalizing tpm")
    tpm_spread <- spread_abundance_by(obs_raw, "tpm",
      sample_to_covariates$sample)
    tpm_sf <- norm_fun_tpm(tpm_spread[filter_bool, , drop = FALSE])
    tpm_norm <- as_df(t(t(tpm_spread) / tpm_sf))
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
    tpm_norm$sample <- as.character(tpm_norm$sample)

    msg('merging in metadata')
    # put everyone in the same order to avoid a slow join
    obs_norm <- dplyr::arrange(obs_norm, target_id, sample)
    tpm_norm <- dplyr::arrange(tpm_norm, target_id, sample)

    stopifnot(all.equal(obs_raw$target_id, obs_norm$target_id) &&
      all.equal(obs_raw$sample, obs_norm$sample))

    suppressWarnings({
      if (!all.equal(dplyr::select(obs_norm, target_id, sample),
          dplyr::select(tpm_norm, target_id, sample))) {
        stop('Invalid column rows. In principle, can simply join. Please report error.')
      }

      # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
      #   by = c('target_id', 'sample'))
      obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, tpm))
    })

    # add in eff_len and len
    obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(obs_raw, eff_len, len))


    obs_norm <- as_df(obs_norm)
    ret$obs_norm <- obs_norm
    ret$est_counts_sf <- est_counts_sf
    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
    ret$obs_norm_filt <- dplyr::semi_join(obs_norm, filter_df, by = 'target_id')
    ret$tpm_sf <- tpm_sf

    #### This code through the for loop is a candidate for moving to another function
    path <- kal_dirs[1]
    kal_path <- get_kallisto_path(path)
    target_id <- as.character(rhdf5::h5read(kal_path$path, "aux/ids"))
    num_transcripts <- length(target_id)
    ret$bs_quants <- list()

    which_target_id <- ret$filter_df$target_id

    if (ret$gene_mode) {
      msg(paste0("aggregating by column: ", aggregation_column))
      # Get list of IDs to aggregate on (usually genes)
      # Also get the filtered list and update the "filter_df" and "filter_bool"
      # variables for the sleuth object
      target_mapping <- data.table::data.table(target_mapping)
      target_mapping[target_mapping[[aggregation_column]] == "",
                     aggregation_column] <- NA
      agg_id <- unique(target_mapping[, aggregation_column, with = FALSE])
      agg_id <- agg_id[[1]]
      agg_id <- agg_id[!is.na(agg_id)]
      mappings <- dplyr::select_(target_mapping, "target_id", aggregation_column)
      mappings <- data.table::as.data.table(mappings)
      which_tms <- which(mappings$target_id %in% which_target_id)
      which_agg_id <- unique(mappings[which_tms, aggregation_column, with = FALSE])
      which_agg_id <- which_agg_id[[1]]
      which_agg_id <- which_agg_id[!is.na(which_agg_id)]
      filter_df <- adf(target_id = which_agg_id)
      filter_bool <- agg_id %in% which_agg_id

      msg(paste0(length(which_agg_id), " genes passed the filter"))

      # Taken from gene_summary; scale normalized observed counts to "reads/base"
      norm_by_length <- TRUE
      tmp <- data.table::as.data.table(ret$obs_raw)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = TRUE)
      scale_factor <- tmp[, scale_factor := median(eff_len),
                          by=list(sample,eval(parse(text=aggregation_column)))]
      obs_norm_gene <- reads_per_base_transform(ret$obs_norm,
          scale_factor, aggregation_column, mappings, norm_by_length)
      # New code: get gene-level TPM (simple sum of normalized transcript TPM)
      tmp <- data.table::as.data.table(tpm_norm)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = T)
      if (any(is.na(tmp[[aggregation_column]]))) {
        rows_to_remove <- is.na(tmp[[aggregation_column]])
        num_missing <- length(unique(tmp[rows_to_remove, target_id]))
        warning(num_missing, " target_ids are missing annotations for the aggregation_column: ",
                aggregation_column, ".\nThese target_ids will be dropped from the gene-level analysis.",
                "\nIf you did not expect this, check your 'target_mapping' table for missing values.")
        tmp <- tmp[!rows_to_remove]
      }
      tpm_norm_gene <- tmp[, j = list(tpm = sum(tpm)),
                           by = list(sample, eval(parse(text = aggregation_column)))]
      data.table::setnames(tpm_norm_gene, 'parse', 'target_id')
      tpm_norm_gene <- as_df(tpm_norm_gene)

      # Same steps as above to add TPM column to "obs_norm" table
      obs_norm_gene <- dplyr::arrange(obs_norm_gene, target_id, sample)
      tpm_norm_gene <- dplyr::arrange(tpm_norm_gene, target_id, sample)

      stopifnot(all.equal(dplyr::select(obs_norm_gene, target_id, sample),
            dplyr::select(tpm_norm_gene, target_id, sample)))
      suppressWarnings({
        if ( !all.equal(dplyr::select(obs_norm_gene, target_id, sample),
            dplyr::select(tpm_norm_gene, target_id, sample), check.attributes = FALSE) ) {
              stop('Invalid column rows. In principle, can simply join. Please report error.')
            }

        # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
        #   by = c('target_id', 'sample'))
        obs_norm_gene <- dplyr::bind_cols(obs_norm_gene, dplyr::select(tpm_norm_gene, tpm))
      })

      # These are the updated gene-level variables
      ret$filter_df <- adf(target_id = which_agg_id)
      ret$filter_bool <- agg_id %in% which_agg_id
      ret$obs_norm <- obs_norm_gene
      ret$obs_norm_filt <- dplyr::semi_join(obs_norm_gene, filter_df, by = 'target_id')

      rm(obs_norm, tpm_norm, obs_norm_gene, tpm_norm_gene)

      # This is the gene-level version of the matrix
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_agg_id),
                                     ncol = length(ret$kal))
      which_ids <- which_agg_id
    } else {
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_target_id),
                                     ncol = length(ret$kal))
      which_ids <- which_target_id
    }

    msg('summarizing bootstraps')
    apply_function <- if (num_cores == 1) {
      lapply
    } else {
      function(x, y) parallel::mclapply(x, y, mc.cores = num_cores)
    }
    bs_results <- apply_function(seq_along(kal_dirs), function(i) {
      samp_name <- sample_to_covariates$sample[i]
      kal_path <- get_kallisto_path(kal_dirs[i])
      process_bootstrap(i, samp_name, kal_path,
                        num_transcripts, est_counts_sf[[i]],
                        read_bootstrap_tpm, ret$gene_mode,
                        extra_bootstrap_summary,
                        target_id, mappings, which_ids, ret$gene_column,
                        ret$transform_fun_counts, ret$transform_fun_tpm,
                        max_bootstrap)
    })

    # if mclapply results in an error (a warning is shown), then print error and stop
    error_status <- sapply(bs_results, function(x) is(x, "try-error"))
    if (any(error_status)) {
      error_msgs <- sapply(which(error_status), function(i) {
        bad_run <- bs_results[[i]]
        samp_name <- sample_to_covariates$sample[i]
        trace <- .traceback(bad_run)
        paste0("Sample '", samp_name, "' had this error message: ", trace[1])
      })
      formatted_error <- paste(error_msgs, collapse = "")
      message(formatted_error)
      stop("At least one core from mclapply had an error. See the above error message(s) for more details.")
    }

    # mclapply is expected to retun the bootstraps in order; this is a sanity check of that
    indices <- sapply(bs_results, function(result) result$index)
    stopifnot(identical(indices, order(indices)))

    sigma_q_sq_tpm <- NULL
    if(read_bootstrap_tpm) {
      all_sample_tpm <- sapply(bs_results, function(result) result$bootstrap_tpm_result)
      rownames(all_sample_tpm) <- which_ids
      sigma_q_sq_tpm <- rowMeans(all_sample_tpm)
      sigma_q_sq_tpm <- sigma_q_sq_tpm[order(names(sigma_q_sq_tpm))]
      ret$bs_quants <- lapply(bs_results, function(result) result$bs_quants)
      names(ret$bs_quants) <- sample_to_covariates$sample
    } else if(extra_bootstrap_summary) {
      ret$bs_quants <- lapply(bs_results, function(result) result$bs_quants)
      names(ret$bs_quants) <- sample_to_covariates$sample
    }

    all_sample_bootstrap <- sapply(bs_results, function(result) result$bootstrap_result)
    rownames(all_sample_bootstrap) <- which_ids

    # end summarize bootstraps
    msg('')

    sigma_q_sq <- rowMeans(all_sample_bootstrap)
    names(sigma_q_sq) <- which_ids

    # This is the rest of the gene_summary code
    if (ret$gene_mode) {
      obs_counts <- obs_to_matrix(ret, "scaled_reads_per_base")[which_agg_id, , drop = FALSE]
      obs_tpm <- obs_to_matrix(ret, "tpm")[which_agg_id, , drop = FALSE]
    } else {
      obs_counts <- obs_to_matrix(ret, "est_counts")[which_target_id, , drop = FALSE]
      obs_tpm <- obs_to_matrix(ret, "tpm")[which_target_id, , drop = FALSE]
    }

    sigma_q_sq <- sigma_q_sq[order(names(sigma_q_sq))]
    obs_counts <- ret$transform_fun_counts(obs_counts)
    obs_counts <- obs_counts[order(rownames(obs_counts)),]
    obs_tpm <- ret$transform_fun_tpm(obs_tpm)
    obs_tpm <- obs_tpm[order(rownames(obs_tpm)),]

    ret$bs_summary <- list(obs_counts = obs_counts, sigma_q_sq = sigma_q_sq)
    ret$bs_summary <- list(obs_counts = obs_counts, obs_tpm = obs_tpm,
                           sigma_q_sq = sigma_q_sq,
                           sigma_q_sq_tpm = sigma_q_sq_tpm)
  } else {
    # The filter_bool and filter_df can be done with the raw counts
    # Everything else was skipped because they depend on the normalization step
    # Those items are set to empty vectors, data frames, and lists.
    # Setting the normalization and transformation functions to 'NA' to indicate
    # that these steps were skipped
    ret$norm_fun_counts <- ret$norm_fun_tpm <- NA
    ret$transform_fun_counts <- ret$transform_fun_tpm <- NA
    ret$obs_norm <- data.frame()
    ret$est_counts_sf <- vector()

    est_counts_spread <- spread_abundance_by(obs_raw, "est_counts",
      sample_to_covariates$sample)
    if(!is.null(filter_target_id)) {
       msg("A list of target IDs for filtering was found. Using this for filtering")
       target_ids <- rownames(est_counts_spread)
       filter_bool <- target_ids %in% filter_target_id
       names(filter_bool) <- target_ids
    } else {
      filter_bool <- apply(est_counts_spread, 1, filter_fun)
    }
    filter_true <- filter_bool[filter_bool]

    if (sum(filter_bool) == 0) {
      stop("Zero targets passed the filter you used. Please double check the filter used.")
    }

    msg(paste0(sum(filter_bool), ' targets passed the filter'))
    filter_df <- adf(target_id = names(filter_true))

    if (ret$gene_mode) {
      msg(paste0("aggregating by column: ", aggregation_column))
      which_target_id <- ret$filter_df$target_id
      # Get list of IDs to aggregate on (usually genes)
      # Also get the filtered list and update the "filter_df" and "filter_bool"
      # variables for the sleuth object
      target_mapping <- data.table::data.table(target_mapping)
      target_mapping[target_mapping[[aggregation_column]] == "",
                     aggregation_column] <- NA
      agg_id <- unique(target_mapping[, aggregation_column, with = FALSE])
      agg_id <- agg_id[[1]]
      agg_id <- agg_id[!is.na(agg_id)]
      mappings <- dplyr::select_(target_mapping, "target_id", aggregation_column)
      mappings <- data.table::as.data.table(mappings)
      which_tms <- which(mappings$target_id %in% which_target_id)
      which_agg_id <- unique(mappings[which_tms, aggregation_column, with = FALSE])
      which_agg_id <- which_agg_id[[1]]
      which_agg_id <- which_agg_id[!is.na(which_agg_id)]
      filter_df <- adf(target_id = which_agg_id)
      filter_bool <- agg_id %in% which_agg_id

      msg(paste0(length(which_agg_id), " genes passed the filter"))
    }

    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
    ret$obs_norm_filt <- data.frame()
    ret$tpm_sf <- vector()
    ret$bs_quants <- list()
    ret$bs_summary <- list()
  }

  class(ret) <- 'sleuth'

  ret
}

check_norm_status <- function(obj) {
  if (!is(obj$norm_fun_counts, 'function') || !is(obj$norm_fun_tpm, 'function')) {
    stop("This sleuth object was prepared without normalization. If you wish to do this step,",
         " repeat 'sleuth_prep' with 'normalize' set to TRUE (the default).")
  } else {
    return(TRUE)
  }
}

# check versions of kallistos and num bootstraps, etc
check_kal_pack <- function(kal_list) {

  versions <- sapply(kal_list, function(x) attr(x, 'kallisto_version'))
  u_versions <- unique(versions)
  if (length(u_versions) > 1) {
    warning('More than one version of kallisto was used: ',
      paste(u_versions, collapse = "\n"))
  }

  ntargs <- sapply(kal_list, function(x) attr(x, 'num_targets'))
  u_ntargs <- unique(ntargs)
  if (length(u_ntargs) > 1 ) {
    stop('Inconsistent number of transcripts. Please make sure you used the same index everywhere.')
  }

  list(versions = u_versions)
}

# this function is mostly to deal with annoying ENSEMBL transcript names that
# have a trailing .N to keep track of version number
#
# this also checks to see if there are duplicate entries for any target IDs
# and issues a warning if sleuth prep is in transcript mode, but stops if
# sleuth prep is in gene mode, since duplicate entries creates problems when
# doing the aggregation
#
# finally, this method forces all of the columns in the target_mapping to be
# character columns, to prevent any issues with factors or other data types
# interfering with downstream uses of the target_mapping table
#
# @return the target_mapping if an intersection is found. a target_mapping that
# matches \code{t_id} if no matching is found
check_target_mapping <- function(t_id, target_mapping, gene_mode) {
  t_id <- data.table::as.data.table(t_id)
  t_id$target_id <- as.character(t_id$target_id)
  target_mapping <- apply(target_mapping, 2, as.character)
  target_mapping <- data.table::as.data.table(target_mapping)

  tmp_join <- dplyr::inner_join(t_id, target_mapping, by = 'target_id')

  if (!nrow(tmp_join)) {
    t_id <- dplyr::mutate(t_id, target_id_old = sub('\\.[0-9]+', '', target_id))
    target_mapping_upd <- dplyr::rename(target_mapping, target_id_old = target_id)

    tmp_join <- dplyr::inner_join(t_id, target_mapping_upd,
      c('target_id_old'))

    if (!nrow(tmp_join)) {
      stop("couldn't solve nonzero intersection")
    }

    target_mapping <- dplyr::left_join(target_mapping,
      unique(dplyr::select(tmp_join, target_id_new = target_id,
        target_id = target_id_old),
        by = NULL),
      by = 'target_id')
    # if we couldn't recover, use the old one
    target_mapping <- dplyr::mutate(target_mapping,
      target_id_new = ifelse(is.na(target_id_new), target_id, target_id_new))
    target_mapping <- dplyr::select(target_mapping, -target_id)
    target_mapping <- dplyr::rename(target_mapping, target_id = target_id_new)

    warning(paste0('intersection between target_id from kallisto runs and ',
      'the target_mapping is empty. attempted to fix problem by removing .N ',
      'from target_id, then merging back into target_mapping.',
      ' please check obj$target_mapping to ensure this new mapping is correct.'))
  }

  if(any(duplicated(target_mapping$target_id))) {
    indices <- which(duplicated(target_mapping$target_id))
    duplicated_ids <- target_mapping$target_id[indices]
    formatted_ids <- paste(duplicated_ids, collapse = ", ")
    if(gene_mode) {
      stop("There is at least one duplicated target ID in the target mapping. ",
           "Since sleuth prep is in gene aggregation mode, any duplicated ",
           "entries will cause errors during gene aggregation. Here is a list ",
           "of all the duplicated target IDs:\n", formatted_ids)
    } else {
      warning("There is at least one duplicated target ID in the target ",
              "mapping. Since sleuth prep is not in gene aggregation mode, ",
              "duplicated entries will not cause errors during preparation, ",
              "but may cause unexpected behavior when making any plots or ",
              "tables. Here is a list of all the duplicated target IDs:\n",
              formatted_ids)
    }
  }
  adf(target_mapping)
}

#' Normalization factors
#'
#' More or less the DESeq size factors
#'
#' @param mat a matrix of estimated counts
#' @return a vector of length \code{ncol(mat)} with a sequencing depth factor
#' for each sample
#' @export
norm_factors <- function(mat) {
  nz <- apply(mat, 1, function(row) !any(round(row) == 0))
  mat_nz <- mat[nz, , drop = FALSE]
  p <- ncol(mat)
  geo_means <- exp(apply(mat_nz, 1, function(row) mean(log(row))))
  s <- sweep(mat_nz, 1, geo_means, `/`)

  sf <- apply(s, 2, median)
  scaling <- exp( (-1 / p) * sum(log(sf)))

  sf * scaling
}

# Summarize many bootstrap objects
#
# Summarize all the bootstrap samples from a kallisto run. The summarized
# values are then all put into a data.frame and stored in the \code{sleuth}
# object.
#
# @param obj a \code{sleuth} object
# @param force if \code{FALSE}, then will only compute the summary if it has
# not yet been set.
# @param verbose if \code{TRUE}, print verbosely
# @return a \code{kallisto} object with member \code{bootstrap_summary}
# updated and containing a data frame.
sleuth_summarize_bootstrap <- function(obj, force = FALSE, verbose = FALSE) {
  # TODO: make this into an S3 function 'summarize_bootstrap'
  stopifnot(is(obj, "sleuth"))

  if (!is.na(obj$bootstrap_summary) && !force) {
    if (verbose) {
      cat("summarize_bootstrap.sleuth: Already computed summary -- will not recompute.")
    }

    return(obj)
  }

  tpm_bs <- sleuth_summarize_bootstrap_col(obj, "tpm")
  counts_bs <- sleuth_summarize_bootstrap_col(obj, "est_counts")

  s_bs <- inner_join(
    data.table::data.table(tpm_bs),
    data.table::data.table(data.table(counts_bs)),
    #by = c("target_id", "sample", "condition")
    by = c("target_id", "sample")
    )
  s_bs <- adf(s_bs)
  obj$bootstrap_summary <- s_bs

  obj
}

sleuth_summarize_bootstrap_col <- function(obj, col, transform = identity) {
  res <- lapply(seq_along(obj$kal), function(i) {
      cur_samp <- obj$sample_to_covariates$sample[i]

      dplyr::mutate(summarize_bootstrap(obj$kal[[i]], col, transform),
        sample = cur_samp)
    })
  dplyr::bind_rows(res)
}

# Spread abundance by a column
#
# Take a data.frame from a sleuth object (e.g. \code{obs_raw}) and cast it
# into a matrix where the rows are the target_ids and the columns are the
# sample ids. The values are the variable you are "spreading" on.
# @param abund the abundance \code{data.frame} from a \code{sleuth} object
# @param var a character array of length one. The variable for which to get "spread" on (e.g. "est_counts").
# @export
spread_abundance_by <- function(abund, var, which_order) {
  # var <- lazyeval::lazy(var)
  abund <- data.table::as.data.table(abund)
  var_spread <- data.table::dcast(abund, target_id ~ sample, value.var = var)
  # there is a discrepancy between data table's sorting of character vectors
  # and how tidyr previously (or the order function) sorts character vectors
  # so next step is needed to make sure the order is correct
  var_spread <- var_spread[order(var_spread$target_id), ]
  var_spread <- as.data.frame(var_spread, stringsAsFactors = FALSE)
  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL
  result <- as.matrix(var_spread)

  result[, which_order, drop = FALSE]
}

#' @importFrom dplyr %>%
#' @export
melt_bootstrap_sleuth <- function(obj) {
  # TODO: make this into a S3 function
  lapply(seq_along(obj$kal), function(i) {
      cur_samp <- obj$sample_to_covariates$sample[i]
      cur_cond <- obj$sample_to_covariates$condition[i]

      melt_bootstrap(obj$kal[[i]]) %>%
        rename(bs_sample = sample) %>%
        mutate(sample = cur_samp, condition = cur_cond)
    }) %>% rbind_all()
}

#' Get kallisto abundance table from a sleuth object
#'
#' Get back the kallisto abundance table from a sleuth object.
#'
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, return the filtered data
#' @param normalized if TRUE, return the normalized data. Otherwise, return
#' the raw data
#' @param include_covariates if TRUE, join the covariates
#' @return a \code{data.frame} with at least columns
#' @export
kallisto_table <- function(obj,
  use_filtered = FALSE,
  normalized = TRUE,
  include_covariates = TRUE) {

  res <- NULL
  if (use_filtered && normalized) {
    res <- obj$obs_norm_filt
  } else if (normalized) {
    res <- obj$obs_norm
  } else {
    res <- obj$obs_raw
  }

  if (use_filtered && !normalized) {
    res <- dplyr::semi_join(
      data.table::as.data.table(res),
      data.table::as.data.table(obj$filter_df),
      by = 'target_id')
  }

  if (include_covariates) {
    res <- dplyr::left_join(
      data.table::as.data.table(res),
      data.table::as.data.table(obj$sample_to_covariates),
      by = 'sample')
  }

  as_df(res)
}

# TODO: deprecate?
# observations to a matrix
#
# observations to a matrix
#
# @param obj is a sleuth object
# @param value_name either "est_counts" or "tpm"
# @return a matrix with the appropriate names
obs_to_matrix <- function(obj, value_name) {

  obs_counts <- reshape2::dcast(obj$obs_norm, target_id ~ sample,
    value.var = value_name)

  obs_counts <- as.data.frame(obs_counts)
  rownames(obs_counts) <- obs_counts$target_id
  obs_counts$target_id <- NULL
  obs_counts <- as.matrix(obs_counts)
  obs_counts <- obs_counts[, obj$sample_to_covariates$sample, drop = FALSE]

  obs_counts
}


# TODO: check if works -- currently untested

#' Get a data.frame from all kallisto objects
#'
#' Build a data.frame from all kallisto objects given a column name
#'
#' @param obj a sleuth object
#' @param col_name a column to extract
#' @importFrom dplyr %>%
#' @return a data.frame with columns \code{sample} and column \code{col_name}
get_col <- function(obj, ...) {
  n_samples <- nrow(obj$design_matrix)
  n_trans <- nrow(obj$kal[[1]]$abundance)

  #which_cols <- as.character(...)
  lapply(seq_along(obj$kal),
    function(i) {
      which_sample <- obj$sample_to_covariates$sample[i]
      dplyr::select_(obj$kal[[i]]$abundance, "target_id",
        .dots = lazyeval::lazy_dots(...)) %>%
          mutate(sample = which_sample)
    }) %>%
      data.table::rbindlist() %>%
      as.data.frame()
}

#' @export
summary.sleuth <- function(obj, covariates = TRUE) {
  mapped_reads <- sapply(obj$kal, function(k) attr(k, 'num_mapped'))
  n_bs <- sapply(obj$kal, function(k) attr(k, 'num_bootstrap_found'))
  n_bs_read <- sapply(obj$kal, function(k) attr(k, 'num_bootstrap_used'))
  n_processed <- sapply(obj$kal, function(k) attr(k, 'num_processed'))

  res <- adf(sample = obj$sample_to_covariates[['sample']],
    reads_mapped = mapped_reads,
    reads_proc = n_processed,
    frac_mapped = round(mapped_reads / n_processed, 4),
    bootstraps_present = n_bs,
    bootstraps_used = n_bs_read
    )
  if (covariates) {
    res <- dplyr::left_join(res, obj$sample_to_covariates, by = 'sample')
  }

  res
}

#' Create a gene table from a sleuth object
#'
#' Take in a sleuth object with added genes, and return a data table in which genes
#' list the most significant transcript mapping to themselves
#'
#' @param obj a \code{sleuth} object
#' @param test a character string denoting which beta to use
#' @param test_type either 'wt' for Wald test or 'lrt' for likelihood ratio test
#' @param which_model a character string denoting which model to use
#' @param which_group a character string denoting which gene group to use
#' @return a \code{data.frame} with the following columns
#'
#' @return gene name; if ext_gene name specified, it will be legible gene name.
#'         If ens_gene name, it will be an Ensemble gene (assuming user followed the vignette)
#' @return most_sig_trancript: Most significant transcript for the given gene
#' @return pval: p-value for the test chosen
#' @return qval: False discovery rate normalized p-value (Benjamini-Hochberg, see: \code{\link{p.adjust}})
#' @return num_transcripts: Total number of transcripts for the given gene
#' @return list_of_transcripts: All transcripts associated with this gene
#' @examples sleuth_genes <- sleuth_gene_table(so, 'conditionIP', test_type ='wt',
#'                                   which_group = 'ext_gene')
#' head(sleuth_genes) # show info for first 5 genes
#' sleuth_genes[1:5, 6] # show transcripts for first 5 genes
#' @export
sleuth_gene_table <- function(obj, test, test_type = 'lrt', which_model = 'full', which_group = 'ens_gene') {

  if (is.null(obj$target_mapping)) {
    stop("This sleuth object doesn't have added gene names.")
  }
  popped_gene_table <- sleuth_results(obj, test, test_type, which_model)

  popped_gene_table <- dplyr::arrange_(popped_gene_table, which_group, ~qval)
  popped_gene_table <- dplyr::group_by_(popped_gene_table, which_group)

  popped_gene_table <- dplyr::summarise_(popped_gene_table,
    target_id = ~target_id[1],
    pval = ~min(pval, na.rm  = TRUE),
    qval = ~min(qval, na.rm = TRUE),
    num_transcripts = ~n(),
    # all_target_ids = ~toString(target_id[1:length(target_id)])
    all_target_ids = ~paste0(target_id[1:length(target_id)], collapse = ',')
    )

  filter_empty <- nchar(popped_gene_table[[which_group]]) == 0 | # empty transcript name
    is.na(popped_gene_table[[which_group]]) | # empty gene name
    is.na(popped_gene_table$qval) # missing q-value
  filter_empty <- !filter_empty
  popped_gene_table <- popped_gene_table[filter_empty, ]

  adf(popped_gene_table)
}


#' Get the names of the transcripts associated to a gene
#'
#' Get the names of the transcripts associated to a gene, assuming genes are added to the input \code{sleuth} object.
#'
#' @param obj a \code{sleuth} object
#' @param test a character string denoting which beta to use
#' @param test_type either 'wt' for wald test or 'lrt' for likelihood ratio test
#' @param which_model a character string denoting which model to use
#' @param gene_colname the name of the column in which the desired gene apperas gene appears. Once genes have been added to a sleuth
#' object, you can inspect the genes names present in your sleuth object via \code{obj$target_mapping}, assuming 'obj' is the name of your sleuth object.
#' This parameter refers to the name of the column that the gene you are searching for appears in. Checkout the column names using \code{names(obj$target_mapping)}
#' @param gene_name a string containing the name of the gene you are interested in
#' @return a vector of strings containing the names of the transcripts that map to a gene
#' @export
transcripts_from_gene <- function(obj, test, test_type,
  which_model, gene_colname, gene_name) {

  if (obj$gene_mode) {
    stop("this sleuth object is in gene mode. Please use 'gene_from_gene' instead.")
  }

  table <- sleuth_results(obj, test, test_type, which_model, pval_aggregate = FALSE)
  table <- dplyr::select_(table, ~target_id, gene_colname, ~qval)
  table <- dplyr::arrange_(table, gene_colname, ~qval)
  if (!(gene_name %in% table[, 2])) {
      stop("Couldn't find gene ", gene_name)
  }
  table$target_id[table[, 2] == gene_name]
}

#' Get the gene ID using other gene identifiers
#'
#' Get the \code{target_id} of a gene using other gene identifiers.
#' The identifiers found under the \code{obj$gene_column} are often
#' difficult to remember (e.g. ensembl gene ID, ENSG00000111640).
#' This function allows a user to find that difficult-to-remember
#' identifier using more-easily-remembered identifiers, such as
#' gene symbol (e.g. "GAPDH").
#'
#' @param obj a \code{sleuth} object
#' @param gene_colname the name of the column containing 'gene_name'.
#'   This parameter refers to the name of the column that the gene you are searching for appears in.
#'   Check the column names using \code{colnames(obj$target_mapping)}.
#' @param gene_name a string containing the name of the gene you are interested in.
#' @return a character vector containing the \code{target_id} of the gene, found under
#'   \code{obj$gene_column} within \code{obj$target_mapping}.
#'   If the column name provided is the same as \code{obj$gene_column}, and the
#'   gene_name used is found, that gene_name will be returned.
#' @examples
#'   \dontrun{gene_from_gene(obj, "gene_symbol", "GAPDH")}
#' @export
gene_from_gene <- function(obj, gene_colname, gene_name) {

  if (is.null(obj$gene_column)) {
    stop("this sleuth object is in transcript mode. Please use 'transcripts_from_gene' instead.")
  }

  table <- as.data.frame(obj$target_mapping)
  if (gene_colname == obj$gene_column) {
    if (!(gene_name %in% table[, obj$gene_column])) {
      stop("Couldn't find gene ", gene_name)
    } else {
      return(gene_name)
    }
  }

  table <- unique(dplyr::select_(table, obj$gene_column, gene_colname))
  if (!(gene_name %in% table[, 2])) {
    stop("Couldn't find gene ", gene_name)
  }
  hits <- unique(table[table[,2] == gene_name, 1])
  if (length(hits) > 1) {
    warning("there was more than one gene ID that matched this identifier; taking the first one")
  }
  hits[1]
 }

#' Change sleuth transform counts function
#'
#' Replace the transformation function of a sleuth object for estimated counts
#' NOTE: if you change the transformation function after having done a fit,
#' the fit(s) will have to be redone using the new transformation.
#' @examples transform_fun_counts(x) <- function(x) identity(x)
#' @export
`transform_fun_counts<-` <- function(obj, fxn) {
  stopifnot(is.function(fxn))
  if(!is.null(obj$fits)) {
    obj$transform_fun_counts <- fxn
    warning(paste("Your sleuth object has fits based on the old transform function.",
                  "Please rerun sleuth_prep and sleuth_fit."))
    obj$fits <- lapply(obj$fits, function(x) {
                         x$transform_synced <- FALSE
                         x
                       })
  } else {
    stop("Your sleuth object was prepared using the old transform function. Please rerun",
         " sleuth_prep using the new transform function.")
  }
  obj
}

#' Change sleuth transform TPM function
#'
#' Replace the transformation function of a sleuth object for TPM values
#'
#' NOTE: if you change the transformation function after having done a fit,
#' the fit(s) will have to be redone using the new transformation.
#' @examples transform_fun_tpm(x) <- function(x) identity(x)
#' @export
`transform_fun_tpm<-` <- function(obj, fxn) {
  stopifnot(is.function(fxn))
  if(!is.null(obj$fits)) {
    obj$transform_fun_tpm <- fxn
    warning(paste("Your sleuth object has fits based on the old transform function.",
                  "Please rerun sleuth_prep and sleuth_fit."))
    obj$fits <- lapply(obj$fits, function(x) {
                         x$transform_synced <- FALSE
                         x
                       })
  } else {
    stop("Your sleuth object was prepared using the old transform function. Please rerun",
         " sleuth_prep using the new transform function.")
  }
  obj
}

#' Extend internal '$<-' for sleuth object
#'
#' @details: This extension is mainly to address two cases:
#' \itemize{
#'   \item case 1: when \code{transform_fun_counts} or
#'   \code{transform_fun_tpm}is changed by user. This informs
#'   the user that \code{sleuth_prep} and \code{sleuth_fit} need
#'   to be rerun.
#'   \item case 2: when gene_mode or pval_aggregate are modified.
#'   This warns the user if there is a conflict between these two
#'   mutually exclusive modes, and overrides the one that was left unmodified
#'   by the user. This also warns the user if unexpected behavior is likely to occur.
#' }
#' Otherwise it acts normally.
#' @examples obj$transform_fun_counts <- function(x) log2(x+0.5)
#' @export
`$<-.sleuth` <- function(obj, name, value) {
  if(name == "norm_fun_counts" || name == "norm_fun_tpm") {
    stop("This sleuth object was prepared using the old normalization function. Please rerun",
         " 'sleuth_prep' using the new normalization function.")
  }
  if(name == "transform_fun_counts" || name == "transform_fun_tpm") {
    if(!is.null(obj$fits)) {
      obj[[name]] <- value
      warning(paste("Your sleuth object has fits based on the old transform function.",
                    "Please rerun sleuth_prep and sleuth_fit."))
      obj$fits <- lapply(obj$fits, function(x) {
                           x$transform_synced <- FALSE
                           x
                         })
    } else {
      stop("This sleuth object was prepared using the old transform function. Please rerun",
         " 'sleuth_prep' using the new transform function.")
    }
  }
  if (name == "gene_mode" && value && !obj[[name]] && !is.null(obj$gene_column)) {
    warning("You set 'gene_mode' to TRUE. If you have not already used 'sleuth_prep' with ",
            "'gene_mode' set to TRUE, this will cause unexpected behavior and may break downstream steps.")
    if (obj$pval_aggregate) {
      warning("'pval_aggregate' is TRUE. Setting to FALSE as p-value aggregation mode is mutually exclusive.")
      obj[[name]] <- value
      obj$pval_aggregate <- FALSE
    } else {
      obj[[name]] <- value
    }
  } else if (name == "pval_aggregate" && value && !is.null(obj$gene_column)) {
    stop("You set 'pval_aggregate' to TRUE, but no 'gene_column' is set. Please set a 'gene_column' first.")
  } else if (name == "pval_aggregate" && value && obj$gene_mode) {
    warning("You set 'pval_aggregate' to TRUE, but 'gene_mode' is also TRUE. Setting 'gene_mode' to FALSE.",
            "If 'sleuth_prep' was run using 'gene_mode' set to TRUE, this will lead to unexpected behavior ",
            "and may break downstream steps.")
    obj[[name]] <- value
    obj$gene_mode <- FALSE
  } else {
    obj[[name]] <- value
  }
  obj
}
