propagate_transcript_filter <- function(filter_df, target_mapping,
  grouping_column) {

  filtered_target_mapping <- dplyr::inner_join(as.data.table(filter_df), # nolint
    as.data.table(target_mapping), by = 'target_id') # nolint

  filtered_target_mapping <- dplyr::select_(filtered_target_mapping,
    grouping_column)

  data.table::setnames(filtered_target_mapping, grouping_column, 'target_id')
  filtered_target_mapping <- dplyr::distinct(filtered_target_mapping)

  filtered_target_mapping
}

check_quant_mode <- function(obj, units) {
  stopifnot( is(obj, 'sleuth') )
  if (obj$gene_mode & units == 'est_counts') {
    warning(paste("your sleuth object is in gene mode,",
                  "but you selected 'est_counts'. Selecting 'scaled_reads_per_base'..."))
    units <- 'scaled_reads_per_base'
  } else if (!obj$gene_mode & units == 'scaled_reads_per_base') {
    warning(paste("your sleuth object is not in gene mode,",
                  "but you selected 'scaled_reads_per_base'. Selecting 'est_counts'..."))
    units <- 'scaled_reads_per_base'
  }

  units
}
