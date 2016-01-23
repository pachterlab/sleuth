propagate_transcript_filter <- function(filter_df, target_mapping,
  grouping_column) {

  filtered_target_mapping <- dplyr::inner_join(as.data.table(filter_df),
    as.data.table(target_mapping), by = 'target_id')

  filtered_target_mapping <- dplyr::select_(filtered_target_mapping,
    grouping_column)

  data.table::setnames(filtered_target_mapping, grouping_column, 'target_id')
  filtered_target_mapping <- dplyr::distinct(filtered_target_mapping)

  filtered_target_mapping
}
