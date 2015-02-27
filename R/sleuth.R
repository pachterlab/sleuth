
#' @export
new_sleuth <- function(kal_list, sample_names, condition_names) {
  if (length(kal_list) != length(sample_names)) {
    stop(paste0("'", substitute(kal_list), "' must be the same length as '",
        substitute(sample_names), "'"))
  }

  if (length(kal_list) != length(condition_names)) {
    stop(paste0("'", substitute(kal_list), "' must be the same length as '",
        substitute(condition_names), "'"))
  }

  # append sample ane condition columns to data
  kal_list <- lapply(seq_along(kal_list), function(it)
    {
      kal_list[[it]]$abundance <- kal_list[[it]]$abundance %>%
        mutate(sample = sample_names[it], condition = condition_names[it])
      kal_list[[it]]
    })


  obs_abundance <- rbind_all(lapply(kal_list, function(k) k$abundance))


  structure(list(kal = kal_list,
      obs_abundance = obs_abundance
      ),
    class = "sleuth")
}
