#' Merge results from expression estimation tools
#'
#' Merge the results from expression estimation tools into one data.table that
#' is easy to analyze and plot
#'
#' @param exp_list a list of expression results, each which contain columns:
#' target_id, tpm, est_counts
#' @param exp_labels a character vector of the same length as exp_list with
#' labels for each method
#' @param a data.frame with columns: target_id, counts, tpm
#' @export
merge_results <- function(exp_list, exp_labels, oracle) {
  stopifnot( is(exp_list, "list") )
  stopifnot( length(exp_list) == length(exp_labels) )

  exp_list <- lapply(seq_along(exp_list),
    function(i)
    {
      res <- exp_list[[i]] %>%
        select(target_id, tpm, est_counts) %>%
        data.table::data.table()
      data.table::setnames(res, "tpm", paste0("tpm_", exp_labels[i]))
      data.table::setnames(res, "est_counts", paste0("est_counts_", exp_labels[i]))
    })

  oracle <- oracle %>%
    rename(tpm_oracle = tpm, est_counts_oracle = counts)

  all_res <- Reduce(function(x, y) inner_join(x, y, by = c("target_id")), exp_list)


  melt_by <- function(data, unit_by) {
    m_unit <- data %>%
      select(target_id, starts_with(unit_by)) %>%
      melt(id.vars = "target_id", variable.name = "method")
    ret <- oracle %>%
      select(target_id, starts_with(unit_by)) %>%
      inner_join(m_unit, by = "target_id") %>%
      rename(estimate = value)
    setnames(ret, paste0(unit_by, "_oracle"), "oracle")
    ret
  }

  m_tpm <- melt_by(all_res, "tpm")
  m_est_counts <- melt_by(all_res, "est_counts")

  all_res <- all_res %>%
    inner_join(oracle, by = "target_id")

  structure(list(all_data = all_res, m_tpm = m_tpm, m_est_counts = m_est_counts),
    class = "merged_res")
}

#' Compute correlation of a merged results with an oracle
#'
#' Compute correlation of a merged results with an oracle
#'
#' @param mres a \code{merged_res} object from \code{merged_result}
#' @return a named list with "tpm" and "est_counts"
#' @export
compute_cor_oracle <- function(mres) {
  stopifnot( is(mres, "merged_res") )

  both_res <- lapply(list(mres$m_tpm, mres$m_est_counts),
    function(res)
    {
      res %>%
        group_by(method) %>%
        summarise(
          pearson = cor(oracle, estimate, method = "pearson"),
          spearman = cor(oracle, estimate, method = "spearman")
          )
    })

  setNames(both_res, c("tpm", "est_counts"))
}

#' Compute all pairwise correlations of a unit
#'
#' Compute all pairwise correlations of a unit
#' @param mres a \code{merged_res} object from \code{merged_result}
#' @param a character vector of length 1 that is either "tpm" or "est_counts"
#' @export
pairwise_cor <- function(mres, unit) {
  stopifnot(is(mres, "merged_res"))
  stopifnot(is(unit, "character"))
  stopifnot(length(unit) == 1)

  unit_data <- mres$all_data %>%
    select(starts_with(unit))

  pcor <- unit_data %>%
    cor(method = "pearson") %>%
    reshape2::melt(varnames = c("method_a", "method_b")) %>%
    filter(method_a != method_b) %>%
    distinct(value) %>%
    mutate(method_a = gsub(paste0(unit, "_"), "", method_a)) %>%
    mutate(method_b = gsub(paste0(unit, "_"), "", method_b)) %>%
    arrange(method_a, method_b)

  scor <- unit_data %>%
    cor(method = "spearman") %>%
    reshape2::melt(varnames = c("method_a", "method_b")) %>%
    filter(method_a != method_b) %>%
    distinct(value) %>%
    mutate(method_a = gsub(paste0(unit, "_"), "", method_a)) %>%
    mutate(method_b = gsub(paste0(unit, "_"), "", method_b)) %>%
    arrange(method_a, method_b)

  list(pearson = pcor, spearman = scor)
}

#' @export
plot_pair <- function(mres, unit, method_a, method_b) {
  stopifnot( is(mres, "merged_res") )
  stopifnot( unit == "tpm" || unit == "est_counts" )

  col_a <- paste0(unit, "_", method_a)
  col_b <- paste0(unit, "_", method_b)

  sub_data <- mres$all_data %>%
    select(target_id,
      matches(col_a),
      matches(col_b))

  ggplot(sub_data, aes_string(col_a, col_b)) +
    geom_point(alpha = 0.08) +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    stat_smooth(method = "lm", color = "blue", alpha= 0.8) +
    theme_bw()
}

#' Read eXpress data
#'
#' @param fname the file name of the express 'results.xprs' file.
#' @return a data.table
#' @export
read_xprs <- function(fname) {
    xprs <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = TRUE)
    xprs
}

#' Read kallisto output for merge_results
#'
#' Read a kallisto data into a format that place nicely with
#' \code{merge_results}. Basically equivalent to
#' \code{read_kallisto("dir")$abundance}, but doesn't do extra stuff with
#' metadata that \code{read_kallisto_does}
#'
#' @param fname the path to "expression.txt" from kallisto
#' @return a data.table
#' @export
read_kallisto_rename <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = TRUE)
    kal
}

#' @export
read_sailfish <- function(fname) {
  # TODO: fix sailfish TPM
    sf <- fread(fname, header = FALSE, skip = 5, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(sf) <- c("target_id", "length", "tpm", "rpkm", "kpkm", "est_counts")
    # sf %>%
    #     rename(tpm_sailfish = tpm, counts_sf = est_counts) %>%
    #     arrange(target_id)
    sf %>%
      select(target_id, tpm, est_counts)
}

#' @export
read_salmon <- function(fname) {
  # TODO: fix salmon TPM
    salmon <- fread(fname, header = FALSE, skip = 13, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(salmon) <- c("target_id", "length", "tpm", "fpkm", "est_counts")
    # salmon %>%
    #     rename(tpm_salmon = tpm, counts_salmon = est_counts) %>%
    #     arrange(target_id)
    salmon %>%
      select(target_id, tpm, est_counts)
}

read_kallisto_py <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    # kal %>%
    #     rename(target_id = name, tpm_kal_py = tpm, counts_kal_py = est_counts) %>%
    #     arrange(target_id)
    kal %>%
      rename(target_id = name)
}

#' @export
read_oracle <- function(fname, targ_to_eff_len) {
    oracle <- fread(paste0("sed 's/^ *//g' ", fname), header = FALSE,
        data.table = FALSE)

    targ_to_eff_len <- select(targ_to_eff_len, target_id, eff_length)

    oracle %>%
        rename(target_id = V1, counts = V2) %>%
        left_join(targ_to_eff_len, ., by = "target_id") %>%
        mutate(counts = replace(counts, is.na(counts), 0.0)) %>%
        mutate(tpm = counts / eff_length) %>%
        mutate(tpm = replace(tpm, eff_length == 0, 0.0)) %>%
        mutate(tpm = 1e6 * tpm / sum(tpm)) %>%
        rename(tpm_oracle = tpm) %>%
        arrange(target_id)
}

