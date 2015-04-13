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
      reshape2::melt(id.vars = "target_id", variable.name = "method")
    ret <- data.table::data.table(oracle) %>%
      select(target_id, starts_with(unit_by)) %>%
      inner_join(data.table::data.table(m_unit), by = "target_id") %>%
      rename(estimate = value)
    data.table::setnames(ret, paste0(unit_by, "_oracle"), "oracle")
    ret
  }

  m_tpm <- melt_by(all_res, "tpm")
  m_est_counts <- melt_by(all_res, "est_counts")

  all_res <- all_res %>%
    inner_join(data.table::data.table(oracle), by = "target_id")

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
          spearman = cor(oracle, estimate, method = "spearman"),
          med_scaled_err = median(abs(scaled_error(estimate, oracle)),
              na.rm = TRUE))
    })

  setNames(both_res, c("tpm", "est_counts"))
}

#' @export
filtered_no_summary <- function(mres, filter_exp) {
  stopifnot( is(mres, "merged_res") )
  do_filter <- if (missing(filter_exp)) {
    FALSE
  } else {
    filter_exp <- deparse(substitute(filter_exp))
    filtered_ids <- mres$all_data %>%
      filter_(.dots = list(filter_exp)) %>%
      select(target_id)
    TRUE
  }

  both_res <- lapply(list(mres$m_tpm, mres$m_est_counts),
    function(res)
    {
      if (do_filter) {
        res <- data.table(res) %>%
          inner_join(data.table(filtered_ids), by = c("target_id"))
      }

      res %>%
        group_by(method) %>%
        mutate(
          scaled_err = abs(scaled_error(estimate, oracle)),
          per_err = abs(percent_error(estimate, oracle))
          )
    })

  setNames(both_res, c("tpm", "est_counts"))
}

#' @export
filtered_summary <- function(mres, filter_exp) {
  stopifnot( is(mres, "merged_res") )
  do_filter <- if (missing(filter_exp)) {
    FALSE
  } else {
    filter_exp <- deparse(substitute(filter_exp))
    filtered_ids <- mres$all_data %>%
      filter_(.dots = list(filter_exp)) %>%
      select(target_id)
    TRUE
  }

  both_res <- lapply(list(mres$m_tpm, mres$m_est_counts),
    function(res)
    {
      if (do_filter) {
        res <- data.table(res) %>%
          inner_join(data.table(filtered_ids), by = c("target_id"))
      }

      res %>%
        group_by(method) %>%
        summarise(
          pearson = cor(estimate, oracle, method = "pearson"),
          spearman = cor(estimate, oracle, method = "spearman"),
          med_scaled_err = median(abs(scaled_error(estimate, oracle)),
              na.rm = TRUE),
          med_per_err = median(abs(percent_error(estimate, oracle)))
          )
    })

  setNames(both_res, c("tpm", "est_counts"))
}

scaled_error <- function(estimate, truth) {
  2 *(estimate - truth)  / (estimate + truth)
}

percent_error <- function(estimate, truth) {
  (estimate - truth) / truth
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
    mutate(ah = ifelse(as.character(method_a) < as.character(method_b), paste0(method_a, method_b), paste0(method_b, method_a))) %>%
    distinct(ah) %>%
    select(-ah) %>%
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

#' @export
plot_ranks <- function(mres, unit, method) {
  stopifnot( is(mres, "merged_res") )
  stopifnot( unit == "tpm" || unit == "est_counts" )

  col_meth <- paste0(unit, "_", method)
  col_oracle <- paste0(unit, "_", "oracle")
  sub_data <- mres$all_data %>%
    select_("target_id", col_meth, col_oracle) %>%
    as.data.frame()
  sub_data$meth_rank <- rank(sub_data[,col_meth])
  sub_data$oracle_rank <- rank(sub_data[,col_oracle])
  sub_data$d <- with(sub_data, oracle_rank - meth_rank)
  # sub_data <- as.data.frame(sub_data, stringsAsFactors = FALSE) %>%
  #   mutate(meth_rank = lazyeval::interp(~rank(x), x = col_meth))

  # ggplot(sub_data, aes_string("oracle_rank", "d")) +
  #   geom_point(alpha = 0.2)
  # print(head(sub_data))
  plt <- ggvis::ggvis(sub_data, x = ~oracle_rank, y = ~d) %>%
    ggvis::layer_points(fill := "black", opacity := 0.2) %>%
    ggvis::add_tooltip(function(dat){
      print(dat)
      print(sub_data[(sub_data$oracle_rank == dat$oracle_rank) & (sub_data$d == dat$d),])
      }, "hover")

  list(plt = plt, data = sub_data)
  plt
}

#' @export
plot_pe <- function(mres, unit, method) {
  stopifnot( is(mres, "merged_res") )
  stopifnot( unit == "tpm" || unit == "est_counts" )

  col_meth <- paste0(unit, "_", method)
  col_oracle <- paste0(unit, "_", "oracle")
  sub_data <- mres$all_data %>%
    select_("target_id", col_meth, col_oracle) %>%
    as.data.frame()

  sub_data$pe <- pe(sub_data[,col_meth], sub_data[,col_oracle])
  ggplot(sub_data, aes_string(col_oracle, "pe")) +
    geom_point(alpha = 0.02)
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
    # kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
    #     data.table = TRUE)
    # kal
  read.table(fname, stringsAsFactors=FALSE, header = TRUE)
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

#' @export
read_rsem <- function(fname) {
  read.table(fname, header = TRUE, stringsAsFactors = FALSE) %>%
    rename(target_id = transcript_id,
      est_counts = expected_count,
      tpm = TPM, eff_len = effective_length) %>%
    select(-gene_id, -FPKM, -IsoPct)
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


#' @export
read_cufflinks <- function(fname, mean_frag_len) {
  data <- data.table::fread(fname, header = TRUE)

  data <- data %>%
    mutate(tpm = FPKM * 1e6/ sum(FPKM),
      est_counts = tpm_to_alpha(tpm, length - mean_frag_len))

  data %>%
    select(target_id = tracking_id, tpm, est_counts)
}


tpm_to_alpha <- function(tpm, eff_len) {
  stopifnot(length(tpm) == length(eff_len))

  alpha <- (tpm * eff_len) / 1e6
  alpha <- alpha / sum(alpha)

  alpha
}

#' @export
alpha_load <- function(dir_name, oracle_targs) {
  alpha <- data.table::fread(file.path(dir_name, "alpha.txt"), data.table = FALSE, header = TRUE)
  alpha <- t(as.matrix(alpha))
  alpha <- alpha[oracle_targs,]

  ll <- data.table::fread(file.path(dir_name, "ll.txt"), header = FALSE, data.table = FALSE)

  list(alpha = alpha, ll = ll[,1])
}

#' @export
alpha_to_oracle_cor <- function(alpha, oracle_counts) {
  all_spearman <- lapply(1:ncol(alpha),
    function(i)
    {
      if (i %% 100 == 0) print(i)
      cor(alpha[,i], oracle_counts, method = "spearman")
    })

  unlist(all_spearman)
}

#' Simulate counts table
#'
#' Simulate a counts table, given TPM and effective length
#'
#' @param mixture a data.frame containing columns "tpm", "target_id" and "eff_len"
#' @param total_counts the total number of reads
