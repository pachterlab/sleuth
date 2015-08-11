#' Merge results from differential expression tools
#'
#' Merge results from differential expression tools to benchmark calls. Plots and stuff can be made from this object as well
#'
#' @param de_list a list of de results
#' @param de_label a character vector of the same length as de_list with labels for each method
#' @param oracle a data.frame with the columns target_id, is_de and other optional columns (TBD)
#' @export
new_de_benchmark <- function(de_list, de_labels, oracle) {
  stopifnot( is(de_list, "list") )
  stopifnot( length(de_list) == length(de_labels) )

  # extract relevant columns and rename
  de_list <- lapply(seq_along(de_list),
    function(i)
    {
      res <- de_list[[i]] %>%
        select(target_id, pval, qval) %>%
        data.table::data.table()

      c('pval', 'qval') %>%
        data.table::setnames(res, ., paste0(., '_', de_labels[i]))

      res
    })

  all_res <- Reduce(function(x, y) dplyr::inner_join(x, y, by = c("target_id")),
    de_list)

  oracle <- oracle %>%
    select(target_id, is_de, log_fc)

  melt_by <- function(data, unit_by) {
    m_unit <- data %>%
      select(target_id, starts_with(unit_by)) %>%
      reshape2::melt(id.vars = "target_id", variable.name = "method")
    ret <- data.table::data.table(oracle) %>%
      inner_join(data.table::data.table(m_unit), by = "target_id") %>%
      rename(estimate = value)
    # data.table::setnames(ret, paste0(unit_by, "_oracle"), "oracle")
    ret
  }
  m_pval <- melt_by(all_res, "pval")
  m_qval <- melt_by(all_res, "qval")

  all_res <- all_res %>%
    inner_join(data.table::data.table(oracle), by = "target_id")

  ret <- list(all_data = all_res,
    m_pval = m_pval,
    m_qval = m_qval,
    labels = de_labels)
  class(ret) <- "de_benchmark"
  ret
}

#' @export
fdr_tpr_plot <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

  de_bench$m_qval %>%
    group_by(method) %>%
    arrange(estimate) %>%
    mutate(tpr = cummean(is_de)) %>%
    ggplot(aes(estimate, tpr, group = method)) +
      geom_line(aes(colour = method)) +
      xlab("FDR") +
      ylab("TPR") +
      xlim(0, 1) +
      ylim(0, 1)
}


#' @export
fdr_nde_plot <- function(de_bench, estimate = TRUE) {
  stopifnot( is(de_bench, "de_benchmark") )

  n_true_de <- sum(de_bench$all_data$is_de)
  print(n_true_de)

  pvals <- mutate(de_bench$m_pval, method = sub("pval_", "", method))
  qvals <- select(de_bench$m_qval, target_id, method, estimate)
  qvals <- mutate(qvals, method = sub("qval_", "", method))
  qvals <- rename(qvals, qval = estimate)

  pvals <- inner_join(pvals, qvals, by = c("target_id", "method"))

  plt <- pvals %>%
    mutate(method = sub("pval_", "", method)) %>%
    group_by(method) %>%
    arrange(estimate) %>%
    mutate(nde = 1:n(), tFDR = cummean(!is_de)) %>%
    ggplot(aes(nde, estimate, group = method))
  if (estimate) {
    plt <- plt + geom_line(aes(nde, qval, colour = method), linetype = 3)
  }

  plt <- plt +
      geom_line(aes(nde, tFDR, colour = method, linetype = method), size = 0.8,
        alpha = 0.8) +
      geom_vline(xintercept = n_true_de, linetype = 3) +
      geom_hline(yintercept = 0.10, linetype = 3) +
      #theme_bw() +
      xlab("Number of features called DE") +
      ylab("FDR") +
      ylim(0, 1)

  plt
}

#' @export
pval_distribution <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

  de_bench$m_pval %>%
    ggplot(aes(estimate, y = ..density..)) +
    geom_histogram(binwidth = 0.02) +
    facet_wrap( ~ method )
}

#' @export
de_rank_scatter <- function(de_bench, cutoff = 5000) {
  stopifnot( is(de_bench, "de_benchmark") )
  rank_diff <- function(x) {
    x$true_rank <- rank(-abs(x$log_fc))
    x$est_rank <- rank(x$estimate)
    tmp <- head(arrange(x, est_rank), cutoff)
    mutate(tmp, relative_rank = rank(true_rank))
  }

  grped <- as.data.frame(de_bench$m_qval) %>%
    group_by(method)

  tmp <- grped %>%
    filter(method == 'qval_DESeq2') %>%
    ungroup()

  do(grped, rank_diff(.))
}
