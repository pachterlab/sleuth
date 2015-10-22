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

#' Mean-variance relationship
#'
#' Plot the mean-variance relationship of transcripts as modeled in \code{sleuth}. Each
#' dot represents a transcript. The blue dots represent the transcripts that
#' are used in the shrinkage estimation. The fitted curve represents the smooth fit.
#'
#' The x-axis represents the mean expression of transcripts pooled across all samples. The
#' y-axis represents the 'biological' variance after the technical variance has
#' been removed.
#'
#' @param obj a \code{sleuth} object
#' @param which_model which fit to use
#' @param point_alpha the alpha (opacity) of the points (0, 1)
#' @param point_size the size of the points
#' @param smooth_alpha the alpha (opacity) of the line
#' @param smooth_size the size of the line
#' @param smooth_color the color of the smooth line
#' @return a \code{ggplot} object
#' @export
plot_mean_var <- function(obj,
  which_model = 'full',
  point_alpha = 0.4,
  point_size = 2,
  point_colors = c('black', 'dodgerblue'),
  smooth_alpha = 1,
  smooth_size = 0.75,
  smooth_color = 'red'
  ) {

  if ( !model_exists(obj, which_model) ) {
    stop("Model '", which_model, "' does not exists. Try looking at models(",
      substitute(obj), ')')
  }

  df <- obj$fits[[which_model]]$summary

  p <- ggplot(df, aes(mean_obs, sqrt(sqrt(sigma_sq_pmax))))
  p <- p + geom_point(aes(colour = iqr),
    alpha = point_alpha, size = point_size)
  p <- p + geom_line(aes(mean_obs, sqrt(sqrt(smooth_sigma_sq))),
    size = smooth_size, alpha = smooth_alpha, colour = smooth_color)
  p <- p + scale_colour_manual(values = c("black", "dodgerblue"))
  p <- p + theme(legend.position = "none")
  p <- p + xlab("mean( log( counts + 0.5 ) )")
  p <- p + ylab("sqrt( sigma )")

  p
}

#' Plot PCA
#'
#' Plot projections of samples onto the principal components for a set of RNA-Seq experiments
#'
#' @param obj a \code{sleuth} object
#' @param pc_x integer denoting the principle component to use for the x-axis
#' @param pc_y integer denoting the principle component to use for the y-axis
#' @param use_filtered if TRUE, use filtered data. otherwise, use all data
#' @param units either 'est_counts' or 'tpm'
#' @param text_labels if TRUE, use text labels instead of points
#' @param color_by a variable to color by. if NA, then will leave all as 'black'
#' @return a gpplot object
#' @export
plot_pca <- function(obj,
  pc_x = 1L,
  pc_y = 2L,
  use_filtered = TRUE,
  units = 'est_counts',
  text_labels = FALSE,
  color_by = NULL,
  point_size = 3,
  point_alpha = 0.8,
  ...) {
  stopifnot( is(obj, 'sleuth') )

  mat <- NULL
  if (use_filtered) {
    mat <- spread_abundance_by(obj$obs_norm_filt, units)
  } else {
    mat <- spread_abundance_by(obj$obs_norm, units)
  }

  pca_res <- prcomp(mat)

  pcs <- as_df(pca_res$rotation[, c(pc_x, pc_y)])
  pcs$sample <- rownames(pcs)
  rownames(pcs) <- NULL

  pc_x <- paste0('PC', pc_x)
  pc_y <- paste0('PC', pc_y)

  pcs <- dplyr::left_join(pcs, obj$sample_to_covariates,
    by = 'sample')

  p <- NULL
  if ( !is.null( color_by ) ) {
    p <- ggplot(pcs, aes_string(pc_x, pc_y, colour = color_by))
  } else {
    p <- ggplot(pcs, aes_string(pc_x, pc_y))
  }

  if ( text_labels ) {
    p <- p + geom_text(aes(label = sample))
  } else {
    p <- p + geom_point(size = point_size, alpha = point_alpha)
  }

  p
}


#' Plot Loadings and Interpretations 
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, use filtered data. otherwise, use all data
#' @param PC  principal component to view genes contribution to that PC
#' @param gene user input on which gene and which PC's contribute the most
#' @param bool, scale or not
#' @param absolute, default true, to see all PC's magnitude (recommended)
#' @return a ggplot object
#' @export
plot_loadings <- function(obj, 
  use_filtered = TRUE,
  gene = NULL, #string please
  PC = NULL, #apparently if u type in a string or an integer (corresponding to the PC), they're both ok
  pc_count = NULL,
  bool = FALSE,
  absolute = TRUE, 
  ...) {
  stopifnot( !is.null(gene) && !is.null(PC) )# make sure proper arguments are given
  stopifnot( is(obj, 'sleuth') )

  mat <- NULL
  if (use_filtered) {
    mat <- spread_abundance_by(obj$obs_norm_filt, units)
  } else {
    mat <- spread_abundance_by(obj$obs_norm, units)
  }
  
  pca_calc <- prcomp(mat, scale = bool)
  #transpose
  loadings <- pca_calc

  #given a gene
  if (!is.null(gene)) {
    loadings <- pca_calc$x[gene,]
    if (absolute) {
      loadings <- abs(loadings)
    }
    loadings <- sort(loadings, decreasing = TRUE)
    names <- names(loadings)
  }

  #given a PC, which samples contribute the most?
  if (!is.null(PC)) {
    loadings <- pca_calc$x[,PC]
    if (absolute) {
      loadings <- abs(loadings)
    }
    loadings <- sort(loadings, decreasing = TRUE)
    names <- names(loadings)
  }

  if (!is.null(pc_count)) {
      loadings <- loadings[1:pc_count]
      names <- names[1:pc_count]
    } else {
      loadings <- loadings[1:5]
      names <- names[1:5]
    }

  dat <- data.frame(pc = names, loadings = loadings)
  dat$pc <- factor(dat$pc, levels = unique(dat$pc))

  p <- ggplot(dat, aes(x = pc, y = loadings)) 
  p <- p + geom_bar(stat = "identity")
  p <- p + xlab("Principal Components") + ylab("Contribution Scores")

  #logistics of graph
  print(is.numeric(PC))
  if (is.numeric(PC)) {
    PC <- paste0("PC ", PC)
  }

  if (!is.null(gene)) {
    p <- p + ggtitle(gene)
  } else {
    p <- p + ggtitle(PC)
  }

  p

}

#' Plot PC variances retained by percentage with option to compare specified PC
#'
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, use filtered data. otherwise, use all data
#' @param pca_number user input on how many PC to display, otherwise default is 5
#' @param bool determines scaling
#' @param (integer) PC_relative gives the option to compare subsequent principal
#' components and their contributions
#' @return a ggplot object
#' @export
plot_pc_variance <- function(obj, 
  use_filtered = TRUE,
  pca_number = NULL,
  bool = FALSE,
  PC_relative = NULL, #this is an integer
  ...) {

  mat <- NULL
  if (use_filtered) {
    mat <- spread_abundance_by(obj$obs_norm_filt, units)
  } else {
    mat <- spread_abundance_by(obj$obs_norm, units)
  }

  pca_calc <- prcomp(mat, scale = bool) #PCA calculations 

  #computation
  eigenvalues <- (pca_calc$sdev)^2  
  var <- eigenvalues*100/sum(eigenvalues)
  var2 <- eigenvalues*100/sum(eigenvalues) #because i suck at coding
  
  #from here to ....
  if (!is.null(pca_number)) {
    colsize <- pca_number
    var <- var[1:pca_number]
  } else {
    colsize <- 5 #default 5
    var <- var[1:5] #default 5
  }
  pc_asdf <- data.frame(PC_count = 1:colsize, var = var) #order here matters
  #...here is for comparison of variance for subsequent PC's

  #here is comparison per given principal component
  if(!is.null(PC_relative)) {
    pc_asdf <- data.frame(PC_count = 1:length(eigenvalues), var = var2) #because i suck at coding (var2 stupid shit)
    pc_asdf <- pc_asdf[PC_relative:nrow(pc_asdf),] #new data frame

    #if user wants to give some PCA count to graph (default is 5 or until the end)
    if (!is.null(pca_number) && (PC_relative + pca_number <= length(eigenvalues))) { #check if it does not overflow data frame
      pc_asdf <- pc_asdf[1:pca_number,] #new data frame if user wants to give a pca count number
    } else if (PC_relative + 5 >= length(eigenvalues)) {
      pc_asdf <- pc_asdf[1:nrow(pc_asdf),] 
    }
  } 

  p <- ggplot(pc_asdf, aes(x = PC_count, y = var)) + geom_bar(stat = "identity")
  p <- p + scale_x_continuous(breaks = 1:length(eigenvalues))
  p <- p + ylab("% of Variance") + xlab("Principal Components")

  p
  
}


#' Plot density
#'
#' Plot the density of a some grouping
#'
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, use filtered data. otherwise use all data
#' @param units either 'est_counts' or 'tpm'
#' @param trans a string pointing to a function to use for the transformation.
#' @param grouping a string from the columns of \code{sample_to_covariates} in
#' the sleuth object for which to group and color by
#' @param offset the offset so that transformations such as log don't compute
#' -Inf. If NULL, then will not add an offset
#' @return a \code{ggplot2} object
#' @export
plot_group_density <- function(obj,
  use_filtered = TRUE,
  units = 'est_counts',
  trans = 'log',
  grouping = setdiff(colnames(obj$sample_to_covariates), 'sample'),
  offset = 1
  ) {

  res <- kallisto_table(obj, use_filtered = use_filtered, include_covariates = TRUE)
  # res <- NULL
  # if (use_filtered) {
  #   res <- obj$obs_norm_filt
  # } else {
  #   res <- obj$obs_norm
  # }

  gdots <- list(target_id = ~target_id)
  gdots[[grouping]] <- as.formula(paste0('~', grouping))
  res <- dplyr::group_by_(res, .dots = gdots)

  mean_str <- paste0('mean(', units, ' )')
  if (!is.null(offset) && offset != 0L) {
    mean_str <- paste0(mean_str, ' + ', offset)
  }

  if (!is.null(trans)) {
    mean_str <- paste0( trans, '( ', mean_str, ' )' )
  }

  res <- dplyr::summarise_(res, .dots = setNames(mean_str, 'expression'))

  p <- ggplot(res, aes(expression))
  p <- p + geom_density(aes_string(colour = grouping, fill = grouping), alpha = 0.2)
  p <- p + xlab(mean_str)

  p
}

#' Plot sample density
#'
#' Plot the density of one particular sample
#'
#' @param obj a \code{sleuth} object

#' @param which_sample a character string matching a sample in
#' \code{obj$sample_to_covariates}
#' @param use_filtered if TRUE, use filtered data. Otherwise use all data
#' @param units either \code{'est_counts'} or \code{'tpm'}
#' @param trans a string pointing to a function to use for the transformation.
#' @param offset the offset so that transformations such as log don't compute
#' -Inf. If NULL, then will not add an offset
#' @return a \code{ggplot2} object
#' @export
plot_sample_density <- function(obj,
  which_sample = obj$sample_to_covariates$sample[1],
  use_filtered = TRUE,
  units = 'est_counts',
  trans = 'log',
  offset = 1
  ) {
  res <- NULL
  if (use_filtered) {
    res <- obj$obs_norm_filt
  } else {
    res <- obj$obs_norm
  }


  res <- dplyr::filter(res, sample == which_sample)
  trans_str <- units
  if (!is.null(offset) && offset != 0L) {
    trans_str <- paste0(trans_str, ' + ', offset)
  }
  if (!is.null(trans)) {
    trans_str <- paste0(trans, '( ', trans_str, ') ' )
  }

  p <- ggplot(res, aes_string(trans_str))
  p <- p + geom_density(fill = 'dodgerblue', alpha = 0.4)

  p
}

#' Sample to sample scatter plot
#'
#' Make a scatter plot of transcripts from two samples. to assess correlation
#'
#' @param obj a \code{sleuth} object
#' @param sample_x the string corresponding to the sample name in \code{obj$sample_to_covariates}
#' @param sample_y same as \code{sample_x} but for the y-axis
#' @param use_filtered if TRUE, use filtered data. otherwise, use all data
#' @param units either 'est_counts' or 'tpm'
#' @param offset a linear offset to help deal with zeroes if transforming the abundances
#' @param point_alpha the alpha on the points
#' @param xy_line if TRUE, plot the xy_line
#' @param xy_line_color a string denoting the color for the xy line
#' @param trans a string pointing to a function to use for the transformation.
#' This function must exist in the global namespace. This means you should be
#' able to call \code{eval('myfun')} and get a function back.
#' @param xlim a numeric vector of length two denoting the x limits
#' @param ylim same as xlim but for the y-axis
#' @return a ggplot object for the scatterplot
#' @export
plot_scatter <- function(obj,
  sample_x = obj$sample_to_covariates$sample[1],
  sample_y = obj$sample_to_covariates$sample[2],
  use_filtered = TRUE,
  units = 'est_counts',
  offset = 1,
  point_alpha = 0.2,
  xy_line = TRUE,
  xy_line_color = 'red',
  trans = 'log',
  xlim = NULL,
  ylim = NULL
  ) {

  abund <- NULL
  if (use_filtered) {
    abund <- spread_abundance_by(obj$obs_norm_filt, units)
  } else {
    abund <- spread_abundance_by(obj$obs_norm, units)
  }
  abund <- abund + offset
  abund <- as_df(abund)

  abund <- dplyr::mutate(abund, target_id = rownames(abund))

  if (!is.null(trans)) {
    sample_x <- paste0( trans, '( ', sample_x)
    sample_y <- paste0( trans, '( ', sample_y)
  }

  if ( (!is.null(offset) && !is.na(offset)) && offset != 0 ) {
    off <- deparse(eval(offset))
    sample_x <- paste0(sample_x, ' + ', off)
    sample_y <- paste0(sample_y, ' + ', off)
  }

  if (!is.null(trans)) {
    sample_x <- paste0(sample_x, ' )')
    sample_y <- paste0(sample_y, ' )')
  }

  p <- ggplot(abund, aes_string(sample_x, sample_y))
  p <- p + geom_point(alpha = point_alpha)

  if (xy_line) {
    p <- p + geom_abline(intercept = 0, slope = 1, colour = xy_line_color)
  }

  if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }

  if (!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }

  p
}

#' Plot technical variance versus observed variance
#'
#' Plot technical variance versus observed variance
#'
#' @param obj a \code{sleuth} object
#' @param which_model a character string denoting which model to use for the
#' test
#' @param point_alpha the alpha for the points
#' @param xy_line if TRUE, plot the xy_line
#' @param xy_line_color a string denoting the color for the xy line
#' @param highlight a \code{data.frame} with one column, \code{target_id}.
#' These points will be highlighted in the plot. if \code{NULL}, no points will be highlighted.
#' @param highlight_color the color to highlight points.
#' @return a \code{ggplot2} object
#' @export
plot_vars <- function(obj,
  which_beta = NULL,
  which_model = 'full',
  sig_level = 0.10,
  point_alpha = 0.2,
  sig_color = 'red',
  xy_line = TRUE,
  xy_line_color = 'red',
  highlight = NULL,
  highlight_color = 'green'
  ) {
  stopifnot( is(obj, 'sleuth') )

  cur_summary <- NULL

  if (is.null(which_beta)) {
    cur_summary <- obj$fits[[which_model]][['summary']]
    cur_summary <- dplyr::mutate(cur_summary,
      obs_var = sigma_sq + sigma_q_sq)
  } else {
    cur_summary <- sleuth_results(obj, which_beta, which_model,
      rename_cols = FALSE, show_all = FALSE)
    cur_summary <- dplyr::mutate(cur_summary,
      obs_var = sigma_sq + sigma_q_sq,
      significant = qval < sig_level)
  }

  p <- ggplot(cur_summary, aes(sqrt(obs_var), sqrt(sigma_q_sq)))

  if (is.null(which_beta)) {
    p <- p + geom_point(alpha = point_alpha)
  } else {
    p <- p + geom_point(aes(colour = significant), alpha = point_alpha)
    p <- p + scale_colour_manual(values = c('black', sig_color))
  }

  if (xy_line) {
    p <- p + geom_abline(intercept = 0, slope = 1, colour = xy_line_color)
  }

  if (!is.null(highlight)) {
    highlight_in <- data.frame()
    suppressWarnings({
      highlight_in <- dplyr::semi_join(cur_summary, highlight, by = 'target_id')
    })
    if (nrow(highlight) > 0) {
      p <- p + geom_point(aes(sqrt(obs_var), sqrt(sigma_q_sq)), data = highlight_in, colour = highlight_color)
    } else {
      warning("Couldn't find any transcripts from highlight set in this sleuth_test. They were probably filtered out.")
    }

    if ( nrow(highlight_in) > 0 && (nrow(highlight_in) != nrow(highlight)) ) {
      warning("Couldn't find any transcripts from highlight set in this sleuth_test. They were probably filtered out.")
    }
  }

  p <- p + xlab('raw standard deviation')
  p <- p + ylab('bootstrap standard deviation')

  p
}

#' MA plot
#'
#' Make an 'MA plot' for a given test. MA plots display, for each transcript, the mean of abundances across samples on the
#' x-axis and fold change on the y-axis.
#'
#' @param obj a \code{sleuth} object
#' @param which_beta a character string denoting which beta to use for
#' highlighting the transcript
#' @param which_model a character string denoting which model to use for the
#' test
#' @param point_alpha the alpha for the points
#' @return a \code{ggplot2} object
#' @export
plot_ma <- function(obj, which_beta, which_model = 'full',
  sig_level = 0.10,
  point_alpha = 0.2,
  sig_color = 'red',
  highlight = NULL,
  highlight_color = 'green'
  ) {
  stopifnot( is(obj, 'sleuth') )

  res <- sleuth_results(obj, which_beta, which_model, rename_cols = FALSE,
    show_all = FALSE)
  res <- dplyr::mutate(res, significant = qval < sig_level)

  p <- ggplot(res, aes(mean_obs, b))
  p <- p + geom_point(aes(colour = significant), alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color))
  p <- p + xlab('mean( log( counts + 0.5 ) )')
  p <- p + ylab(paste0('beta: ', which_beta))

  if (!is.null(highlight)) {
    suppressWarnings({
      highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
    })
    if (nrow(highlight) > 0) {
      p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
    } else {
      warning("Couldn't find any transcripts from highlight set in this test. They were probably filtered out.")
    }
  }

  p
}


#' Plot bootstrap summary
#'
#' Get a d
#' @param bs_df a bootstrap summary data.frame from \code{summarize_bootstrap}
#' @param ... additional arguments to `geom_point`
#' @return a ggplot2 object
#' @export
plot_bootstrap <- function(obj,
  transcript,
  units = 'est_counts',
  color_by = setdiff(colnames(obj$sample_to_covariates), 'sample'),
  x_axis_angle = 50
  ) {

  df <- get_bootstraps(obj, transcript)

  if (nrow(df) == 0) {
    stop("Couldn't find transcript ", transcript)
  }
  p <- ggplot(df, aes_string('sample', units))
  p <- p + geom_boxplot(aes_string(fill = color_by))
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle, hjust = 1))
  p <- p + ggtitle(transcript)
  p
}

#' Plot sample heatmap
#'
#' Plot sample heatmap using the Jensen-Shannon divergence
#'
#' @param obj a \code{sleuth} object
#' @param use_filtered if TRUE, use filtered data. otherwise, use everything
#' @param color_high the 'high' color
#' @param color_low the 'low' color
#' @param x_axis_angle the angle at which to put the x-axis labels
#' @return a \code{ggplot2} object
#' @export
plot_sample_heatmap <- function(obj,
  use_filtered = TRUE,
  color_high = 'white',
  color_low = 'dodgerblue',
  x_axis_angle = 50
  ) {
  abund <- NULL
  if (use_filtered) {
    abund <- spread_abundance_by(obj$obs_norm_filt, 'tpm')
  } else {
    abund <- spread_abundance_by(obj$obs_norm, 'tpm')
  }
  all_pairs <- apply_all_pairs(abund, jsd)

  all_pairs <- reshape2::melt(all_pairs, varnames = c('sample_x', 'sample_y'),
    value.name = 'jsd')

  p <- ggplot(all_pairs, aes(sample_x, sample_y))
  p <- p + geom_tile(aes(fill = jsd))
  p <- p + geom_text(aes(label = round(jsd, 3)))
  p <- p + scale_fill_gradient(high = color_high, low = color_low)
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle, hjust = 1))
  p <- p + xlab('')
  p <- p + ylab('')

  p
}

#' Plot volcano plot
#'
#' Plot a volcano plot. A volcano plot is a plot of beta value (regression coefficient)
#' vs. log(significance). Ideally, it looks like a volcano; more significance typically
#' results in higher beta
#' @param obj a  \code{sleuth} object
#' @param which_beta a character string denoting which beta to use for
#' highlighting the transcript
#' @param which_model a character string denoting which model to use for the
#' test
#' @param sig_level the significance level for Fdr
#' @param point_alpha the alpha for the points
#' @param sig_color what color to make the 'significant' transcripts
#' @param highlight a \code{data.frame} with one column, \code{target_id}.
#' These points will be displayed below in a table.
#' @return a \code{ggplot} object
#' @export


plot_volcano = function(obj, which_beta, which_model = 'full',
    sig_level = 0.10,
    point_alpha = 0.2,
    sig_color = 'red',
    highlight = NULL
    ) {
    stopifnot( is(obj, 'sleuth') )
    
    res <- sleuth_results(obj, which_beta, which_model, rename_cols = FALSE,
        show_all = FALSE)
    res <- dplyr::mutate(res, significant = qval < sig_level)
    
   
        p = ggplot(res, aes(b, -log10(qval)))
        p <- p + geom_point(aes(colour = significant), alpha = point_alpha)
        p <- p + scale_colour_manual(values = c('black', sig_color))
        p <- p + xlab('beta_value')
        p <- p + ylab('-log10(qval)')
        p <- p + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')
    
    
    p
}

#' QQ norm plot
#'
#' Create a Q-Q norm plot of the Wald statistics. The x-axis has the
#' theoretical quantile you would expect from a standard normal distribution.
#' The y-axis has the observed quantiles. It is a \code{ggplot2} version of
#' what you would get from \code{\link{qqnorm}} and \code{\link{qqline}}.
#'
#' @param obj a \code{sleuth} object
#' @param which_beta a character string denoting which beta to use for
#' highlighting the transcript
#' @param which_model a character string denoting which model to use for the
#' test
#' @param sig_level the significance level for Fdr
#' @param point_alpha the alpha for the points
#' @param sig_color what color to make the 'significant' transcripts
#' @param highlight a \code{data.frame} with one column, \code{target_id}.
#' These points will be highlighted in the plot. if \code{NULL}, no points will
#' be highlighted.
#' @param highlight_color the color to highlight points.
#' @param line_color what color to make the QQ line
#' @return a \code{ggplot2} object
#' @export
plot_qqnorm <- function(obj, which_beta, which_model = 'full',
  sig_level = 0.10,
  point_alpha = 0.2,
  sig_color = 'red',
  highlight = NULL,
  highlight_color = 'green',
  line_color = 'blue'
  ) {
  stopifnot( is(obj, 'sleuth') )

  res <- sleuth_results(obj, which_beta, which_model, rename_cols = FALSE,
    show_all = FALSE)
  res <- dplyr::mutate(res, wald_stat = b / se_b)
  res <- dplyr::filter(res, !is.na(wald_stat))
  res <- dplyr::mutate(res, significant = qval < sig_level)
  pnts <- stats::qqnorm(res$wald_stat, plot.it = FALSE)
  res <- dplyr::mutate(res, theoretical = pnts[['x']], observed = pnts[['y']])

  y <- quantile(res$observed, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))

  slope <- diff(y) / diff(x)
  intercept <- y[1L] - slope * x[1L]

  p <- ggplot(res, aes(theoretical, observed))
  p <- p + geom_point(aes(colour = significant), alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color))
  p <- p + xlab('theoretical quantile')
  p <- p + ylab(paste0('observed quantile: ', which_beta))
  p <- p + geom_abline(intercept = intercept, slope = slope, color = line_color)

  if (!is.null(highlight)) {
    suppressWarnings({
      highlight <- dplyr::semi_join(res, highlight, by = 'target_id')
    })
    if (nrow(highlight) > 0) {
      p <- p + geom_point(aes(mean_obs, b), data = highlight, colour = highlight_color)
    } else {
      warning("Couldn't find any transcripts from highlight set in this test.
        They were probably filtered out.")
    }
  }

  p
}

#' Plot clustered heatmap
#'
#' Plot a clustered heatmap. The clustering is done by the hclust function.
#'
#' @param transcripts a vector of strings containing a list of transcripts to be plotted in a heatmap
#' @param obj a \code{sleuth} object
#' @param units a string specifying which units to use, either tpm or est_counts
#' @param trans a string specifying a function to transform the data by
#' @return a \code{ggplot} object
#' @export


plot_transcript_heatmap <- function(transcripts, obj, units = 'tpm', trans = 'log')
{
    if(!all(transcripts %in% obj$obs_norm$target_id))
    {
        stop("Couldn't find the following transcripts: ", paste(transcripts[!(transcripts %in% so$obs_norm$target_id)], collapse = ", "))
    }
    
    
    tabd_df = obj$obs_norm[obj$obs_norm$target_id %in% transcripts,]
    
    if(units == 'tpm')
    {
        tabd_df = dplyr::select(tabd_df, target_id, sample, tpm)
        tabd_df = reshape2::dcast(tabd_df, target_id ~sample, value.var = 'tpm')
    }
    else if (units == 'est_counts')
    {
        tabd_df = dplyr::select(tabd_df, target_id, sample, est_counts)
        tabd_df = reshape2::dcast(tabd_df, target_id ~sample, value.var = 'est_counts')
    }
    else
    {
        stop("Didn't recognize the following unit: ", units)
    }
    
    rownames(tabd_df) = tabd_df$target_id
    tabd_df$target_id = NULL
        
    if(nchar(trans) > 0 && !is.null(trans)) {
        tFunc = eval(parse(text = trans))
        
        ggPlotExpression(as.matrix(tFunc(tabd_df)), clustRows = FALSE)
        #gplots::heatmap.2(as.matrix(tFunc(tabd_df)), Colv = FALSE, dendrogram='row', trace='none', key.xlab ='abundance', margins = c(10,30), keysize = hm_keysize, lwid = c(1,4), col = heat.colors(15))
    }
    else {
        
        ggPlotExpression(as.matrix(tabd_df), clustRows = FALSE)
        #Change the following to not rely on gplots:
        #gplots::heatmap.2(as.matrix(tabd_df), Colv = FALSE, dendrogram='row', trace='none', key.xlab ='abundance', margins = c(10,30), keysize = hm_keysize, lwid = c(1,4), col = heat.colors(15))
    }

}


#' Heatmap of expression
#'
#' Plot all of the points in an expression matrix
#'
#' @param exMat the expression matrix
#' @param clustRows if TRUE, cluster the rows by hierarchical clustering.
#' @param clustCols if TRUE, cluster the columns by hierarchical clustering.
#' @param rowNames if TRUE, print the row names on the plot
#' @param colNames if TRUE, print the column names on the plot
#' @return a ggplot object
ggPlotExpression <- function(exMat, clustRows = TRUE, clustCols = TRUE,
                             rowNames = TRUE, colNames = TRUE)
{
    if (is(exMat, 'matrix')) {
        exMat <- as.matrix(exMat)
        stopifnot(class(exMat) == 'matrix')
    }
    exMat = t(exMat)
    rowOrder <- 1:nrow(exMat)
    colOrder <- 1:ncol(exMat)
    if (clustRows)
        rowOrder <- orderByDendrogram(exMat)
    if (clustCols)
        colOrder <- orderByDendrogram(t(exMat))
    exMat <- exMat[rowOrder, colOrder]
    meltMat <- reshape2::melt(exMat, varnames = c("x", "y"))
    breaksM <- round(seq(min(meltMat$value, na.rm = T), max(meltMat$value, na.rm = T), 
                         length.out = 10), 3)
                         #print(rownames(exMat))
    if (is.null(colnames(exMat)))
        colnames(exMat) <- 1:ncol(exMat)
    meltMat$y <- factor(meltMat$y, levels = colnames(exMat))
    meltMat$x <- factor(meltMat$x, levels = rownames(exMat))
    p <- ggplot(meltMat, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradientn(colours = heat.colors(20),
                                                guide = guide_legend(title = "Expression: ",
                                                                     reverse = T, size = 14)) 
    p <- p + theme_bw() + theme(legend.text = element_text(size = 14),
                                                                   legend.title = element_text(size = 14),
                                                           legend.direction = 'vertical',
                                                           legend.position = 'top',
                                                           legend.background = element_rect(fill = "gray95", colour = "black", size = 0.5, linetype = 1),
                                                           axis.title=element_blank())
    if (rowNames)
        p <- p + theme(axis.text.x=element_text(angle = 90, size=14))
    else
        p <- p + theme(axis.text.x=element_text(size=0))

    if (colNames)
        p <- p + theme(axis.text.y=element_text(size=14))
    else
        p <- p + theme(axis.text.y=element_text(size=0))

    p
    #list(plot = p, rowOrder = rowOrder, colOrder = colOrder)
}

#' Order by dendrogram
#'
#' @param mat a matrix where the rows are observations and the columns are different dimensions on the matrix
#' @return a vector of label orderings

orderByDendrogram <- function(mat)
{
    hc <- hclust(dist(mat))
    dc <- as.dendrogram(hc)
    order.dendrogram(dc)
}
