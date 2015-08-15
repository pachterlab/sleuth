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
#' @param text_labels if TRUE, use text labels instead of points
#' @param color_by a variable to color by. if NA, then will leave all as 'black'
#' @param ... additional arguments passed to \code{\link{geom_point}} or
#' \code{\link{geom_text}}
#' @return a gpplot object
#' @export
plot_pca <- function(obj,
  pc_x = 1L,
  pc_y = 2L,
  text_labels = FALSE,
  color_by = NULL,
  center = TRUE,
  scale = FALSE,
  ...) {
  stopifnot( is(obj, 'sleuth') )

  mat <- tidyr::spread(
    dplyr::select(obj$obs_norm, target_id, sample, est_counts),
    sample,
    est_counts)
  rownames(mat) <- mat$target_id
  mat$target_id <- NULL
  mat <- as.matrix(mat)

  pca_res <- prcomp(mat, center = center, scale = scale)

  pcs <- sleuth:::as_df(pca_res$rotation[, c(pc_x, pc_y)])
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
    p <- p + geom_point()
  }

  p
}

#' Sample to sample scatter plot
#'
#' Make a scatter plot of transcripts from two samples. to assess correlation
#'
#' @param obj a \code{sleuth} object
#' @param sample_x the string corresponding to the sample name in \code{obj$sample_to_covariates}
#' @param sample_y same as \code{sample_x} but for the y-axis
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
  offset = 1,
  point_alpha = 0.2,
  xy_line = TRUE,
  xy_line_color = 'red',
  trans = 'log',
  xlim = NULL,
  ylim = NULL) {

  abund <- spread_abundance_by(obj$obs_norm, 'est_counts')
  abund <- abund + offset
  abund <- as_df(abund)

  abund <- dplyr::mutate(abund, target_id = rownames(abund))

  if (!is.null(trans)) {
    sample_x <- paste0( trans, '( ', sample_x)
    sample_y <- paste0( trans, '( ', sample_y)
  }

  if ( offset != 0 ) {
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
  sig_color = 'red'
  ) {
  stopifnot( is(obj, 'sleuth') )

  res <- sleuth_results(obj, which_beta, which_model)
  res <- dplyr::mutate(res, significant = ifelse( qval < sig_level, TRUE, FALSE ))

  p <- ggplot(res, aes(mean_obs, b))
  p <- p + geom_point(aes(colour = significant), alpha = point_alpha)
  p <- p + scale_colour_manual(values = c('black', sig_color))
  p <- p + xlab('mean( log( counts + 0.5 ) )')
  p <- p + ylab(paste0('beta: ', which_beta))

  p
}

#' Plot bootstrap summary
#'
#' Get a d
#' @param bs_df a bootstrap summary data.frame from \code{summarize_bootstrap}
#' @param ... additional arguments to `geom_point`
#' @return a ggplot2 object
#' @export
plot_bootstrap <- function(bs_df, ...)
{
    mean_col <- grep("mean_", colnames(bs_df), value = TRUE)
    cv_col <- grep("cv_", colnames(bs_df), value = TRUE)

    ggplot(bs_df, aes_string(mean_col, cv_col)) +
        geom_point(...)
}

#' @export
prep_bs_plots <- function(obj) {
  stopifnot( is(obj, "sleuth") )

  bs <- lapply(obj[["kal"]], dcast_bootstrap, "est_counts")
  obj$dcast_bs <- bs

  obj
}

#' @export
plot_transcript_variability <- function(obj, target_id, group_string = NULL) {
  target_exp <- lapply(seq_along(obj$dcast_bs),
    function(i)
    {
      m_est <- reshape2::melt(obj$dcast_bs[[i]][target_id,],
        value.name = "est_counts")
      suppressWarnings(cbind(m_est, obj$sample_to_condition[i,]))
    }) %>%
      rbind_all()

  plt <- ggplot(target_exp, aes(sample, est_counts))

  if (is.null(group_string)) {
    plt <- plt + geom_boxplot()
  } else {
    plt <- plt + geom_boxplot(aes_string(fill = group_string))
  }

  if (!is(target_id, "character")) {
    target_id <- rownames(obj$dcast_bs[[1]])[target_id]
  }

  list(plt = plt + ggtitle(target_id), target_id = target_id)
}

#' @export
plot_transcript <- function(obj, trans_name, group_string = NULL) {
  stopifnot( is(obj, "sleuth") )

  if ( !is(trans_name, "character") ) {
    stopifnot( is(trans_name, "numeric") || is(trans_name, "integer"))
    trans_name <- rownames(obj$dcast_bs[[1]])[trans_name]
  }

  target_exp <- obj$obs_norm %>%
    filter(target_id == trans_name)

  plt <- ggplot(target_exp, aes(sample, est_counts))

  if (is.null(group_string)) {
    plt <- plt + geom_boxplot()
  } else {
    plt <- plt + geom_boxplot(aes_string(fill = group_string))
  }

  plt + ggtitle(trans_name)
}
