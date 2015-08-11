#' Mean-variance relationship
#'
#' Plot the mean-variance relationship using ggvis
#'
#' @param obj a \code{sleuth} object
#' @param which_model which fit to use
#' @param point_alpha the alpha of the points (0, 1)
#' @return a \code{ggvis} object
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
#' Plot PCA for a set of RNA-Seq experiments
#'
#' @param obj a \code{sleuth} object
#' @param pc_x integer denoting which principle component to take for the x-axis
#' @param pc_y integer denoting which principle component to take for the y-axis
#' @param text_labels if TRUE, use text labels instead of points
#' @param color_by a variable to color by. if NA, then will leave all as 'black'
#' @param center center the samples before running PCA
#' @param scale scale the samples before running PCA
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

#' Interactive sleuth visualization with Shiny
#'
#' Interactive sleuth visualization with Shiny. To exit, type \code{ESC} in R.
#' @param obj a \code{sleuth} object already processed and has run
#' \code{\link{sleuth_fit}} and \code{\link{sleuth_test}}
#' @param ... additional parameters sent to ploltting functions
#' @return a \code{\link{shinyApp}} result
#' @export
#' @seealso \code{\link{sleuth_fit}}, \code{\link{sleuth_test}}
sleuth_interact <- function(obj, ...) {
  stopifnot( is(obj, 'sleuth') )
  if ( !require('shiny') ) {
    stop("'sleuth_interact()' requires 'shiny'. Please install it using
      install.packages('shiny')")
  }

  poss_covars <- dplyr::setdiff(
    colnames(obj$sample_to_covariates),
    'sample')

  p_layout <- navbarPage(
    'sleuth',
    tabPanel('analysis',
      plotOutput('mv_plt')),
    tabPanel('maps',
      fluidPage(
        selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
          selected = 1),
        selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
          selected = 2),
        selectInput('text_labels', label = 'text labels: ',
          choices = c(TRUE, FALSE), selected = TRUE),
        selectInput('color_by', label = 'color by: ',
          choices = c(NULL, poss_covars), selected = NULL),
        selectInput('center', label = 'center: ',
          choices = c(TRUE, FALSE), selected = TRUE),
        selectInput('scale', label = 'scale: ',
          choices = c(FALSE, TRUE), selected = FALSE),
        plotOutput('pca_plt')
      ))
    )

  server_fun <- function(input, output) {
    output$pca_plt <- renderPlot({

      color_by <- ifelse(is.null(input$color_by), NULL,
        as.character(input$color_by))

      plot_pca(obj,
        pc_x = as.integer(input$pc_x),
        pc_y = as.integer(input$pc_y),
        text_labels = as.logical(input$text_labels),
        color_by = color_by,
        center = as.logical(input$center),
        scale = as.logical(input$scale)
        )

    })

    output$mv_plt <- renderPlot({
      plot_mean_var(obj)
    })
  }

  shinyApp(ui = p_layout, server = server_fun)
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
