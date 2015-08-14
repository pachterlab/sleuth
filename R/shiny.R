#' Interactive sleuth visualization with Shiny
#'
#' Interactive sleuth visualization with Shiny. To exit, type \code{ESC} in R.
#'
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
  samp_names <- obj$sample_to_covariates[['sample']]
  poss_models <- names(models(obj))
  cat('these are the samp names:', samp_names, '\n')

  p_layout <- navbarPage(
    'sleuth',

    tabPanel('differential analysis',
      fluidRow(
        column(1,
          numericInput('max_fdr', label = 'Fdr cutoff:', value = 0.10,
            min = 0, max = 1, step = 0.01)),
        column(3,
          selectInput('which_model', label = 'fit: ',
            choices = poss_models,
            selected = poss_models[1])
          ),
        column(2,
          uiOutput('which_beta_ctrl')
            ),
        column(1,
          numericInput('ma_alpha', label = 'point alpha:', value = 0.2,
            min = 0, max = 1, step = 0.01))
        ),
      plotOutput('ma')),

    tabPanel('analysis',
      fluidRow(
        column(3,
          selectInput('sample_x', label = 'x-axis: ',
            choices = samp_names,
            selected = samp_names[1])
          ),
        column(3,
          selectInput('sample_y', label = 'y-axis: ',
            choices = samp_names,
            selected = samp_names[2])
          ),
        column(1,
          textInput('trans', label = 'transformation: ',
            value = 'log')),
        column(1,
          numericInput('scatter_alpha', label = 'point alpha:', value = 0.2,
            min = 0, max = 1, step = 0.01))
        ),
      plotOutput('scatter')),

    tabPanel('diagnostics',
      plotOutput('mv_plt')),

    tabPanel('maps',
      fluidRow(
        column(1,
          selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
            selected = 1)
          ),
        column(1,
          selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
            selected = 2)
          ),
        column(1,
          selectInput('text_labels', label = 'text labels: ',
            choices = c(TRUE, FALSE), selected = TRUE)
          ),
        column(4,
          selectInput('color_by', label = 'color by: ',
            choices = c(NULL, poss_covars), selected = NULL)
          ),
        column(1,
          selectInput('center', label = 'center: ',
            choices = c(TRUE, FALSE), selected = TRUE)
          ),
        column(1,
          selectInput('scale', label = 'scale: ',
          choices = c(FALSE, TRUE), selected = FALSE))),
        fluidRow(plotOutput('pca_plt'))
      ))

  server_fun <- function(input, output) {

    output$scatter <- renderPlot({
      plot_scatter(obj, input$sample_x, input$sample_y,
        trans = input$trans, point_alpha = input$scatter_alpha)
    })

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

    output$which_beta_ctrl <- renderUI({
      cat('hi: ', input$which_model, '\n')
      poss_tests <- tests(models(obj)[[input$which_model]])
      print(poss_tests)
      selectInput('which_beta', 'beta: ', choices = poss_tests)
    })

    output$ma <- renderPlot({
      cat('which_beta: ', input$which_beta, '\n')
      val <- input$which_beta
      if ( is.null(val) ) {
        poss_tests <- tests(models(obj)[[input$which_model]])
        val <- poss_tests[1]
      }
      plot_ma(obj, val, input$which_model, sig_level = input$max_fdr)
    })
  }

  shinyApp(ui = p_layout, server = server_fun)
}
