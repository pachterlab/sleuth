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
    a('sleuth', href = 'http://pimentel.github.io/sleuth', style = 'color: black;'),

    tabPanel('summaries',
      ####
      dataTableOutput('summary_dt')
      ),

    tabPanel('maps',

      ####
      fluidRow(
        column(3,
          selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
            selected = 1)
          ),
        column(3,
          selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
            selected = 2)
          ),
        column(2,
          selectInput('text_labels', label = 'text labels: ',
            choices = c(TRUE, FALSE), selected = TRUE)
          ),
        column(4,
          selectInput('color_by', label = 'color by: ',
            choices = c(NULL, poss_covars), selected = NULL)
          )
        ),
      fluidRow(
        column(2,
          selectInput('pca_units', label = 'units: ',
            choices = c('est_counts', 'tpm'),
            selected = 'est_counts')),
        column(2,
          checkboxInput('pca_filt', label = 'filter: ',
            value = TRUE))
        ),
        fluidRow(plotOutput('pca_plt'))
      ),

    navbarMenu('analyses',

      ####
      tabPanel('MA plots',
        fluidRow(
          column(2,
            numericInput('max_fdr', label = 'max Fdr:', value = 0.10,
              min = 0, max = 1, step = 0.01)),
          column(4,
            selectInput('which_model', label = 'fit: ',
              choices = poss_models,
              selected = poss_models[1])
            ),
          column(4,
            uiOutput('which_beta_ctrl')
            ),
          column(2,
            numericInput('ma_alpha', label = 'opacity:', value = 0.2,
              min = 0, max = 1, step = 0.01))
          ),
        fluidRow(plotOutput('ma', click = 'ma_click', hover = 'ma_hover',
            brush = 'ma_brush')),
        fluidRow(plotOutput('vars')),
        fluidRow(verbatimTextOutput('ma_brush_out'))
        ),

      ####
      tabPanel('mean-variance plot',
        plotOutput('mv_plt')
        ),

      ####
      tabPanel('scatter plots',
        fluidRow(
          column(4,
            selectInput('sample_x', label = 'x-axis: ',
              choices = samp_names,
              selected = samp_names[1])
            ),
          column(4,
            selectInput('sample_y', label = 'y-axis: ',
              choices = samp_names,
              selected = samp_names[2])
            ),
          column(2,
            textInput('trans', label = 'transform: ',
              value = 'log')),
          column(2,
            numericInput('scatter_alpha', label = 'opacity:', value = 0.2,
              min = 0, max = 1, step = 0.01))
          ),
        fluidRow(
          column(2,
            selectInput('scatter_units', label = 'units: ',
              choices = c('est_counts', 'tpm'),
              selected = 'est_counts')),
          column(2,
            checkboxInput('scatter_filt', label = 'filter: ',
              value = TRUE))
          ),
        plotOutput('scatter'))
      )
    ) # navbarPage

  server_fun <- function(input, output) {

    output$summary_dt <- renderDataTable(summary(obj))

    ###
    output$scatter <- renderPlot({
      plot_scatter(obj, input$sample_x, input$sample_y,
        trans = input$trans, point_alpha = input$scatter_alpha,
        units = input$scatter_units,
        use_filtered = input$scatter_filt)
    })

    ###
    output$pca_plt <- renderPlot({

      color_by <- ifelse(is.null(input$color_by), NULL,
        as.character(input$color_by))

      plot_pca(obj,
        pc_x = as.integer(input$pc_x),
        pc_y = as.integer(input$pc_y),
        text_labels = as.logical(input$text_labels),
        color_by = color_by,
        use_filtered = input$pca_filt,
        units = input$pca_units
        )

    })

    ###
    output$mv_plt <- renderPlot({
      plot_mean_var(obj)
    })

    ### MA
    active_ma <- reactiveValues(highlight = data.frame(target_id = character()))

    output$which_beta_ctrl <- renderUI({
      poss_tests <- tests(models(obj)[[input$which_model]])
      print(poss_tests)
      selectInput('which_beta', 'beta: ', choices = poss_tests)
    })

    output$ma <- renderPlot({
      val <- input$which_beta
      if ( is.null(val) ) {
        poss_tests <- tests(models(obj)[[input$which_model]])
        val <- poss_tests[1]
      }
      plot_ma(obj, val, input$which_model, sig_level = input$max_fdr,
        point_alpha = input$ma_alpha)
    })

    output$vars <- renderPlot({
      plot_vars(obj,
        which_model = input$which_model,
        point_alpha = input$ma_alpha,
        highlight = active_ma$highlight
        )
    })

    output$ma_brush_out <- renderPrint({
      wb <- input$which_beta
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj)[[input$which_model]])
        wb <- poss_tests[1]
      }
      res <- sleuth_results(obj, wb, input$which_model)
      if (!is.null(input$ma_brush)) {
        res <- enclosed_brush(res, input$ma_brush)
      } else {
        res <- NULL
      }
      active_ma$highlight <- res

      res
    })
  }

  shinyApp(ui = p_layout, server = server_fun)
}

enclosed_brush <- function(df, brush) {
  xvar <- brush$mapping$x
  yvar <- brush$mapping$y

  xbool <- brush$xmin <= df[[xvar]] & df[[xvar]] <= brush$xmax
  ybool <- brush$ymin <= df[[yvar]] & df[[yvar]] <= brush$ymax

  df[xbool & ybool,]
}
