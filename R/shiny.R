#' Interactive sleuth visualization with Shiny
#'
#' Interactive sleuth visualization with Shiny. To exit, type \code{ESC} in R.
#'
#' @param obj a \code{sleuth} object already processed and has run
#' \code{\link{sleuth_fit}} and \code{\link{sleuth_test}}
#' @param ... additional parameters sent to plotting functions
#' @return a \code{\link{shinyApp}} result
#' @export
#' @seealso \code{\link{sleuth_fit}}, \code{\link{sleuth_test}}
sleuth_live <- function(obj, ...) {
  stopifnot( is(obj, 'sleuth') )
  if ( !require('shiny') ) {
    stop("'sleuth_interact()' requires 'shiny'. Please install it using
      install.packages('shiny')")
  }

  poss_covars <- dplyr::setdiff(
    colnames(obj$sample_to_covariates),
    'sample')
  samp_names <- obj$sample_to_covariates[['sample']]
  poss_models <- names(models(obj, verbose = FALSE))

  p_layout <- navbarPage(
    a('sleuth', href = 'http://pachterlab.github.io/sleuth', target = '_blank',
      style = 'color: black;'),
    windowTitle = 'sleuth',

    tabPanel('diagnostics',
      ####
      fluidRow(
        column(12,
          p(h3('scatter plot '), "Display scatter plot for any two samples and then select a set of transcripts to explore their variance across samples.")
          ),
          offset = 1),
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
            numericInput('scatter_offset', label = 'offset: ', value = 1))
          ),
        fluidRow(
          column(2,
            selectInput('scatter_units', label = 'units: ',
              choices = c('est_counts', 'tpm'),
              selected = 'est_counts')),
          column(2,
            checkboxInput('scatter_filt', label = 'filter',
              value = TRUE)),
          column(2,
            numericInput('scatter_alpha', label = 'opacity:', value = 0.2,
              min = 0, max = 1, step = 0.01))
          ),
        fluidRow(plotOutput('scatter', brush = 'scatter_brush')),
        fluidRow(plotOutput('scatter_vars')),
        fluidRow(dataTableOutput('scatter_brush_table'))
      ),

    navbarMenu('summaries',
      ####
      tabPanel('experimental data',
      fluidRow(
        column(12,
          p(h3('experimental data'), "Names of samples, number of mapped reads, number of boostraps performed by kallisto, and sample to covariate mappings.")
          ),
          offset = 1),
        fluidRow(dataTableOutput('summary_dt'))
        ),

      ####
      tabPanel('densities',
      fluidRow(
        column(12,
          p(h3('distribution of abundances'), "Distributions of abundances of individual samples or groupings by covariates.")
          ),
          offset = 1),
        fluidRow(
          column(4,
            selectInput('cond_dens_grp', 'grouping: ',
              choices = poss_covars,
              selected = poss_covars[1])
            ),
          column(2,
            selectInput('cond_dens_units', label = 'units: ',
              choices = c('tpm', 'est_counts'),
              selected = 'tpm')),
          column(2,
            checkboxInput('cond_dens_filt', label = 'filter',
              value = TRUE)),
          column(2,
            textInput('cond_dens_trans', label = 'transform: ',
              value = 'log')),
          column(2,
            numericInput('cond_dens_offset', label = 'offset: ', value = 1))
          ),
        fluidRow(plotOutput('condition_density')),
        fluidRow(
          column(4,
            selectInput('samp_dens', 'sample: ',
              choices = samp_names,
              selected = samp_names[1]))),
        fluidRow(plotOutput('sample_density'))
        )
      ),

    navbarMenu('maps',
      ###
      tabPanel('sample heatmap',
      fluidRow(
        column(12,
          p(h3('sample heatmap'), "Jensen-Shannon divergence between pairs of samples.")
          ),
          offset = 1),
        fluidRow(checkboxInput('samp_heat_filt', label = 'filter', value = TRUE)),
        fluidRow(plotOutput('samp_heat_plt'))
        ),
      ####
      tabPanel('PCA',
      fluidRow(
        column(12,
          p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components.")
          ),
          offset = 1),
        fluidRow(
          column(3,
            selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
              selected = 1)
            ),
          column(3,
            selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
              selected = 2)
            ),
          column(4,
            selectInput('color_by', label = 'color by: ',
              choices = c(NULL, poss_covars), selected = NULL)
            ),
          column(2,
            numericInput('pca_point_size', label = 'size: ', value = 3))
          ),
        fluidRow(
          column(2,
            selectInput('pca_units', label = 'units: ',
              choices = c('est_counts', 'tpm'),
              selected = 'est_counts')),
          column(3,
            checkboxInput('pca_filt', label = 'filter',
              value = TRUE),
            checkboxInput('text_labels', label = 'text labels',
              value = TRUE)
            )
          ),
        fluidRow(plotOutput('pca_plt'))
        )
      ),

    navbarMenu('analyses',

      ####
      tabPanel('MA plot',
      fluidRow(
        column(12,
          p(h3('MA plot'), "Plot of abundance versus fixed effect (e.g. fold change). Select a set of transcripts to explore their variance across samples. ")
          ),
          offset = 1),
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
        # fluidRow(column(8,
        #     selectizeInput('ma_trans', 'highlight transcripts: ',
        #       choices = obj$filter_df[['target_id']],
        #       multiple = TRUE, width = '100%'))
        #   ),
        fluidRow(plotOutput('ma', brush = 'ma_brush')),
        #fluidRow(plotOutput('vars', brush = 'vars_brush')),
        fluidRow(plotOutput('vars')),
        fluidRow(dataTableOutput('ma_brush_out'))
        ),

      ####
      tabPanel('mean-variance plot',
      fluidRow(
        column(12,
          p(h3('mean-variance plot'), "Plot of abundance versus square root of standard deviation which is used for shrinkage estimation. The blue dots are in the interquartile range and the red curve is the fit used by sleuth." )
          ),
          offset = 1),
        fluidRow(plotOutput('mv_plt'))
        ),


      ####
      tabPanel('transcript table',
      fluidRow(
        column(12,
          p(h3('transcript table'), "Table of transcript names, gene names (if supplied), sleuth parameter estimates, tests, and summary statistics." )
          ),
          offset = 1),
        fluidRow(
          column(4,
            selectInput('which_model_de', label = 'fit: ',
              choices = poss_models,
              selected = poss_models[1])
            ),
          column(4,
            uiOutput('which_beta_ctrl_de')
            )
          ),
        dataTableOutput('de_dt')
        ),

      tabPanel('transcript view',
      fluidRow(
        column(12,
          p(h3('transcript view'), "Boxplots of transcript abundances showing technical variation in each sample." )
          ),
          offset = 1),
        fluidRow(column(4,
            textInput('bs_var_input', label = 'transcript: ', value = '')
            ),
          column(4,
            selectInput('bs_var_color_by', label = 'color by: ',
              choices = c(NULL, poss_covars), selected = NULL)
            ),
          column(3,
            selectInput('bs_var_units', label = 'units: ',
              choices = c('est_counts', 'tpm'),
              selected = 'est_counts'))
          ),
          fluidRow(HTML('&nbsp;&nbsp;&nbsp;'), actionButton('bs_go', 'view')),
          fluidRow(plotOutput('bs_var_plt'))
          )
      )
    ) # navbarPage

  server_fun <- function(input, output) {

    output$summary_dt <- renderDataTable(summary(obj))

    output$condition_density <- renderPlot({
      plot_group_density(obj,
        grouping = input$cond_dens_grp,
        units = input$cond_dens_units,
        use_filtered = input$cond_dens_filt,
        trans = input$cond_dens_trans,
        offset = input$cond_dens_offset)
    })

    output$sample_density <- renderPlot({
      plot_sample_density(obj,
        which_sample = input$samp_dens,
        units = input$cond_dens_units,
        use_filtered = input$cond_dens_filt,
        trans = input$cond_dens_trans,
        offset = input$cond_dens_offset
        )
    })

    ###
    rv_scatter <- reactiveValues(highlight = NULL, data = NULL)

    output$scatter <- renderPlot({
      p <- plot_scatter(obj, input$sample_x, input$sample_y,
        trans = input$trans, point_alpha = input$scatter_alpha,
        units = input$scatter_units,
        use_filtered = input$scatter_filt,
        offset = as.numeric(input$scatter_offset))
      # get the data in the form that it is displayed in the plot
      x <- eval(p$mapping$x, envir = p$data)
      y <- eval(p$mapping$y, envir = p$data)
      rv_scatter$data <- data.frame(target_id = p$data$target_id, x = x, y = y,
        stringsAsFactors = FALSE)

      p
    })

    output$scatter_vars <- renderPlot({
      wb <- input$which_beta
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
        wb <- poss_tests[1]
      }
      plot_vars(obj,
        which_beta = wb,
        which_model = input$which_model,
        point_alpha = input$scatter_alpha,
        highlight = rv_scatter$highlight_vars,
        sig_level = input$max_fdr
        )
    })


    output$scatter_brush_table <- renderDataTable({
      res <- NULL
      wb <- input$which_beta
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
        wb <- poss_tests[1]
      }
      sr <- sleuth_results(obj, wb, input$which_model, rename_cols = FALSE,
        show_all = TRUE)
      if (!is.null(input$scatter_brush)) {
        cur_brush <- input$scatter_brush
        cur_brush$mapping$x <- 'x'
        cur_brush$mapping$y <- 'y'
        res <- enclosed_brush(rv_scatter$data, cur_brush)
        rv_scatter$highlight_vars <- res
      }  else {
        res <- NULL
      }

      # TODO: total hack -- fix this correctly eventually
      if (is(res, 'data.frame')) {
        res <- dplyr::inner_join(
          data.table::as.data.table(sr),
          data.table::as.data.table(dplyr::select(res, target_id)),
          by = 'target_id')
        res <- dplyr::rename(res,
          mean = mean_obs,
          var = var_obs,
          tech_var = sigma_q_sq,
          final_sigma_sq = smooth_sigma_sq_pmax)
        res <- as_df(res)
      }

      res
    })

    ###
    output$samp_heat_plt <- renderPlot({
      plot_sample_heatmap(obj, use_filtered = input$samp_heat_filt)
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
        units = input$pca_units,
        point_size = input$pca_point_size
        )

    })

    ### MV plot
    output$mv_plt <- renderPlot({
      plot_mean_var(obj)
    })

    ### MA
    rv_ma <- reactiveValues(
      highlight_vars = NULL
      )

    output$which_beta_ctrl <- renderUI({
      poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
      selectInput('which_beta', 'beta: ', choices = poss_tests)
    })

    output$ma <- renderPlot({
      val <- input$which_beta
      if ( is.null(val) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
        val <- poss_tests[1]
      }
      plot_ma(obj, val,
        input$which_model,
        sig_level = input$max_fdr,
        point_alpha = input$ma_alpha
        )
    })

    # observe({
    #   if (!is.null(input$ma_trans)) {
    #     res <- data.frame(target_id = input$ma_trans, stringsAsFactors = FALSE)
    #     rv_ma$highlight_ma <- res
    #     rv_ma$highlight_vars <- res
    #     wb <- input$which_beta
    #     if ( is.null(wb) ) {
    #       poss_tests <- tests(models(obj)[[input$which_model]])
    #       wb <- poss_tests[1]
    #     }
    #     sres <- sleuth_results(obj, wb, input$which_model, rename_cols = TRUE, show_all = TRUE)
    #     output$ma_brush_out <- renderDataTable({
    #       dplyr::semi_join(sres, res, by = 'target_id' )
    #     })
    #   }
    # })

    output$vars <- renderPlot({
      wb <- input$which_beta
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
        wb <- poss_tests[1]
      }
      plot_vars(obj,
        which_beta = wb,
        which_model = input$which_model,
        point_alpha = input$ma_alpha,
        highlight = rv_ma$highlight_vars,
        sig_level = input$max_fdr
        )
    })

    #observe
    output$ma_brush_out <- renderDataTable({
      wb <- input$which_beta
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model]])
        wb <- poss_tests[1]
      }
      res <- sleuth_results(obj, wb, input$which_model, rename_cols = FALSE,
        show_all = FALSE)
      if (!is.null(input$ma_brush)) {
        res <- enclosed_brush(res, input$ma_brush)
        rv_ma$highlight_vars <- res
      }  else {
        res <- NULL
      }

      # TODO: total hack -- fix this correctly eventually
      if (is(res, 'data.frame')) {
        res <- dplyr::rename(res,
          mean = mean_obs,
          var = var_obs,
          tech_var = sigma_q_sq,
          final_sigma_sq = smooth_sigma_sq_pmax)
      }

      res
    })

    ### DE table
    output$which_beta_ctrl_de <- renderUI({
      poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model_de]])
      selectInput('which_beta_de', 'beta: ', choices = poss_tests)
    })

    output$de_dt <- renderDataTable({
      wb <- input$which_beta_de
      if ( is.null(wb) ) {
        poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model_de]])
        wb <- poss_tests[1]
      }
      sleuth_results(obj, wb, input$which_model_de)
    })

    bs_var_text <- eventReactive(input$bs_go, {
        input$bs_var_input
      })

    ### bootstrap var
    output$bs_var_plt <- renderPlot({
      plot_bootstrap(obj, bs_var_text(),
        units = input$bs_var_units,
        color_by = input$bs_var_color_by)
    })
  }

  shinyApp(ui = p_layout, server = server_fun)
}

#' @export
enclosed_brush <- function(df, brush) {
  df <- as_df(df)
  xvar <- brush$mapping$x
  yvar <- brush$mapping$y

  xbool <- brush$xmin <= df[[xvar]] & df[[xvar]] <= brush$xmax
  ybool <- brush$ymin <= df[[yvar]] & df[[yvar]] <= brush$ymax

  df[xbool & ybool,]
}
