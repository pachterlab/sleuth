navbarMenu('maps',

	####
	tabPanel('PC variance',
	fluidRow(   
		column(12,
			p(h3('principal component variance'), "plot variances retained by each principal component")
		),
		offset = 1),
		fluidRow(
			column(3, #figure out dynamic PC parts
				selectInput('PC_relative', label = 'starting PC', choices = "", #this may pose a problem later
					selected = 1)
				),
			column(3
				numericInput('Number of Principal Components', label = 'pca_number', choices = 3:10,
					selected = 2)
				),
			),
			fluidRow(
				 column(3,
            checkboxInput('pc_filt', label = 'filter',
              value = TRUE),
            checkboxInput('bool', label = 'scale',
              value = TRUE)
            )
				),
			fluidRow(plotOutput('plt_pc_var')),
			),

		tabPanel('loadings',
			fluidRow(
				column(12,
					p(h3('loadings'), "plot loadings and gene contributions to principal components")
					),
				offset = 1),
			fluidRow(
				column(3, 
					selectInput('Number of principal components', label = 'pc_count', choices = 3:10,
						selected = 1)
					),
				column(3,
					selectInput('Which gene', label = 'gene', choices = '',
						selected = 1)
					),
			fluidRow(
			 column(3,
          checkboxInput('pcfilt', label = 'filter',
            value = TRUE),
          checkboxInput('absl', label = 'absolute value',
            value = TRUE),
          checkboxInput('bool', label = 'scale',
            value = TRUE)
          			)
				),
			 fluidRow(plotOutput('plt_pc_loadings')),
			)
		)
	)

### Plot pc loadings
output$plt_pc_loadings <- renderPlot({
	plot_loadings(obj,
		use_filtered = input$pcfilt,
		pc_count = input$pc_count,
		bool = input$bool,
		gene = input$gene,
		absolute = input$absl
		)
	})


###plot pc variance
output$plt_pc_var <- renderPlot({
	plot_pc_variance(obj,
		use_filtered = input$pc_filt,
		pca_number = input$pca_number,
		bool = input$bool,
		PC_relative = input$PC_relative
		)
	})
