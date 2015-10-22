library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Hello Shiny!"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("bins",
                  "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30)
    ),
    navbarMenu('maps',

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
              choices = 1:5
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
          )
        )
        )
      ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))