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

#' Interactive sleuth visualization with Shiny
#'
#' Interactive sleuth visualization with Shiny. To exit, type \code{ESC} in R.
#'
#' @param obj a \code{sleuth} object already processed and has run
#' \code{\link{sleuth_fit}} and \code{\link{sleuth_wt}} or \code{\link{sleuth_lrt}}
#' @param settings see the function \code{\link{sleuth_live_settings}} for options
#' @param options additional options which are sent to shiny
#' @param ... additional parameters sent to plotting functions
#' @return a \code{\link{shinyApp}} result
#' @export
#' @seealso \code{\link{sleuth_fit}}, \code{\link{sleuth_live_settings}}, \code{\link{sleuth_deploy}}
sleuth_live <- function(obj, settings = sleuth_live_settings(),
  options = list(port = 42427), ...) {
  stopifnot( is(obj, 'sleuth') )
  if ( !require('shiny') ) {
    stop("'sleuth_live()' requires 'shiny'. Please install it using
      install.packages('shiny')")
  }

  if (obj$gene_mode) {
    counts_unit <- "scaled_reads_per_base"
    bsg_choices <- setdiff(names(obj$target_mapping), c('target_id'))
    index <- which(bsg_choices == obj$gene_column)
    bsg_choices <- c(bsg_choices[index], bsg_choices[-index])
    gene_mode_choice <- "true"
  } else {
    counts_unit <- "est_counts"
    gene_mode_choice <- "false"
    bsg_choices <- c()
  }

  # set up for the different types of tests
  poss_covars <- dplyr::setdiff(
    colnames(obj$sample_to_covariates),
    'sample')
  samp_names <- obj$sample_to_covariates[['sample']]
  poss_models <- names(models(obj, verbose = FALSE))

  poss_wt <- list_tests(obj, 'wt')
  poss_lrt <- list_tests(obj, 'lrt')
  valid_test_types <- if (!is.null(poss_wt)) {
    c('Wald' = 'wt')
  } else {
    c()
  }
  if (!is.null(poss_lrt)) {
    valid_test_types <- c(valid_test_types, c('likelihood ratio' = 'lrt'))
  }

  # if (length(valid_test_types) == 0) {
  #   stop("We found no valid tests. Please add some tests and rerun sleuth_live()")
  # }

  p_layout <- navbarPage(
    a('sleuth', href = 'http://pachterlab.github.io/sleuth', target = '_blank',
      style = 'color: black;'),
    windowTitle = 'sleuth',
    tabPanel('overview',
      fluidRow(
       div(h3('sleuth live'), align = 'center')
       ),
     fluidRow(
       column(10, offset = 1,
          p('This Shiny app is designed for exploratory data analysis of
          kallisto-sleuth processed RNA-Seq data. There are four menu tabs
          that can be used to choose plots and tables to view.'),

          p(strong('sleuth live features:')),
             tags$ul(
               tags$li(strong('v0.29.0'),
                       ':',
                       'Integration of gene mode, and enhancements to heat maps',
                       'by Warren McGee.'
               ),
                    tags$li(strong('v0.28.0'),
                      ':',
                      'Download buttons for plots and tables by Alex Tseng.',
                      'PCA variance explained and loadings by Daniel Li.',
                      'Fragment length distribution plot, bias table,',
                        'and integration of likelihood ratio test by Harold Pimentel.'
                      ),
                     tags$li(strong('v0.27.3'),
                       ':',
                        'Gene table, gene viewer, transcript heatmap,',
                          'and volcano plot by Pascal Sturmfels.'
                      ),
                     tags$li(strong('v0.27.2'),
                       ':',
                       'Design matrix, kallisto table, transcript view,',
                       'and QQplot by Harold Pimentel.'
                     ),
                     tags$li(strong('v0.27.1'),
                      ':',
                      'Densities, MA plot, mean-variance plot, PCA,',
                      'processed data, sample heatmap, scatter plots,',
                      'and test table by Harold Pimentel.')
                     )
          ))),

    navbarMenu('analyses',

      tabPanel('gene view',
        conditionalPanel(
          condition = 'input.settings_gene_mode == "true"',
          fluidRow(
            column(12,
                   p(
                     h3('gene view'),
                     'Boxplots of gene abundances showing technical variation in each sample.')
            ),
            offset = 1),
          fluidRow(
            column(3,
                   textInput('bsg_var_input',
                      label = HTML('gene: ',
                                   '<button onclick="bsg_var_input()">?</button>',
                                   '<script> function bsg_var_input() {',
                                   'alert("Enter the target_id of a gene here to view a boxplot of its technical variation.',
                                   'You can find target_ids in the test table under Analyses.',
                                   'If you select a different column from "genes from",',
                                   'enter a gene identifier that matches that column.',
                                   'The most significant gene by q-value is already entered.");',
                                   '} </script>'),
                      value = (obj$tests[[1]][[1]])[
                        order(obj$tests[[1]][[1]]$qval),]$target_id[1])
            ),
            column(3,
                   selectInput('bsg_var_color_by',
                      label = HTML('color by: ',
                                   '<button onclick="bsg_var_color_by()">?</button>',
                                   '<script> function bsg_var_color_by() {',
                                   'alert("Color the box plots by different column names from your design file.',
                                   'Doing so can help explain which traits of a given gene best explain the',
                                   'difference across samples.");',
                                   '} </script>'),
                      choices = c(NULL, poss_covars), selected = NULL)
            ),
            column(3,
                   selectInput('bsg_var_units',
                               label = HTML('units: ',
                                            '<button onclick="bsg_var_units()">?</button>',
                                            '<script> function bsg_var_units() {',
                                            'alert("tpm: Transcripts Per Million\\n',
                                            'scaled_reads_per_base: The average number of reads mapping to each base across the whole gene\\n',
                                            '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                                            '} </script>'),
                               choices  = (if (length(obj$bs_quants) == 0) { c('N/A') }
                                           else { names(obj$bs_quants[[1]]) } ),
                               selected = (if (length(obj$bs_quants) == 0) { 'N/A' }
                                           else { names(obj$bs_quants[[1]])[1] })
                   )),
            column(3,
                   selectInput('bsg_gene_colname',
                               label = HTML('genes from: ',
                                            '<button onclick="bsg_gene_colname()">?</button>',
                                            '<script> function bsg_gene_colname() {',
                                            'alert("Choose a column name from your target mapping data frame.',
                                            'This allows you to use alternative gene identifiers (e.g. HGNC symbol)',
                                            'to select a gene instead of the target ID.',
                                            'The column you used to aggregate at the gene level is the default");',
                                            '} </script>'),
                               choices = bsg_choices)
            )
          ),
          fluidRow(HTML('&nbsp;&nbsp;&nbsp;'), actionButton('bsg_go', 'view')),
          fluidRow(plot <- (if (length(obj$bs_quants) == 0) {
            HTML('&nbsp&nbsp&nbsp&nbsp You need to run sleuth with at least ',
                 'one of extra_bootstrap_summary or ',
                 'read_bootstrap_tpm to use this feature.<br>') }
            else { plotOutput('bsg_var_plt') }
          )),
          fluidRow(
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_bsg_var_plt", "Download Plot"))
          )
        ),
        conditionalPanel(
          condition = 'input.settings_gene_mode == "false"',
          fluidRow(
            column(12,
                   p(h3('gene view'),
                     'Boxplots of abundances of transcripts mapping to a given gene,',
                     'and their technical variation.',
                     'This step can take a while, especially with many plots.')
            ),
            offset = 1),
          fluidRow(column(3,
                          textInput('gv_var_input',
                                    label = HTML('gene: ',
                                                 '<button onclick="gv_var_input()">?</button>',
                                                 '<script> function gv_var_input() {',
                                                 'alert("Enter the gene identifier (determined by what is selected for \'genes from\')',
                                                 'here to view a boxplot of the technical variation of its top transcripts ranked by q-value.',
                                                 'The number is determined by \'# of plots\'.',
                                                 'You can find gene identifiers in the test table under Analyses.");',
                                                 '} </script>'), value = '')
          ),
          column(3,
                 selectInput('gv_var_color_by',
                             label = HTML('color by: ',
                                          '<button onclick="gv_var_color_by()">?</button>',
                                          '<script> function gv_var_color_by() {',
                                          'alert("Color the box plots by different column names from your design file.',
                                          'Doing so can help explain which traits of a given transcript best explain the',
                                          'difference in counts/TPMs of that transcript across samples.");',
                                          '} </script>'),
                             choices = c(NULL, poss_covars),
                             selected = NULL)),
          column(3,
                 selectInput('gv_var_units',
                             label = HTML('units: ',
                                          '<button onclick="gv_var_units()">?</button>',
                                          '<script> function gv_var_units() {',
                                          'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
                                          '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                                          '} </script>'),
                             choices = c(counts_unit, 'tpm'),
                             selected = counts_unit)),
          column(3,
                 uiOutput('gv_gene_column')
          )
          ),
          fluidRow(
            column(3, actionButton('gv_go', 'view')),
            column(3, numericInput('gv_maxplots', label = '# of plots (max 15): ', value = 3,
                                   min = 1, max = 15, step = 1))
          ),
          fluidRow(uiOutput('no_genes_message')),
          fluidRow(plot <- (if (length(obj$bs_quants) == 0) {
            HTML(paste('&nbsp&nbsp&nbsp&nbsp You need to run sleuth with at ',
                       'least one of extra_bootstrap_summary or',
                       'read_bootstrap_tpm to use this feature.<br>'))
          } else { uiOutput('gv_var_plts') }
          ))
        )
      ),

      ####
      tabPanel('heat map',
        fluidRow(
          column(12,
              p(h3('heat map'),
              HTML('Plot of select abundances in a clustered heat map using the',
              'R <a target="_blank" href="https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html">',
              'hclust function</a>. The heat map iteratively groups together transcripts with',
              'similar expression patterns, highlighting similarities among the transcripts themselves.',
              'For a more in-depth discussion of hierarchical clustering, see',
              '<a target="_blank" href="http://link.springer.com/article/10.1007/BF02289588">10.1007/BF02289588</a>.',
              '</br>Enter space-separated values (The ten transcripts with the lowest q-values',
              'are already entered. Click view to see them!).'))
              ),
          offset = 1
        ),
        fluidRow(
          column(4,
            textInput('hm_transcripts', label = HTML('enter target ids: ',
            '<button onclick="hm_transcripts()">?</button>',
            '<script> function hm_transcripts() {',
            'alert("Enter a space-separated list of transcript names here to view a hierarchical clustering of those transcripts.");',
            '} </script>'),
            value = paste((obj$tests[[1]][[1]])[
              order(obj$tests[[1]][[1]]$qval),]$target_id[1:10], collapse = " "))
              ),
          column(8, style = "margin-top:15px;",
            checkboxGroupInput('hm_covars',
                               label = HTML('covariates:',
                                            '<button onclick="hm_covars()">?</button>',
                                            '<script> function hm_covars() {',
                                            'alert("Add annotation information from different column names from your design file.',
                                            'Doing so can help explain which traits of a given sample best explain the',
                                            'clustering across samples.");',
                                            '} </script>'),
                               choices = as.list(poss_covars), inline = TRUE)
                 )),
        fluidRow(
          column(3,
            selectInput('hm_units', label = HTML('units: ',
            '<button onclick="hm_units()">?</button>',
            '<script> function hm_units() {',
            'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
            'scaled_reads_per_base (gene mode only): The average number of reads mapping to each base across the whole gene\\n',
            '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
            '} </script>'),
            choices = c(counts_unit, 'tpm'), selected = 'tpm')
              ),
          column(2,
            textInput('hm_trans', label = HTML('transform: ',
            '<button onclick="hm_trans()">?</button>',
            '<script> function hm_trans() {',
            'alert("A transformation to be applied to the raw data before clustering.");',
            '} </script>'), value = 'log')
              ),
          column(2,
            numericInput('hm_offset', label = HTML('offset: ',
            '<button onclick="hm_offset()">?</button>',
            '<script> function hm_offset() {',
            'alert("A constant amount to be added to the raw data before the transformation and clustering.");',
            '} </script>'), value = 1)),
          column(2, style = "margin-top: 10px;",
            checkboxInput('hm_cluster', label = HTML('cluster transcripts',
            '<button onclick="hm_cluster()">?</button>',
            '<script> function hm_trans() {',
            'alert("Check this box to use hierarchical clustering on the selected transcripts.\\n',
            'Note that hierarchical clustering of the samples is already turned on by default.");',
            '} </script>'),
            value = TRUE)),
          column(1,
            actionButton('hm_go', 'view')
          )
        ),
        tags$style(type = 'text/css', "#hm_go {margin-top: 25px}"),
        fluidRow(plotOutput('hm_plot')),
        fluidRow(uiOutput("download_hm_plt_button"))
    ),


      ####
      tabPanel('MA plot',
      fluidRow(
        column(12,
          p(h3('MA plot'),
            'Plot of abundance versus fixed effect (e.g. fold change).',
            'Select a set of transcripts to explore their variance across samples.')
          ),
          offset = 1),
          conditionalPanel(condition = 'input.settings_test_type == "wt"',
            fluidRow(
              column(2,
                numericInput('max_fdr', label = 'max Fdr:', value = 0.10,
                  min = 0, max = 1, step = 0.01)),
              column(4,
                selectInput('which_model_ma', label = 'fit: ',
                  choices = poss_models,
                  selected = poss_models[1])
                ),
              column(4,
                uiOutput('which_beta_ctrl_ma')
                ),
              column(2,
                numericInput('ma_alpha', label = 'opacity:', value = 0.2,
                  min = 0, max = 1, step = 0.01))
              ),

            fluidRow(plotOutput('ma', brush = 'ma_brush')),
            fluidRow(
              div(align = "right", style = "margin-right:15px",
                  downloadButton("download_ma_plt", "Download Plot"))),
            fluidRow(plotOutput('vars')),
            fluidRow(
              div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                  downloadButton("download_ma_var_plt", "Download Plot"))),
            fluidRow(dataTableOutput('ma_brush_out')),
            fluidRow(uiOutput("download_ma_table_button"))
            ),

          conditionalPanel(condition = 'input.settings_test_type == "lrt"',
            strong(paste("Only supported for 'setting' Wald tests.",
                         "Go to the settings tab to change the 'test type' setting."))
          )
        ),



      ####
      tabPanel('test table',
        conditionalPanel(condition = 'input.settings_test_type == "wt"',
          fluidRow(column(12, p(h3('test table'),
                     HTML('Table of transcript names, gene names (if supplied),',
                          'sleuth parameter estimates, tests, and summary statistics.',
                          '<button onclick="tableType()">What do the column names mean?</button>',
                          '<script> function tableType() {',
                          'alert("Transcript Table Columns:\\n',
                          'target_id: Transcript name (dependent on transcriptome used in kallisto)\\n',
                          'pval: p-value of the transcript for the selected test\\n',
                          'qval: False discovery rate adjusted p-value, using Benjamini-Hochberg\\n',
                          'b: the "beta" value (effect size). Technically a biased estimator of the fold change\\n',
                          'se_b: the standard error of the beta\\n',
                          'mean_obs: Mean of natural log counts of observations\\n',
                          'var_obs: Variance of observation\\n',
                          'tech_var: Technical variance of observation from the bootstraps\\n',
                          'sigma_sq: Raw estimator of the variance once the technical variance has been removed\\n',
                          'sigma_sq_pmax: max(sigma_sq, 0)\\n',
                          'smooth_sigma_sq: Smooth regression fit for the shrinkage estimation\\n',
                          'final_simga_sq: max(sigma_sq, smooth_sigma_sq) – used for covariance estimation of beta\\n',
                          '");',
                          '} </script>'))
          ), offset = 1),
          fluidRow(
            column(3,
              selectInput('which_model_de', label = 'fit: ',
                choices = poss_models,
                selected = poss_models[1])
              ),
            column(3,
              uiOutput('which_beta_ctrl_de')
              ),
            column(3,
              uiOutput('table_type')
              ),
            column(3,
              uiOutput('group_by')
              )
            ),
            dataTableOutput('de_dt'),
            fluidRow(
              div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                  downloadButton("download_test_table",  "Download Table")))
          ),
        conditionalPanel(condition = 'input.settings_test_type == "lrt"',
            fluidRow(column(12, p(h3('test table'),
                       HTML('Table of transcript names, gene names (if supplied),',
                            'sleuth parameter estimates, tests, and summary statistics.',
                            '<button onclick="tableType()">What do the column names mean?</button>',
                            '<script> function tableType() {',
                            'alert("Transcript Table Columns:\\n',
                            'target_id: Transcript name (dependent on transcriptome used in kallisto)\\n',
                            'pval: p-value of the transcript for the selected test\\n',
                            'qval: False discovery rate adjusted p-value, using Benjamini-Hochberg\\n',
                            'test_stat: Chi-squared test statistic (likelihood ratio test)\\n',
                            'rss: residual sum of squares under the null model\\n',
                            'degrees_free: the degrees of freedom for the test_that\\n',
                            'mean_obs: Mean of natural log counts of observations\\n',
                            'var_obs: Variance of observation\\n',
                            'tech_var: Technical variance of observation from the bootstraps\\n',
                            'sigma_sq: Raw estimator of the variance once the technical variance has been removed\\n',
                            'sigma_sq_pmax: max(sigma_sq, 0)\\n',
                            'smooth_sigma_sq: Smooth regression fit for the shrinkage estimation\\n',
                            'final_simga_sq: max(sigma_sq, smooth_sigma_sq) – used for covariance estimation of beta\\n',
                            '");',
                            '} </script>'))
                         ), offset = 1),
          fluidRow(
            column(3,
              uiOutput('test_control_de')
              ),
            column(3,
              uiOutput('table_type_lrt')
              ),
            column(3,
              uiOutput('group_by_lrt')
              )
            ),
          dataTableOutput('lrt_de_dt')
          )
        ),

        ####
        tabPanel('transcript view',
          uiOutput('transcript_view')
        ),

        ####
        tabPanel('volcano plot',
          fluidRow(
            column(12,
              p(
                h3('volcano plot'),
                'Plot of beta value (regression) versus log of significance.',
                'Select a set of transcripts to explore their variance across samples.'
                )
              ),
            offset = 1),
          conditionalPanel(condition = 'input.settings_test_type == "wt"',
            fluidRow(
              column(2,
                numericInput('max_fdr_vol', label = 'max Fdr:', value = 0.10,
                    min = 0, max = 1, step = 0.01)),
              column(4,
                selectInput('which_model_vol', label = 'fit: ',
                  choices = poss_models,
                  selected = poss_models[1])
              ),
              column(4,
                uiOutput('which_beta_ctrl_vol')
                ),
              column(2,
                numericInput('vol_alpha', label = 'opacity:', value = 0.2,
                  min = 0, max = 1, step = 0.01))
                ),
            fluidRow(plotOutput('vol', brush = 'vol_brush')),
            fluidRow(
              div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_volcano_plt", "Download Plot"))),
            fluidRow(dataTableOutput('vol_brush_out')),
            fluidRow(uiOutput("download_volcano_table_button"))
          ),

          conditionalPanel(condition = 'input.settings_test_type == "lrt"',
            strong('Only supported for "setting" Wald test.')
          )
        )
      ),

    navbarMenu('maps',

    ####
      tabPanel('PCA',
      fluidRow(
        column(12,
          p(
            h3('principal component analysis'),
            'PCA projections of sample abundances onto any pair of components.',
              ' PCA is computed on the transcript expression.',
              ' Each sample is a vector with dimension equal to the number of transcripts.',
              ' See',
            a(href="https://liorpachter.wordpress.com/2014/05/26/what-is-principal-component-analysis/",
              target="_blank", 'this blog post'),
            'for an overview of PCA.')
          ),
          offset = 1),
        fluidRow(
          column(2,
            selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
              selected = 1)
            ),
          column(2,
            selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
              selected = 2)
            ),
          column(3,
            selectInput('color_by',
              label = HTML('color by: ',
                           '<button onclick="color_by()">?</button>',
                           '<script> function color_by() {',
                           'alert("Color the PCA plot by different column names from your design file.',
                           'Doing so can help explain which traits of a given sample best explain the',
                           'clustering of samples on the PCA plot.");',
                           '} </script>'),
              choices = c(NULL, poss_covars), selected = NULL)
            ),
          column(2,
            numericInput('pca_point_size', label = 'size: ', value = 3)),
          column(3,
            selectInput('pca_units',
              label = HTML('units: ',
                           '<button onclick="pca_units()">?</button>',
                           '<script> function pca_units() {',
                           'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
                           'scaled_reads_per_base (gene mode only): The average number of reads mapping to each base across the whole gene\\n',
                           '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                           '} </script>'),
              choices = c(counts_unit, 'tpm'),
              selected = counts_unit))
        ),
        fluidRow(
          column(2,
            checkboxInput('pca_filt',
              label = HTML('filter',
                           '<button onclick="pca_filt()">?</button>',
                           '<script> function pca_filt() {',
                           'alert("Check this box to use only those transripts that passed the filter for the PCA analysis.");',
                           '} </script>'),
              value = TRUE)
              ),
          column(2,
            checkboxInput('text_labels', label = 'text labels',
              value = TRUE)
            )
          ),
        fluidRow(plotOutput('pca_plt')),
        fluidRow(
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_pca_plt", "Download Plot"))),
        fluidRow(
        column(12,
          p(
            h3('loadings'),
            'observe contributions of samples or transcripts to the principal component')
          ),
          offset = 1),
        fluidRow(
          column(3,
            textInput('sample', label = 'transcript: ', value = '',
              )
            ),
          column(2,
            selectInput('pc_input', label = 'principal component: ', choices = 1:5,
              selected = 1)
            ),
          column(3,
            selectInput('pc_count', label = 'number of PCs or transcripts: ', choices = 1:10,
              selected = 5)),
          column(2,
            checkboxInput('pca_loading_abs', label = 'absolute value',
              value = TRUE)
            ),
          column(2,
            checkboxInput('scale', label = 'scale',
              value = FALSE)
            )
          ),

        fluidRow(
          column(12,
            p(h3('variance explained'))
            ),
            offset = 1),
        fluidRow(plotOutput('plt_pc_var')),
        fluidRow(
          column(12,
            p(h3('loadings'))
            ),
            offset = 1),
        fluidRow(plotOutput('plt_pc_loadings'))
        ),

      ###
      tabPanel('sample heatmap',
      fluidRow(
        column(12,
          p(h3('sample heatmap'),
            a(href='https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence',
              target='_blank', 'Jensen-Shannon divergence'),
            'between pairs of samples computed at the transcript level.',
            ' The larger this value, the more dissimilar the samples are.')
          ),
          offset = 1),
        fluidRow(column(1, style = "margin-top: 15px;",
                   checkboxInput('samp_heat_filt',
                                 label = HTML('filter',
                                              '<button onclick="samp_heat_filt()">?</button>',
                                              '<script> function samp_heat_filt() {',
                                              'alert("Check this box to use only those transripts that passed the filter before clustering.");',
                                              '} </script>'),
                                 value = TRUE)),
                 column(2, style = "margin-top: 15px;",
                   checkboxInput('samp_heat_cluster',
                                 label = HTML('cluster samples',
                                              '<button onclick="samp_heat_cluster()">?</button>',
                                              '<script> function samp_heat_cluster() {',
                                              'alert("Check this box to use hierarchical clustering on the samples.");',
                                              '} </script>'), value = TRUE)
                   ),
                 column(9,
                   checkboxGroupInput('samp_heat_covars',
                                      label = HTML('covariates:',
                                                   '<button onclick="samp_heat_covars()">?</button>',
                                                   '<script> function samp_heat_covars() {',
                                                   'alert("Add annotation information from different column names from your design file.',
                                                   'Doing so can help explain which traits of a given sample best explain the',
                                                   'clustering across samples.");',
                                                   '} </script>'),
                                      choices = as.list(poss_covars), inline = TRUE)
                   )),
        fluidRow(plotOutput('samp_heat_plt')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
            downloadButton("download_samp_heat_plt", "Download Map")))
        )
      ),

    navbarMenu('summaries',

      ####
      tabPanel('densities',
      fluidRow(
        column(12,
          p(
            h3('distribution of abundances'),
            'Distributions of abundances of individual samples or groupings by covariates.'
            )
          ),
          offset = 1),
        fluidRow(
          column(4,
            selectInput('cond_dens_grp', 'grouping: ',
              choices = poss_covars,
              selected = poss_covars[1])
            ),
          column(2,
            selectInput('cond_dens_units',
              label = HTML('units: ',
                           '<button onclick="cond_dens_units()">?</button>',
                           '<script> function cond_dens_units() {',
                           'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
                           'scaled_reads_per_base (gene mode only): The average number of reads mapping to each base across the whole gene\\n',
                           '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                           '} </script>'),
              choices = c('tpm', counts_unit),
              selected = 'tpm')),
          column(2,
            checkboxInput('cond_dens_filt',
              label = HTML('filter',
                           '<button onclick="cond_dens_filt()">?</button>',
                           '<script> function cond_dens_filt() {',
                           'alert("Check this box to use only those transripts',
                           'that passed the filter before plotting.");',
                           '} </script>'),
              value = TRUE)),
          column(2,
            textInput('cond_dens_trans',
              label = HTML('transform: ',
                           '<button onclick="cond_dens_trans()">?</button>',
                           '<script> function cond_dens_trans() {',
                           'alert("A transformation to be applied to the raw data before plotting.");',
                           '} </script>'),
              value = 'log')),
          column(2,
            numericInput('cond_dens_offset',
                         label = HTML('offset: ',
                                      '<button onclick="cond_dens_offset()">?</button>',
                                      '<script> function cond_dens_offset() {',
                                      'alert("A constant amount to be added to the raw data before the transformation and clustering.");',
                                      '} </script>'),
                         value = 1))
          ),
        fluidRow(plotOutput('condition_density')),
        fluidRow(
          div(align = "right", style = "margin-right:15px",
            downloadButton("download_cond_dens_plt", "Download Plot"))),
        fluidRow(
          column(4,
            selectInput('samp_dens', 'sample: ',
              choices = samp_names,
              selected = samp_names[1]))),
        fluidRow(plotOutput('sample_density')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom: 10px",
            downloadButton("download_samp_dens_plt", "Download Plot")))
        ),


      ###
      tabPanel('design matrix',

      fluidRow(
        column(12,
          p(h3('design matrix'), "View the design matrix used to fit each model.")
          ),
          offset = 1),

        fluidRow(
          column(4,
            selectInput('which_model_design', label = 'fit: ',
              choices = poss_models,
              selected = poss_models[1])
            )
          ),

        fluidRow(
          verbatimTextOutput('design_matrix')
          #tableOutput('design_matrix')
          )
        ),

      ####
      tabPanel('fragment length distribution',
        fluidRow(
          column(12,
            p(
              h3('fragment length distribution plot'),
              'Plot fragment length distribution used by kallisto in a particular sample.',
              ' In paired-end data, kallisto learns this fragment length distribution by looking at the beginning coordinate and at the end coordinate of unique pseudoalignments.',
              ' In single-end data, the fragment length distribution is provided by the user.')
            )
        ),
        fluidRow(
          column(4,
            selectInput('fld_sample', label = 'sample: ',
              choices = samp_names,
              selected = samp_names[1]
            )
          )
        ),
        fluidRow(
          plotOutput('fld_plt')
          )
      ),

      ####
      tabPanel('processed data',
      fluidRow(
        column(12,
          p(
            h3('processed data'),
            'Names of samples, number of mapped reads, number of boostraps performed by kallisto,',
            'and sample to covariate mappings.'
            )
          ),
          offset = 1),
        fluidRow(
          column(12,
            p(strong('kallisto version(s): '), obj$kal_versions)),
          offset = 1
          ),
        fluidRow(dataTableOutput('summary_dt')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom: 10px",
            downloadButton("download_summary_table", "Download Table")))
        ),

      ###
      tabPanel('kallisto table',

        fluidRow(
        column(12, p(h3('kallisto abundance table'),
          "All of the abundance estimates pulled in from kallisto results into the
           sleuth object.",
           ' The covariates button will include covariates from the `sample_to_covariates` table defined in `sleuth_prep`.',
           ' `eff_len` and `len` are the effective length and true length of the transcript (`target_id`).',
           ' An overview of these units can be found ',
           a(href = 'https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/',
            target = '_blank', 'here'), '.'
           ))
          ),

        fluidRow(
          column(3,
            checkboxInput('norm_tbl', label = 'normalized ',
              value = TRUE)),
          column(3,
            checkboxInput('filt_tbl', label = 'filter ',
              value = TRUE)),
          column(3,
            checkboxInput('covar_tbl', label = 'covariates ',
              value = FALSE))
          ),

        fluidRow(dataTableOutput('kallisto_table')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
            downloadButton("download_kallisto_table", "Download Table")))
        )
      ),

      navbarMenu('diagnostics',

      ####
      tabPanel('bias weights',
        fluidRow(
          column(12,
            p(h3('bias weights'),
            "View the bias parameters modeled by kallisto.")
            )
          ),

        fluidRow(
          column(4,
            selectInput('bias_sample', label = 'sample: ',
              choices = samp_names,
              selected = samp_names[1]
            )
          )
        ),

        fluidRow(
          dataTableOutput('bias_weights_table')
          )
      ),

      ####
      tabPanel('mean-variance plot',
      fluidRow(
        column(12,
          p(
            h3('mean-variance plot'),
            'Plot of abundance versus square root of standard deviation which is used for shrinkage estimation.',
            'The blue dots are in the interquartile range and the red curve is the fit used by sleuth.',
            ' Any points at y = 0 have inferential variance greater than observed total raw variance.'
            )
          ),
          offset = 1),
        fluidRow(plotOutput('mv_plt')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
            downloadButton("download_mv_plt", "Download Plot")))

        ),

      tabPanel('scatter plots',
        ####
        fluidRow(
          column(12,
            p(
              h3('scatter plot '),
              "Display scatter plot for any two samples and then select",
              "a set of transcripts to explore their variance across samples.")
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
            textInput('trans',
              label = HTML('transform: ',
                           '<button onclick="trans()">?</button>',
                           '<script> function trans() {',
                           'alert("A transformation to be applied to the raw data before plotting.");',
                           '} </script>'),
              value = 'log')),
          column(2,
            numericInput('scatter_offset',
                         label = HTML('offset: ',
                                      '<button onclick="scatter_offset()">?</button>',
                                      '<script> function scatter_offset() {',
                                      'alert("A constant amount to be added to the raw data before the transformation and clustering.");',
                                      '} </script>'),
                         value = 1))
          ),
        fluidRow(
          column(2,
            selectInput('scatter_units',
              label = HTML('units: ',
                           '<button onclick="scatter_units()">?</button>',
                           '<script> function scatter_units() {',
                           'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
                           'scaled_reads_per_base (gene mode only): The average number of reads mapping to each base across the whole gene\\n',
                           '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                           '} </script>'),
              choices = c(counts_unit, 'tpm'),
              selected = counts_unit)),
          column(2,
            checkboxInput('scatter_filt',
              label = HTML('filter',
                           '<button onclick="scatter_filt()">?</button>',
                           '<script> function scatter_filt() {',
                           'alert("Check this box to use only those transripts',
                           'that passed the filter before plotting.");',
                           '} </script>'),
              value = TRUE)),
          column(2,
            numericInput('scatter_alpha', label = 'opacity:', value = 0.2,
              min = 0, max = 1, step = 0.01))
          ),
        fluidRow(plotOutput('scatter', brush = 'scatter_brush')),
        fluidRow(
          div(align = "right", style = "margin-right:15px",
            downloadButton("download_scatter_plt", "Download Plot"))),
        fluidRow(plotOutput('scatter_vars')),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom: 10px",
            downloadButton("download_scatter_var_plt", "Download Plot"))),
        fluidRow(dataTableOutput('scatter_brush_table')),
        fluidRow(uiOutput("download_scatter_table_button"))
        ),

      tabPanel('Q-Q plot',
        ####
        fluidRow(
          column(12,
            p(
              h3('Q-Q plot'),
              "Select the test and view the appropriate quantile-quantile plot.",
              'Points that are color red are considered "significant" at the selected fdr-level.')
            ),
          offset = 1
          ),
        fluidRow(
          column(2,
            numericInput('max_fdr_qq', label = 'max Fdr:', value = 0.10,
              min = 0, max = 1, step = 0.01)),

          conditionalPanel(condition = 'input.settings_test_type == "wt"',
              column(4,
              selectInput('which_model_qq', label = 'fit: ',
                choices = poss_models,
                selected = poss_models[1])
              ),
              column(4,
                uiOutput('which_beta_ctrl_qq')
                )
            ),

          conditionalPanel(condition = 'input.settings_test_type == "lrt"',
            column(4,
              selectInput('test_qq', label = 'test: ',
                choices = poss_lrt, selected = poss_lrt[1]))
          )

          ),
        fluidRow(
          plotOutput('qqplot')
          ),
        fluidRow(
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
            downloadButton("download_qq_plt", "Download Plot")))
        )

      ),
    ####
      tabPanel('settings',
        fluidRow(
          column(2,
            selectInput('settings_test_type',
              label = 'test type:',
              choices = valid_test_types,
              selected = valid_test_types[1])
          ),
          column(2,
            selectInput('settings_gene_mode',
              label = 'gene mode:',
              choices = gene_mode_choice,
              selected = gene_mode_choice[1])
          )
        )
      )

    ) # navbarPage

  server_fun <- function(input, output) {
    # Reactive master object that stores all plots and tables for downloading later
    saved_plots_and_tables <- reactiveValues(
      pca_plt = NULL,
      samp_heat_plt = NULL,
      ma_plt = NULL,
      ma_var_plt = NULL,
      ma_table = NULL,
      test_table = NULL,
      volcano_plt = NULL,
      volcano_table = NULL,
      mv_plt = NULL,
      scatter_plt = NULL,
      scatter_var_plt = NULL,
      scatter_table = NULL,
      qq_plt = NULL,
      cond_dens_plt = NULL,
      samp_dens_plt = NULL,
      sample_table = NULL,
      kallisto_table = NULL,
      hm_plt = NULL,
      bs_var_plt = NULL
      )
    user_settings <- reactiveValues(save_width = 45, save_height = 11)
    # TODO: Once user settings are available, read these values from input

    output$which_beta_ctrl_qq <- renderUI({
      current_ui <- NULL
      poss_tests <- list_tests(obj, input$settings_test_type)
      if (settings$test_type == 'wt') {
        poss_tests <- poss_tests[[input$which_model_qq]]
        current_ui <- selectInput('which_test_qq', 'beta: ',
          choices = poss_tests, selected = poss_tests[1])
      } else {
        # TODO: I believe this code is defunct due to the conditionalPanel()
        current_ui <- selectInput('which_test_qq', 'test: ',
          choices = poss_tests, selected = poss_tests[1])
      }

      current_ui
    })

    output$qqplot <- renderPlot({
      poss_tests <- list_tests(obj, input$settings_test_type)
      current_test <- NULL
      if (input$settings_test_type == 'wt') {
        poss_tests <- poss_tests[[input$which_model_qq]]
      }
      current_test <- poss_tests[1]

      qq_plt <- plot_qq(obj, current_test,
        test_type = input$settings_test_type,
        which_model = input$which_model_qq,
        sig_level = input$max_fdr_qq)
      saved_plots_and_tables$qq_plt <- qq_plt
      qq_plt
    })

    output$download_qq_plt <- downloadHandler(
        filename = function() {
          "qq_plot.pdf"
          },
        content = function(file) {
            ggsave(file, saved_plots_and_tables$qq_plt,
              width = user_settings$save_width,
              height = user_settings$save_height,
              units = "cm")
    })

    output$summary_dt <- renderDataTable({
        saved_plots_and_tables$summary_table <- summary(obj)
        saved_plots_and_tables$summary_table
    })

    output$download_summary_table <- downloadHandler(
      filename = function() {
        "processed_data_table.csv"
      },
      content = function(file) {
         write.csv(saved_plots_and_tables$summary_table, file)
    })

    output$condition_density <- renderPlot({
      saved_plots_and_tables$cond_dens_plt <- plot_group_density(obj,
        grouping = input$cond_dens_grp,
        units = input$cond_dens_units,
        use_filtered = input$cond_dens_filt,
        trans = input$cond_dens_trans,
        offset = input$cond_dens_offset)
      saved_plots_and_tables$cond_dens_plt
    })

    output$download_cond_dens_plt <- downloadHandler(
        filename = function() {
          "condition_density_plot.pdf"
          },
        content = function(file) {
            ggsave(file,
              saved_plots_and_tables$cond_dens_plt,
              width = user_settings$save_width,
              height = user_settings$save_height,
              units = "cm")
    })

    output$sample_density <- renderPlot({
      saved_plots_and_tables$samp_dens_plt <- plot_sample_density(obj,
        which_sample = input$samp_dens,
        units = input$cond_dens_units,
        use_filtered = input$cond_dens_filt,
        trans = input$cond_dens_trans,
        offset = input$cond_dens_offset
        )
      saved_plots_and_tables$samp_dens_plt
    })

    output$download_samp_dens_plt <- downloadHandler(
        filename = function() {
          "sample_density_plot.pdf"
          },
        content = function(file) {
            ggsave(file,
              saved_plots_and_tables$samp_dens_plt,
              width = user_settings$save_width,
              height = user_settings$save_height,
              units = "cm")
    })

    ###
    output$kallisto_table <- renderDataTable({
      saved_plots_and_tables$kallisto_table <- kallisto_table(obj,
        use_filtered = input$filt_tbl,
        normalized = input$norm_tbl,
        include_covariates = input$covar_tbl
        )
      saved_plots_and_tables$kallisto_table
    })

    output$download_kallisto_table <- downloadHandler(
      filename = function() {
        "kallisto_table.csv"
        },
      content = function(file) {
         write.csv(saved_plots_and_tables$kallisto_table, file)
    })

    output$design_matrix <- renderPrint({
      design_matrix(obj, input$which_model_design)
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
      saved_plots_and_tables$scatter_plt <- p
      p
    })

    output$download_scatter_plt <- downloadHandler(
        filename = function() {
          "scatter_plot.pdf"
          },
        content = function(file) {
            ggsave(file,
              saved_plots_and_tables$scatter_plt,
              width = user_settings$save_width,
              height = user_settings$save_height,
              units = "cm")
    })

    output$scatter_vars <- renderPlot({
      # NB: inherence the test from the QQ plot
      test_type <- input$settings_test_type

      current_test <- input$which_test_qq
      if (is.null(current_test)) {
        possible_tests <- list_tests(obj, test_type)
        if (test_type == 'wt') {
          possible_tests <- possible_tests[[input$which_model_qq]]
        }
        current_test <- possible_tests[1]
      }

      scatter_var_plt <- plot_vars(obj,
        current_test,
        test_type,
        which_model = input$which_model_qq,
        point_alpha = input$scatter_alpha,
        highlight = rv_scatter$highlight_vars,
        sig_level = input$max_fdr
        )
      saved_plots_and_tables$scatter_var_plt <- scatter_var_plt
      scatter_var_plt
    })

    output$download_scatter_var_plt <- downloadHandler(
        filename = function() {
          "scatter_var_plot.pdf"
          },
        content = function(file) {
          ggsave(file,
            saved_plots_and_tables$scatter_var_plt,
            width = user_settings$save_width,
            height = user_settings$save_height,
            units = "cm")
    })

    output$scatter_brush_table <- renderDataTable({
      test_type <- input$settings_test_type

      current_test <- input$which_test_qq
      if ( is.null(current_test) ||
        !test_exists(obj, current_test, test_type, input$which_model_qq)) {
        possible_tests <- list_tests(obj, test_type)
        if (test_type == 'wt') {
          possible_tests <- possible_tests[[input$which_model_qq]]
        }
        current_test <- possible_tests[1]
      }

      res <- NULL
      sr <- sleuth_results(obj,
        current_test,
        test_type,
        input$which_model_qq,
        rename_cols = FALSE,
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
      if (!is.null(res)) {
            saved_plots_and_tables$scatter_table <- res
            output$download_scatter_table_button <- renderUI({
            div(align = "right", style = "margin-right:15px; margin-top:10px; margin-bottom:10px",
                    downloadButton("download_scatter_table", "Download Table"))
            })
      }
      res
    })

    output$download_scatter_table <- downloadHandler(
      filename = function() {
        "scatter_table.csv"
        },
      content = function(file) {
         write.csv(saved_plots_and_tables$scatter_table, file)
    })

    ###
    output$samp_heat_plt <- renderPlot({
      samp_heat_plt <- plot_sample_heatmap(obj, use_filtered = input$samp_heat_filt,
                                           annotation_cols = input$samp_heat_covars,
                                           cluster_bool = input$samp_heat_cluster)
      saved_plots_and_tables$samp_heat_plt <- samp_heat_plt #this is a gtable object
      samp_heat_plt
    })

    output$download_samp_heat_plt <- downloadHandler(
        filename = function() {
          "samp_heat_plot.pdf"
          },
        content = function(file) {
            ggsave(file,
              saved_plots_and_tables$samp_heat_plt,
              width = user_settings$save_width,
              height = user_settings$save_height,
              units = "cm")
    })

    ###
    output$pca_plt <- renderPlot({

      color_by <- ifelse(is.null(input$color_by), NULL,
        as.character(input$color_by))

      pca_plt <- plot_pca(obj,
        pc_x = as.integer(input$pc_x),
        pc_y = as.integer(input$pc_y),
        text_labels = as.logical(input$text_labels),
        color_by = color_by,
        use_filtered = input$pca_filt,
        units = input$pca_units,
        point_size = input$pca_point_size
        )

      saved_plots_and_tables$pca_plt <- pca_plt
      pca_plt
    })

    output$download_pca_plt <- downloadHandler(
      filename = function() {
        "pca_plot.pdf"
        },
      content = function(file) {
         ggsave(file, saved_plots_and_tables$pca_plt,
           width = user_settings$save_width,
           height = user_settings$save_height,
           units = "cm")
    })

    output$fld_plt <- renderPlot({
      plot_fld(obj, input$fld_sample)
    })

    output$bias_weights_table <- renderDataTable({
      bias_table(obj, input$bias_sample)
    })

   ### Plot pc loadings
    output$plt_pc_loadings <- renderPlot({
      plot_loadings(obj,
        use_filtered = input$pc_filt,
        pc_count = as.integer(input$pc_count),
        scale = as.logical(input$scale),
        sample = input$sample,
        units = input$pca_units,
        pc_input = as.integer(input$pc_input),
        pca_loading_abs = input$pca_loading_abs
        )
    })

    ###plot pc variance
    output$plt_pc_var <- renderPlot({
      plot_pc_variance(obj,
        use_filtered = input$pc_filt,
        pca_number = 5,
        scale = input$scale,
        units = input$pca_units,
        PC_relative = 1
        )
    })


    ### MV plot
    output$mv_plt <- renderPlot({
      mv_plt <- plot_mean_var(obj)
      saved_plots_and_tables$mv_plt <- mv_plt
      mv_plt
    })

    output$download_mv_plt <- downloadHandler(
      filename = function() {
        "mv_plot.pdf"
        },
      content = function(file) {
         ggsave(file,
           saved_plots_and_tables$mv_plt,
           width = user_settings$save_width,
           height = user_settings$save_height,
           units = "cm")
    })

    ### MA
    rv_ma <- reactiveValues(
      highlight_vars = NULL
      )

    output$which_beta_ctrl_ma <- renderUI({
      possible_tests <- list_tests(obj, input$settings_test_type)
      possible_tests <- possible_tests[[input$which_model_ma]]
      selectInput('which_beta_ma', 'beta: ', choices = possible_tests,
        selected = possible_tests[1])
    })

    output$ma <- renderPlot({

      current_test <- input$which_beta_ma
      if (is.null(current_test)) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_ma]]
        current_test <- possible_tests[1]
      }
      ma_plt <- plot_ma(obj, current_test,
        input$settings_test_type,
        input$which_model_ma,
        sig_level = input$max_fdr,
        point_alpha = input$ma_alpha
        )
      saved_plots_and_tables$ma_plt <- ma_plt
      ma_plt
    })

    output$download_ma_plt <- downloadHandler(
      filename = function() {
        "ma_plot.pdf"
      },
      content = function(file) {
         ggsave(file,
           saved_plots_and_tables$ma_plt,
           width = user_settings$save_width,
           height = user_settings$save_height,
           units = "cm")
      }
    )

    output$vars <- renderPlot({
      wb <- input$which_beta_ma
      if (is.null(wb)) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_ma]]
        wb <- possible_tests[1]
      }
      ma_var_plt <- plot_vars(obj,
        test = wb,
        test_type = input$settings_test_type,
        which_model = input$which_model_ma,
        point_alpha = input$ma_alpha,
        highlight = rv_ma$highlight_vars,
        sig_level = input$max_fdr
        )
      saved_plots_and_tables$ma_var_plt <- ma_var_plt
      ma_var_plt
    })

    output$download_ma_var_plt <- downloadHandler(
      filename = function() {
        "ma_var_plot.pdf"
        },
      content = function(file) {
         ggsave(file,
           saved_plots_and_tables$ma_var_plt,
           width = user_settings$save_width,
           height = user_settings$save_height,
           units = "cm")
      }
    )

    #observe
    output$ma_brush_out <- renderDataTable({
      wb <- input$which_beta_ma
      if ( is.null(wb) ) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_ma]]
        wb <- possible_tests[1]
      }
      res <- sleuth_results(obj,
        wb,
        input$settings_test_type,
        input$which_model_ma, rename_cols = FALSE, show_all = FALSE)
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
      if (!is.null(res)) {
        saved_plots_and_tables$ma_table <- res
        output$download_ma_table_button <- renderUI({
            div(align = "right", style = "margin-right:15px; margin-top:10px; margin-bottom:10px",
                    downloadButton("download_ma_table", "Download Table"))
        })
      }
      res
    })

    output$download_ma_table <- downloadHandler(
      filename = function() {
        "ma_table.csv"
        },
      content = function(file) {
         write.csv(saved_plots_and_tables$ma_table, file)
      }
    )

    ### DE table

    output$which_beta_ctrl_de <- renderUI({
      # poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model_de]])
      possible_tests <- list_tests(obj, input$settings_test_type)
      result <- NULL
      if ( input$settings_test_type == 'wt' ) {
        possible_tests <- possible_tests[[input$which_model_de]]
        result <- selectInput('which_beta_de', 'beta: ', choices = possible_tests)
      } else {
        result <- selectInput('which_beta_de', 'test: ', choices = possible_tests)
      }

      result
    })
    output$table_type <- renderUI({
      if(!is.null(obj$target_mapping)) {
        selectInput('pop_genes', label = HTML('table type: ',
        '<button onclick="tableType()">?</button>',
        '<script> function tableType() {',
        'alert("The transcript table shows in-depth summary statistics and test results for each individual',
        'transcript. The gene table, on the other hand, shows significance statistics for genes, each of which',
        'which may contain many transcripts.',
        '");',
        '} </script>'),
          choices = list('transcript table' = 1, 'gene table' = 2),
          selected = 1)
      }
    })

    output$group_by <- renderUI({
      if (!is.null(input$pop_genes) && input$pop_genes == 2) {
        selectInput('mappingGroup',
          label = HTML('group by: ',
          '<button onclick="mappingGroup()">?</button>',
          '<script> function mappingGroup() {',
          'alert("View genes by their Ensemble gene ID (ens_gene), or by their external gene name (ext_gene).',
          'These gene names are from the gene annotation you used to add genes to your sleuth object.");',
          '} </script>'),
          choices = names(obj$target_mapping)[2:length(names(obj$target_mapping))])
      }
    })

    output$de_dt <- renderDataTable({
      wb <- input$which_beta_de
      if (is.null(wb)) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_de]]
        wb <- possible_tests[1]
      }

      if (!is.null(input$mappingGroup) && (input$pop_genes == 2)) {
          mg <- input$mappingGroup
          test_table <- sleuth_gene_table(obj, wb, input$settings_test_type,
            input$which_model_de, mg)
            saved_plots_and_tables$test_table <- test_table
            test_table
      } else {
          test_table <- sleuth_results(obj, wb, input$which_model_de)
          saved_plots_and_tables$test_table <- test_table
          test_table
      }
    })

    output$download_test_table <- downloadHandler(
      filename = function() {
        "test_table.csv"
        },
      content = function(file) {
         write.csv(saved_plots_and_tables$test_table, file)
      }
    )

    output$test_control_de <- renderUI({
      possible_tests <- list_tests(obj, input$settings_test_type)
      selectInput('which_test_de', 'test: ', choices = possible_tests)
    })

    output$lrt_de_dt <- renderDataTable({
      which_test <- input$which_test_de
      if ( is.null(which_test) ) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        which_test <- possible_tests[1]
      }
      if (!is.null(input$mapping_group_lrt) && (input$pop_genes_lrt == 2)) {
          mg <- input$mapping_group_lrt
          sleuth_gene_table(obj, which_test, input$settings_test_type,
            which_group = mg)
      } else {
          sleuth_results(obj, which_test, input$settings_test_type)
      }
    })

    output$table_type_lrt <- renderUI({
      if(!is.null(obj$target_mapping)) {
        selectInput('pop_genes_lrt', label = HTML('table type: ',
        '<button onclick="tableType()">?</button>',
        '<script> function tableType() {',
        'alert("The transcript table shows in-depth summary statistics and test results for each individual',
        'transcript. The gene table, on the other hand, shows significance statistics for genes, each of which',
        'which may contain many transcripts.");',
        '} </script>'),
          choices = list('transcript table' = 1, 'gene table' = 2),
          selected = 1)
      }
    })
    output$group_by_lrt <- renderUI({
      if (!is.null(input$pop_genes_lrt) && input$pop_genes_lrt == 2) {
        selectInput('mapping_group_lrt',
          label = HTML('group by: ',
          '<button onclick="mappingGroup()">?</button>',
          '<script> function mappingGroup() {',
          'alert("View genes by their Ensemble gene ID (ens_gene), or by their external gene name (ext_gene).',
          'These gene names are from the gene annotation you used to add genes to your sleuth object.");',
          '} </script>'),
          choices = names(obj$target_mapping)[2:length(names(obj$target_mapping))])
      }
    })


    ### Transcript viewer, transcript mode
    output$transcript_view <- renderUI({
      if (!obj$gene_mode) {
        layout <- list(
          fluidRow(
            column(12,
                   p(
                     h3('transcript view'),
                     'Boxplots of transcript abundances showing technical variation in each sample.')
            ),
            offset = 1),
          fluidRow(
            column(4,
                   textInput('bs_var_input',
                             label = HTML('transcript: ',
                                          '<button onclick="bs_var_input()">?</button>',
                                          '<script> function bs_var_input() {',
                                          'alert("Enter the target_id of a transcript here to view a boxplot of its technical variation.',
                                          'You can find target_ids in the test table under Analyses. The most significant transcript',
                                          'by q-value is already entered.");',
                                          '} </script>'),
                             value = (obj$tests[[1]][[1]])[
                               order(obj$tests[[1]][[1]]$qval),]$target_id[1])
            ),
            column(4,
                   selectInput('bs_var_color_by',
                               label = HTML('color by: ',
                                            '<button onclick="bs_var_color_by()">?</button>',
                                            '<script> function bs_var_color_by() {',
                                            'alert("Color the box plots by different column names from your design file.',
                                            'Doing so can help explain which traits of a given transcript best explain the',
                                            'difference in counts of that transcript across samples.");',
                                            '} </script>'),
                               choices = c(NULL, poss_covars), selected = NULL)
            ),
            column(3,
                   selectInput('bs_var_units',
                               label = HTML('units: ',
                                            '<button onclick="bs_var_units()">?</button>',
                                            '<script> function bs_var_units() {',
                                            'alert("tpm: Transcripts Per Million\\nest_counts (transcript mode only): Counts/Number of aligned reads\\n',
                                            'scaled_reads_per_base (gene mode only): The average number of reads mapping to each base across the whole gene\\n',
                                            '\\nFor a better understanding of the units, see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/");',
                                            '} </script>'),
                               choices  = (if (length(obj$bs_quants) == 0) { c('N/A') }
                                           else { names(obj$bs_quants[[1]]) } ),
                               selected = (if (length(obj$bs_quants) == 0) { 'N/A' }
                                           else { names(obj$bs_quants[[1]])[1] })
                   ))
          ),
          fluidRow(HTML('&nbsp;&nbsp;&nbsp;'), actionButton('bs_go', 'view')),
          fluidRow(plot <- (if (length(obj$bs_quants) == 0) {
            HTML('&nbsp&nbsp&nbsp&nbsp You need to run sleuth with at least ',
                 'one of extra_bootstrap_summary or ',
                 'read_bootstrap_tpm to use this feature.<br>') }
            else { plotOutput('bs_var_plt') }
          )),
          fluidRow(
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_bs_var_plt", "Download Plot"))
          )
        )
      } else {
        layout <- list(
          HTML('&nbsp&nbsp&nbsp&nbspThis view is disabled because your sleuth',
               ' object is in gene mode.<br>&nbsp&nbsp&nbsp&nbspTo view',
               ' transcripts, rerun sleuth without the \'aggregation_column\'',
               ' option.<br>&nbsp&nbsp&nbsp&nbspFor more information, see the ',
               '<a href = "https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html">',
               'sleuth getting started guide</a>.')
        )
        layout
      }

    })

    bs_var_text <- eventReactive(input$bs_go, {
        input$bs_var_input
    })

    output$bs_var_plt <- renderPlot({
      saved_plots_and_tables$bs_var_plt <- plot_bootstrap(
        obj, bs_var_text(),
        units = input$bs_var_units,
        color_by = input$bs_var_color_by)
      saved_plots_and_tables$bs_var_plt
    })

    output$download_bs_var_plt <- downloadHandler(
      filename = function() {
        "transcript_bootstrap_vars.pdf"
      },
      content = function(file) {
        ggsave(file,
               saved_plots_and_tables$bs_var_plt,
               width = user_settings$save_width,
               height = user_settings$save_height,
               units = "cm")
      })

    ### Gene viewer, gene mode

    bsg_var_text <- eventReactive(input$bsg_go, {
      if (obj$gene_mode) {
        gene_from_gene(obj,
                       input$bsg_gene_colname,
                       input$bsg_var_input)
      }
    })

    output$bsg_var_plt <- renderPlot({
      saved_plots_and_tables$bs_var_plt <- plot_bootstrap(
        obj, bsg_var_text(),
        units = input$bsg_var_units,
        color_by = input$bsg_var_color_by)
      saved_plots_and_tables$bs_var_plt
    })

    output$download_bsg_var_plt <- downloadHandler(
      filename = function() {
        "gene_bootstrap_vars.pdf"
      },
      content = function(file) {
        ggsave(file,
               saved_plots_and_tables$bs_var_plt,
               width = user_settings$save_width,
               height = user_settings$save_height,
               units = "cm")
      })

    ### Gene Viewer, transcript mode
    # the name of the gene supplied
    gv_var_text <- eventReactive(input$gv_go, {
      if (!is.null(obj$target_mapping) & !obj$gene_mode) {
          input$gv_var_input
      }
    })

    output$gv_gene_column <- renderUI({
      if (!is.null(obj$target_mapping) & !obj$gene_mode) {
        selectInput('gv_gene_colname',
          label = HTML('genes from: ',
                       '<button onclick="gv_gene_colname()">?</button>',
                       '<script> function gv_gene_colname() {',
                       'alert("Choose a column name from your target mapping data frame.',
                       'This allows you to use alternative gene identifiers (e.g. HGNC symbol)',
                       'to select a gene.',
                       'The second column in target mapping is the default");',
                       '} </script>'),,
          choices = setdiff(names(obj$target_mapping), c('target_id')))
      }
    })

    output$which_beta_ctrl_gv <- renderUI({
      poss_tests <- tests(models(obj, verbose = FALSE)[[input$which_model_de]])
      selectInput('which_beta_gv', 'beta: ', choices = poss_tests)
    })


    gv_var_list <- reactive({
      test_type <- input$settings_test_type

      current_test <- input$which_test_qq
      if ( is.null(current_test) ||
        !test_exists(obj, current_test, test_type, input$which_model_qq)) {
        possible_tests <- list_tests(obj, test_type)
        if (test_type == 'wt') {
          possible_tests <- possible_tests[[input$which_model_qq]]
        }
        current_test <- possible_tests[1]
      }

      if (!obj$gene_mode) {
        transcripts_from_gene(obj,
                              current_test,
                              test_type,
                              input$which_model_qq,
                              input$gv_gene_colname,
                              gv_var_text())
      }
    })

    output$gv_var_plts <- renderUI({
      gv_plot_list <- lapply(1:input$gv_maxplots,
        function(i) {
          gv_plotname <- paste("plot", i, sep = "")
          plotOutput(gv_plotname)
        })

        do.call(tagList, gv_plot_list)
      })

      for (i in 1:15) {
        local({
          my_i <- i
          gv_plotname <- paste("plot", my_i, sep = "")
            output[[gv_plotname]] <- renderPlot({
               if (!is.null(obj$target_mapping) && !obj$gene_mode &&
                   !is.na(gv_var_list()[my_i])) {
                 plot_bootstrap(obj,
                   gv_var_list()[my_i],
                   units = input$gv_var_units,
                   color_by = input$gv_var_color_by)
              }
             })
         })
        }

    output$no_genes_message <- renderUI({
      if (is.null(obj$target_mapping)) {
        HTML('&nbsp&nbsp&nbsp&nbspYou need to add genes to your sleuth object to use the gene viewer.<br>',
        '&nbsp&nbsp&nbsp&nbspTo add genes to your sleuth object, see the ',
        '<a href = "https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html">sleuth getting started guide</a>.')
      }
    })

    ### Heat Map
    hm_transcripts <- eventReactive(input$hm_go, {
            unlist(strsplit(input$hm_transcripts, " +"))
    })

    default_hm_plot_height <- 400

    hm_plot_height <- function() {
        if (length(hm_transcripts()) > 5) {
            length(hm_transcripts()) * 60
        } else {
            default_hm_plot_height
        }
    }

    hm_func <- eventReactive(input$hm_go, {
        input$hm_trans
    })

    output$hm_plot <- renderPlot ({
        saved_plots_and_tables$hm_plt <- plot_transcript_heatmap(obj,
          hm_transcripts(),
          input$hm_units, hm_func(),
          cluster_transcripts = input$hm_cluster,
          offset = input$hm_offset,
          annotation_cols = input$hm_covars)

        output$download_hm_plt_button <- renderUI({
          div(
            align = "right",
            style = paste("margin-right:15px; margin-top:",
              (hm_plot_height() - default_hm_plot_height + 15),
              "px; margin-bottom:10px",
              sep = ""),
            downloadButton("download_hm_plt", "Download Plot"))
          })
         saved_plots_and_tables$hm_plt

    }, height = hm_plot_height)


    output$download_hm_plt <- downloadHandler(
      filename = function() {
        "heat_map.pdf"
        },
      content = function(file) {
       ggsave(file,
         saved_plots_and_tables$hm_plt,
         width = user_settings$save_width,
         height = user_settings$save_height,
         units = "cm")
    })
    ### Volcano Plot
    output$which_beta_ctrl_vol <- renderUI({
      possible_tests <- list_tests(obj, input$settings_test_type)
      possible_tests <- possible_tests[[input$which_model_vol]]
      selectInput('which_beta_vol', 'beta: ', choices = possible_tests,
        selected = possible_tests[1])
    })

    output$vol <- renderPlot({
      which_test <- input$which_beta_vol
      if (is.null(which_test)) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_vol]]
        which_test <- possible_tests[1]
      }
      volcano_plt <- plot_volcano(obj, which_test,
        input$settings_test_type,
        input$which_model_vol,
        sig_level = input$max_fdr_vol,
        point_alpha = input$vol_alpha
        )
        saved_plots_and_tables$volcano_plt <- volcano_plt
        volcano_plt
    })

    output$download_volcano_plt <- downloadHandler(
      filename = function() {
        "volcano.pdf"
        },
      content = function(file) {
       ggsave(file,
         saved_plots_and_tables$volcano_plt,
         width = user_settings$save_width,
         height = user_settings$save_height,
         units = "cm")
    })

    #vol_observe
    output$vol_brush_out <- renderDataTable({
      wb <- input$which_beta_vol
      if (is.null(wb)) {
        possible_tests <- list_tests(obj, input$settings_test_type)
        possible_tests <- possible_tests[[input$which_model_vol]]
        wb <- possible_tests[1]
      }

      res <- sleuth_results(obj, wb, input$settings_test_type,
        input$which_model_vol, rename_cols = FALSE, show_all = FALSE)
      res <- dplyr::mutate(res, Nlog10_qval = -log10(qval))
      if (!is.null(input$vol_brush)) {
        res <- brushedPoints(res, input$vol_brush, xvar = 'b', yvar = 'Nlog10_qval')
        res$Nlog10_qval <- NULL
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

      if (!is.null(res)) {
          saved_plots_and_tables$volcano_table <- res
          output$download_volcano_table_button <- renderUI({
          div(align = "right", style = "margin-right:15px; margin-top:10px",
                  downloadButton("download_volcano_table", "Download Table"))
          })
      }

      res
    })

    output$download_volcano_table <- downloadHandler(
      filename = function() {
        "volcano_table.csv"
        },
      content = function(file) {
         write.csv(saved_plots_and_tables$volcano_table, file)
    })

  }

  shinyApp(ui = p_layout, server = server_fun, options = options)
}

#' @export
enclosed_brush <- function(df, brush) {
  df <- as_df(df)
  xvar <- brush$mapping$x
  yvar <- brush$mapping$y

  xbool <- brush$xmin <= df[[xvar]] & df[[xvar]] <= brush$xmax
  ybool <- brush$ymin <= df[[yvar]] & df[[yvar]] <= brush$ymax

  df[xbool & ybool, ]
}

#' settings for sleuth_live
#'
#' This is a helper function for setting preferences in sleuth live.
#' Currently, this is somewhat limited, but it will be expanded in the future.
#'
#' @param test_type either 'wt' for the wald test or 'lrt' for the likelihood ratio test
#' @return a named to list with the options and settings
#' @export
sleuth_live_settings <- function(test_type = 'wt') {
  stopifnot( is(test_type, 'character') )

  result <- list()
  result$test_type <- test_type

  result
}

#' deploy a sleuth object
#'
#' prepare a sleuth object to be deployed in a shiny application.
#'
#' creates a directory \code{path} and creates a valid shiny application.
#'
#' \itemize{
#'  \item saves a sleuth object using \code{\link{sleuth_save}}
#'  \item creates a file \code{app.R} loading the sleuth object and calling \code{\link{sleuth_live}}
#' }
#'
#' @param obj a \code{sleuth} object
#' @param base_dir the base directory in which to save a shiny application
#' @param overwrite if \code{TRUE}, overwrite everything at \code{bayes_dir}
#' @seealso \code{\link{sleuth_save}}, \code{\link{sleuth_load}}, \code{\link{sleuth_live}}
#' @export
sleuth_deploy <- function(obj, base_dir, overwrite = FALSE) {
  obj_name <- 'so.rds'

  if (!overwrite) {
    if (file.exists(file.path(base_dir, obj_name)) &&
      file.exists(file.path(base_dir, 'app.R'))) {
      stop('a sleuth shiny object already exists at this location. Either specify a different directory or set "overwrite = TRUE"')
    }
  }

  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE, mode = '0755')

  sleuth_save(obj, file = file.path(base_dir, obj_name))
  commands <- paste(
    "library(sleuth)",
    "",
    paste0("so <- sleuth_load('", obj_name, "')"),
    "sleuth_live(so)",
    sep = "\n")

  write(commands, file = file.path(base_dir, 'app.R'))

  invisible(NULL)
}
