plot_loadings <- function(obj, 
  use_filtered = TRUE,
  PC = NULL,
  gene_count = NULL,
  bool = FALSE,
  ...) {

  mat <- obj
  
  pca_calc <- prcomp(mat, scale = bool)
  #assume no sample has redundant names
  if (!is.null(PC)){
    #accessor function for row names 
    if (!is.null(gene_count)) {
      loadings <- pca_calc$rotation[1:gene_count, PC]
      colsize <- gene_count
    } else {
      loadings <- pca_calc$rotation[1:5, PC]
      colsize <- 5
    }
  } else {
    if (!is.null(gene_count)) {
      loadings <- pca_calc$rotation[1:gene_count,1]
      colsize <- pca_number
    } else {
      loadings <- pca_calc$rotation[1:5, 1]
      colsize <- 5
    }
  }
  dat <- data.frame(samples = 1:colsize, loadings = loadings)

  p <- ggplot(dat, aes(x = samples, y = loadings)) 
  p <- p + geom_point(shape = 1)
  p <- p + xlab("Samples") + ylab("% Contribution")

  if (!is.null(PC)) {
    p <- p + ggtitle(PC)
  } else {
    p <- p + ggtitle()
  }

  p

}