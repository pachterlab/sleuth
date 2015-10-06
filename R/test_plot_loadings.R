plot_loadings <- function(obj, 
  use_filtered = TRUE,
  gene = NULL,
  pca_number = NULL,
  bool = FALSE,
  ...) {

  mat <- obj
  
  pca_calc <- prcomp(mat, scale = bool)
  #assume no sample has redundant names
  if (!is.null(gene)){
    #accessor function for row names 
    if (!is.null(pca_number)) {
      loadings <- pca_calc$rotation[gene,1:pca_number]
      colsize <- pca_number
    } else {
      loadings <- pca_calc$rotation[gene, 1:5]
      colsize <- 5
    }
  } else {
    if (!is.null(pca_number)) {
      loadings <- pca_calc$rotation[1,1:pca_number]
      colsize <- pca_number
    } else {
      loadings <- pca_calc$rotation[1,1:5]
      colsize <- 5
    }
  }
  dat <- data.frame(PC_count = 1:colsize, loadings = loadings)

  p <- ggplot(dat, aes(x = PC_count, y = loadings)) 
  p <- p + geom_point(shape = 1)
  p <- p + xlab("Principal Components") + ylab("Loadings")

  if (!is.null(gene)) {
    p <- p + ggtitle(gene)
  }

  p

}