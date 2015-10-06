plot_loadings <- function(obj, 
  use_filtered = TRUE,
  gene = NULL,
  pca_number = NULL,
  ...) {

  mat <- obj
  pca_calc <- prcomp(mat, scale = TRUE)
  #assume no sample has redundant names
  if (!is.null(gene)){
    #accessor function for row names 
    if (!is.null(pca_number)) {
      loadings <- pca_calc$rotation[gene,1:pca_number]
    }
    loadings <- pca_calc$rotation[gene, 1:5]
  } 
  if (!is.null(pca_number)) {
    loadings <- pca_calc$rotation[1,1:pca_number]
  }
  loadings <- pca_calc$rotation[1,1:5]

  dat <- as.data.frame(loadings)
  #change aes_string
  p <- ggplot(dat, aes(x = rownames(dat), y = loadings)) + geom_point(shape = 1)
  p <- p + geom_point(shape = 1)
  p <- p + xlab("Principal Components") + ylab("Loadings")
  if (!is.null(gene)) {
    p <- p + ggtitle(gene)
  }

  p

}