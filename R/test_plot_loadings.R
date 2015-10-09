plot_loadings <- function(obj, 
  use_filtered = TRUE,
  PC = NULL,
  gene_count = NULL,
  bool = FALSE,
  ...) {

  mat <- obj
  
  pca_calc <- prcomp(mat, scale = bool)
  #assume no sample has redundant names
  if (!is.null(PC)) {
    loadings <- pca_calc$rotation[, PC]
  } else {
    loadings <- pca_calc$rotation[, 1]
  }

  if (!is.null(gene_count)) {
    loadings <- loadings[1:gene_count]
  } else {
    loadings <- loadings[1:5]
  }

  #sort loadings vector to obtain highest contribution
  loadings <- sort(loadings, decreasing = TRUE)
  names <- names(loadings)

  dat <- data.frame(samples = names, loadings = loadings)
  dat$samples <- factor(dat$samples, levels = unique(dat$samples))

  #debugging
  #print(dat$samples)
  #print(dat)

  p <- ggplot(dat, aes(x = samples, y = loadings)) 
  p <- p + geom_point(shape = 1)
  p <- p + xlab("Samples") + ylab("Contribution")

  if (!is.null(PC)) {
    p <- p + ggtitle(PC)
  } else {
    p <- p + ggtitle("PC1")
  }

  p

}