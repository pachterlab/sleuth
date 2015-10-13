plot_loadings <- function(obj, 
  use_filtered = TRUE,
  gene = NULL,
  pc_count = NULL,
  bool = FALSE,
  absolute = TRUE, #edit this
  ...) {

  mat <- obj
  
  pca_calc <- prcomp(mat, scale = bool)
  #transpose
  loadings <- t(pca_calc$rotation)

  if (!is.null(gene)) {
    loadings <- loadings[, gene]
  } else {
    title <- colnames(loadings)[1]
    loadings <- loadings[, 1]
  }

  if (!is.null(pc_count)) {
    loadings <- loadings[1:pc_count]
  } else {  
    loadings <- loadings[1:5]
  }

  if (absolute) {
    loadings <- abs(loadings)
  }

  #sort loadings vector to obtain highest contribution
  loadings <- sort(loadings, decreasing = TRUE)
  names <- names(loadings)

  dat <- data.frame(pc = names, loadings = loadings)
  dat$pc <- factor(dat$pc, levels = unique(dat$pc))

  p <- ggplot(dat, aes(x = pc, y = loadings)) 
  p <- p + geom_bar(stat = "identity")
  p <- p + xlab("Principal Components") + ylab("Contribution")

  if (!is.null(gene)) {
    p <- p + ggtitle(gene)
  } else {
    p <- p + ggtitle(title)
  }

  p

}