#give a principal component, tells you which contribute the most
#give a gene, tells you which PC's it contributes to the most


plot_loadings <- function(obj, 
  use_filtered = TRUE,
  gene = NULL, #string please
  PC = NULL, #apparently if u type in a string or an integer (corresponding to the PC), they're both ok
  pc_count = NULL,
  bool = FALSE,
  absolute = TRUE, #edit this
  ...) {
  #make sure proper arguments are given
  stopifnot( !is.null(gene) && !is.null(PC) )

  mat <- obj
  
  pca_calc <- prcomp(mat, scale = bool)
  #transpose
  loadings <- pca_calc

  #given a gene
  if (!is.null(gene)) {
    loadings <- pca_calc$x[gene,]
    if (absolute) {
      loadings <- abs(loadings)
    }
    loadings <- sort(loadings, decreasing = TRUE)
    names <- names(loadings)
  }

  #given a PC, which samples contribute the most?
  if (!is.null(PC)) {
    loadings <- pca_calc$x[,PC]
    if (absolute) {
      loadings <- abs(loadings)
    }
    loadings <- sort(loadings, decreasing = TRUE)
    names <- names(loadings)
  }

  if (!is.null(pc_count)) {
      loadings <- loadings[1:pc_count]
      names <- names[1:pc_count]
    } else {
      loadings <- loadings[1:5]
      names <- names[1:5]
    }

  dat <- data.frame(pc = names, loadings = loadings)
  dat$pc <- factor(dat$pc, levels = unique(dat$pc))

  p <- ggplot(dat, aes(x = pc, y = loadings)) 
  p <- p + geom_bar(stat = "identity")
  p <- p + xlab("Principal Components") + ylab("Contribution Scores")

  #logistics of graph
  print(is.numeric(PC))
  if (is.numeric(PC)) {
    PC <- paste0("PC ", PC)
  }

  if (!is.null(gene)) {
    p <- p + ggtitle(gene)
  } else {
    p <- p + ggtitle(PC)
  }

  p

}