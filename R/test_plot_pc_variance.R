plot_pc_variance <- function(obj, 
  use_filtered = TRUE,
  pca_number = NULL,
  bool = FALSE,
  ...) {

  mat <- obj
  pca_calc <- prcomp(mat, scale = bool) #PCA calculations 

  #computation
  eigenvalues <- (pca_calc$sdev)^2  
  variance <- eigenvalues*100/sum(eigenvalues)
  cum_var <- cumsum(variance)

  #did not throw error messages
  if (!is.null(pca_number)) {
    colsize <- pca_number
    cum_var <- cum_var[1:pca_number]
  } else {
    colsize <- 5
    cum_var <- cum_var[1:5]
  }

  pc_asdf <- data.frame(PC_count = 1:colsize, cum_var = cum_var)

  #change to ggplot
  p <- ggplot(pc_asdf, aes(x = PC_count, y = cum_var)) + geom_bar(stat = "identity")
  p <- p + scale_x_continuous(breaks = 1:length(eigenvalues))
  p <- p + ylab("% of Variance") + xlab("Principal Components")

  p

}