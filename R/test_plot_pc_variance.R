plot_pc_variance <- function(obj, 
  use_filtered = TRUE,
  pca_number = NULL,
  ...) {
  stopifnot( is(obj, 'sleuth') ) 

  mat <- NULL
  if (use_filtered) {
    mat <- spread_abundance_by(obj$obs_norm_filt, units)
  } else {
    mat <- spread_abundance_by(obj$obs_norm, units)
  }

  pca_calc <- prcomp(mat, scale = TRUE) #PCA calculations 

  #computation
  eigenvalues <- (pca_calc$sdev)^2  
  variance <- eigenvalues*100/sum(eigenvalues)
  cum_var <- cumsum(variance)

  #add PC columns to take ranges

  pc_asdf <- as_df(eigenvalues = eigenvalues, variance = variance, 
                      cumulative_variance = cum_var) #put PCA loadings into a data frame 


  #change to ggplot
  if (!is.null(pca_number)) {
    p <- barplot(pc_asdf[,2], names.arg = 1:pca_number, #set the x,y graph coordinate names
                    main = "Variances",
                    xlab = "Principal Components",
                    ylab = "% of Variances",
                    col = "dodgerblue")
  } else {
    p <- barplot(pc_asdf[,2], names.arg = 1:5, #set the x,y graph coordinate names
                    main = "Variances",
                    xlab = "Principal Components",
                    ylab = "% of Variances",
                    col = "dodgerblue")
  }

  p
}