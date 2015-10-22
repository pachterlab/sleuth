#how much additional information is revealed

plot_pc_variance <- function(obj, 
  use_filtered = TRUE,
  pca_number = NULL,
  bool = FALSE,
  PC_relative = NULL, #this is an integer
  ...) {


  mat <- obj

  pca_calc <- prcomp(mat, scale = bool) #PCA calculations 

  #computation
  eigenvalues <- (pca_calc$sdev)^2  
  var <- eigenvalues*100/sum(eigenvalues)
  var2 <- eigenvalues*100/sum(eigenvalues) #because i suck at coding
  
  #from here to ....
  if (!is.null(pca_number)) {
    colsize <- pca_number
    var <- var[1:pca_number]
  } else {
    colsize <- 5 #default 5
    var <- var[1:5] #default 5
  }
  pc_asdf <- data.frame(PC_count = 1:colsize, var = var) #order here matters
  #...here is for comparison of variance for subsequent PC's

  #here is comparison per given principal component
  if(!is.null(PC_relative)) {
    pc_asdf <- data.frame(PC_count = 1:length(eigenvalues), var = var2) #because i suck at coding (var2 stupid shit)
    pc_asdf <- pc_asdf[PC_relative:nrow(pc_asdf),] #new data frame

    #if user wants to give some PCA count to graph (default is 5 or until the end)
    if (!is.null(pca_number) && (PC_relative + pca_number <= length(eigenvalues))) { #check if it does not overflow data frame
      pc_asdf <- pc_asdf[1:pca_number,] #new data frame if user wants to give a pca count number
    } else if (PC_relative + 5 >= length(eigenvalues)) {
      pc_asdf <- pc_asdf[1:nrow(pc_asdf),] 
    }
  } 

  p <- ggplot(pc_asdf, aes(x = PC_count, y = var)) + geom_bar(stat = "identity")
  p <- p + scale_x_continuous(breaks = 1:length(eigenvalues))
  p <- p + ylab("% of Variance") + xlab("Principal Components")

  p
  
}