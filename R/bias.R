#make sure this is in the current directory of samples
#make sure bias_table is already loaded in

#test for local machine

#setwd("~/Documents/Research/sleuth/ellahi/results")


#assumptions: 
#1. there are an equal number of hexamers for each base
#2. they are in the order of ACGT

library(ggplot2)


bias_graph <- function(obj,  
	sample = "",
	...) {
	stopifnot( is(obj, 'sleuth'))
  
	#sample df construction
	bt <- bias_table(obj, sample) 
	bt <- data.frame(hexamer = bt$hexamer, bias_weights = bt$bias_weights)
	bt$hexamer <- sapply(bt$hexamer, as.character)
	norm_count <- nrow(bt)
	
	#set up
	A_vector <- c(0,0,0,0,0,0) 
	C_vector <- c(0,0,0,0,0,0) 
	G_vector <- c(0,0,0,0,0,0) 
	T_vector <- c(0,0,0,0,0,0) 
	
	for (i in 1:nrow(bt)) {
		iter <- strsplit(bt$hexamer[i], "")[[1]]
		bias_append <- bt$bias_weights[i]
		#index 1-6
		for (x in 1:6) {
			if (iter[x] == "A") {
				A_vector[x] <- A_vector[x] + bias_append
			}
			if (iter[x] == "C") {
				C_vector[x] <- C_vector[x] + bias_append
			}
			if (iter[x] == "G") {
				G_vector[x] <- G_vector[x] + bias_append
			}
			if (iter[x] == "T") {
				T_vector[x] <- T_vector[x] + bias_append
			}
		}	
	}

	#normalized vectors for bias over each index per base
	A_vector <- A_vector / norm_count
	C_vector <- C_vector / norm_count
	G_vector <- G_vector / norm_count
	T_vector <- T_vector / norm_count
	Av <- c(rep("A", 6))
	Cv <- c(rep("C", 6))
	Gv <- c(rep("G", 6))
	Tv <- c(rep("T", 6))

	index_vector <- c(1,2,3,4,5,6)

	A_df <- data.frame(index = index_vector, bias = A_vector, base = Av)
	C_df <- data.frame(index = index_vector, bias = C_vector, base = Cv)
	G_df <- data.frame(index = index_vector, bias = G_vector, base = Gv)
	T_df <- data.frame(index = index_vector, bias = T_vector, base = Tv)


	final_df <- rbind(A_df, C_df, G_df, T_df)

	p <- ggplot(final_df, aes(index, bias))+ geom_line(aes(colour = base, group = base))
	p <- p + scale_x_continuous(breaks = 1:6)
	p <- p + xlab("hexamer index") + ylab("Ratio (bias weight)")

	p

}

