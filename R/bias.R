#make sure this is in the current directory of samples
#make sure bias_table is already loaded in

#test for local machine

#setwd("~/Documents/Research/sleuth/ellahi/results")


#assumptions: 
#1. there are an equal number of hexamers for each base
#2. they are in the order of ACGT
#PRETTY BRAINDEAD, THERE HAS TO BE A BETTER WAY TO DO THIS

bias_graph <- function(obj,  
	sample = "",
	...) {
	stopifnot( is(obj, 'sleuth'))
  
	#sample df construction
	so <- obj
	bt <- bias_table(so, sample) 
	bt <- data.frame(hexamer = bt$hexamer, bias_weights = bt$bias_weights)
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

	#df construction of actual frame
	index_df <- data.frame(A = A_vector, C = C_vector, G = G_vector, T = T_vector)

	#will return ggplot object of a line graph
	p <- ggplot(data = TEMPORARY, 
		aes(x = index, 
			y = ratio,
			group = bases,
			color = bases)) 
	p <- p + geom_line() + geom_point()
	p <- p + xlab("hexamer index") + ylab("Ratio (bias weight)")

	p

}

