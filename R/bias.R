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
	
	#variable declaration
	A1 <- A2 <- A3 <- A4 <- A5 <- A6 <- C1 <- C2 <- C3 <- C4 <- C5 <- C6 <- NULL
	G1 <- G2 <- G3 <- G4 <- G5 <- G6 <- T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- NULL
	#set up
	for (i in 1:nrow(bt)) {
		iter <- strsplit(bt$hexamer[i], "")[[1]]
		bias_append <- bt$bias_weights[i]
		#index 1-6
		for (i1 in iter[1]) {
			if (i1 == "A") {
				A1 <- A1 + bias_append
			}
			if (i1 == "C") {
				C1 <- C1 + bias_append
			}
			if (i1 == "G") {
				G1 <- G1 + bias_append
			}
			if (i1 == "T") {
				T1 <- T1 + bias_append
			}
		}
		for (i2 in iter[2]) {
			if (i2 == "A") {
				A2 <- A2 + bias_append
			}
			if (i2 == "C") {
				C2 <- C2 + bias_append
			}
			if (i2 == "G") {
				G2 <- G2 + bias_append
			}
			if (i2 == "T") {
				T2 <- T2 + bias_append
			}
		}
		for (i3 in iter[3]) {
			if (i3 == "A") {
				A3 <- A3 + bias_append
			}
			if (i3 == "C") {
				C3 <- C3 + bias_append
			}
			if (i3 == "G") {
				G3 <- G3 + bias_append
			}
			if (i3 == "T") {
				T3 <- T3 + bias_append
			}
		}
		for (i4 in iter[4]) {
			if (i4 == "A") {
				A4 <- A4 + bias_append
			}
			if (i4 == "C") {
				C4 <- C4 + bias_append
			}
			if (i4 == "G") {
				G4  <- G4 + bias_append
			}
			if (i4 == "T") {
				T4 <- T4 + bias_append
			}
		}
		for (i5 in iter[5]) {
			if (i5 == "A") {
				A5 <- A5 + bias_append
			}
			if (i5 == "C") {
				C5 <- C5 + bias_append
			}
			if (i5 == "G") {
				G5 <- G5 + bias_append
			}
			if (i5 == "T") {
				T5 <- T5 + bias_append
			}
		}
		for (i6 in iter[6]) {
			if (i6 == "A") {
				A6 <- A6 + bias_append
			}
			if (i6 == "C") {
				C6 <- C6 + bias_append
			}
			if (i6 == "G") {
				C6 <- C6 + bias_append
			}
			if (i6 == "T") {
				T6 <- T6 + bias_append
			}
		}

	}


	#normalized vectors for bias over each index per base
	A_vector <- c(A1, A2, A3, A4, A5, A6) / norm_count
	C_vector <- c(C1, C2, C3, C4, C5, C6) / norm_count
	G_vector <- c(G1, G2, G3, G4, G5, G6) / norm_count
	T_vector <- c(T1, T2, T3, T4, T5, T6) / norm_count

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

