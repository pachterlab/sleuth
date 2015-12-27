#make sure this is in the current directory of samples
#make sure bias_table is already loaded in

#test for local machine

#setwd("~/Documents/Research/sleuth/ellahi/results")

bias_graph <- function(obj,
	sample = "",
	...) {
	stopifnot( is(obj, 'sleuth') )

	so <- obj
	bt <- bias_table(so, sample) 
	totalSum <- sum(bt$bias_weights)

	#will return ggplot object of a line graph
	p <- ggplot(data = TEMPORARY, 
		aes(x = offset, 
			y = ratio,
			group = bases,
			color = bases)) 
	p <- p + geom_line() + geom_point()
	p <- p + xlab("Offset from 5' fragment end") + ylab("Ratio (bias weight)")

	p

}

