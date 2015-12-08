bias_graph <- function(obj,
	sample_input = "",
	...) {

	stopifnot( is(obj, 'sleuth') )
	so <- obj
	bt <- bias_table(so, sample = sample_input)

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

