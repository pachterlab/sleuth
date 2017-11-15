

#' Lancaster method
#'
#' Weighted p-value aggregation. Code is taken from github.com/lynnyi/aggregation.
#' @param pvalues A vector of p-values (i.e. between 0 and 1). NAs will be filtered out.
#' @param weights A vector of weights, each associated with its respective p-value. Weights must be nonegative. NAs and negative weights will be filtered out with corresponding p-values.
#' @examples
#' lancaster(c(.1, .5), c(2, 4))
lancaster <- function(pvalues, weights)
{
	if(length(weights) != length(pvalues))
	{
		stop('Length of weights not equal to length of pvalues')
	}
	weights <- weights[!is.na(pvalues)]
	pvalues <- pvalues[!is.na(pvalues)]
	pvalues <- pvalues[weights>0]
	weights <- weights[weights>0]
	if(length(pvalues) == 0)
	{
		return(NA)
	}
	if(length(pvalues) == 1)
	{
		return(pvalues)
	}
	if(any(pvalues < 10e-320))
	{
		warning('Extreme p-values around and below 10e-320 will produce a p-value of 0. Replace extreme p-values with 10e-320 to obtain an upper bound for the aggregated p-value.')
	}
	t <- sapply(1:length(pvalues), function(i) lts(pvalues[i], weights[i]))
	t <- sum(t)
	p <- pchisq(t, sum(weights), lower.tail=FALSE)
	p
}

lts <- function(pvalue, weight)
{
	qgamma(pvalue, shape = weight /2, scale = 2, lower.tail=FALSE)
}
