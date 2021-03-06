% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sleuth.R
\name{transcripts_from_gene}
\alias{transcripts_from_gene}
\title{Get the names of the transcripts associated to a gene}
\usage{
transcripts_from_gene(obj, test, test_type, which_model, gene_colname,
  gene_name)
}
\arguments{
\item{obj}{a \code{sleuth} object}

\item{test}{a character string denoting which beta to use}

\item{test_type}{either 'wt' for wald test or 'lrt' for likelihood ratio test}

\item{which_model}{a character string denoting which model to use}

\item{gene_colname}{the name of the column in which the desired gene apperas gene appears. Once genes have been added to a sleuth
object, you can inspect the genes names present in your sleuth object via \code{obj$target_mapping}, assuming 'obj' is the name of your sleuth object.
This parameter refers to the name of the column that the gene you are searching for appears in. Checkout the column names using \code{names(obj$target_mapping)}}

\item{gene_name}{a string containing the name of the gene you are interested in}
}
\value{
a vector of strings containing the names of the transcripts that map to a gene
}
\description{
Get the names of the transcripts associated to a gene, assuming genes are added to the input \code{sleuth} object.
}
