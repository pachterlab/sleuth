% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sleuth.R
\name{transform_fun_tpm<-}
\alias{transform_fun_tpm<-}
\title{Change sleuth transform TPM function}
\usage{
transform_fun_tpm(obj, fxn) <- value
}
\description{
Replace the transformation function of a sleuth object for TPM values
}
\details{
NOTE: if you change the transformation function after having done a fit,
the fit(s) will have to be redone using the new transformation.
}
\examples{
transform_fun_tpm(x) <- function(x) identity(x)
}
