% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alglib.R
\name{havealglib}
\alias{havealglib}
\title{Check whether chebpol has the ALGLIB library}
\usage{
havealglib()
}
\value{
Returns TRUE if ALGLIB is available. Otherwise FALSE.
}
\description{
If ALGLIB was available at compile time, it can be used for compact
support radial basis function interpolations on scattered data.
This function checks whether ALGLIB is available.
}
