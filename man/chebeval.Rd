\name{chebeval}
\alias{chebeval}
\title{Evaluate a Chebyshev interpolation in a point}
\description{
  Given Chebyshev coefficients, evaluate the interpolation in a point.
}
\usage{
chebeval(x, coef, intervals=NULL)
}

\arguments{
\item{x}{The point to evaluate.}
\item{coef}{The Chebyshev coefficients. Typically from a call to
  \code{\link{chebcoef}}, possibly modified.}
\item{intervals}{A list of minimum and maximum values. One for each
  dimension of the hypercube.}
}
\value{
A numeric. The interpolated value.
}

\examples{
# make a function which is known to be unsuitable for Chebyshev approximation
f <- function(x) sign(x)
# make a standard Chebyshev interpolation
ch <- Vectorize(chebappxf(f,50))
# then do a truncated interpolation
val <- evalongrid(f,50)
coef <- chebcoef(val)
# truncate the high frequencies
coef[-(1:10)] <- 0
# make a truncated approximation
tch <- Vectorize(function(x) chebeval(x,coef))
# make a lower degree also
ch2 <- Vectorize(chebappxf(f,10))
# plot the functions
\dontrun{
s <- seq(-1,1,length.out=400)
plot(s,ch(s),col='red',type='l')
lines(s,tch(s),col='blue')
lines(s,f(s))
lines(s,ch2(s),col='green')
}
}