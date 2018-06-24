\name{ucappx}
\alias{ucappx}
\alias{ucappxf}
\title{Interpolation on a uniform grid}
\description{
  A poor-man's approximation on uniform grids.
  If you for some reason can't evaluate your function on a
  Chebyshev-grid, but instead have a uniform grid, you may use this function
  to create an interpolation.
}

\usage{
ucappx(val, intervals=NULL)
ucappxf(fun, dims, intervals=NULL, ...)
}

\arguments{
\item{val}{Array. Function values on a grid.}
\item{intervals}{List of vectors of length two. Specifying the hypercube
extent in each dimension}
\item{fun}{Function to be interpolated.}
\item{dims}{Integer. Number of grid points in each dimension.}
\item{...}{Further arguments to \code{fun}.}

}


\value{
A \code{function(x)} defined on the hypercube, approximating the given
function.
}

\details{
  This does about the same as \code{\link{chebappxg}} for unform grids,
  though no grid map function is constructed, as a fixed such function is
  used.

  A Chebyshev-interpolation \code{ch} is made for \code{val} with
  \code{\link{chebappx}}. Upon evaluation the uniform grid in each
  dimension is mapped differentiably to the Chebyshev-knots so that
  \code{ch} is evaluated in \eqn{sin(\frac{\pi
  x(1-n)}{2n})}{sin(0.5*pi*x*(1-n)/n)} where
  \code{n} is the number of knots in the dimension, possibly after
  \code{x} has been remapped from the hypercube interval to [-1,1].
  
  Thus, the interpolation is not a polynomial.

  For \code{ucappx} the function values are provided, the number of grid
  points in each dimension is to be found in \code{dim(val)}. For
  \code{ucappxf} the function to be interpolated is \code{fun}, and the
  number of grid points is passed in \code{dims}.
  
  As the example shows, this approximation is better than the Chebyshev
  approximation for some functions.
}

\examples{
# Runge function
f <- function(x) 1/(1+25*x^2)
grid <- seq(-1,1,length.out=15)
val <- f(grid)
uc <- Vectorize(ucappx(val))
# and the Chebyshev
ch <- Vectorize(chebappxf(f,15))
# test it at 10 random points
t(replicate(10,{a<-runif(1,-1,1); c(arg=a, uc=uc(a), true=f(a), cheb=ch(a))}))
}