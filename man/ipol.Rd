\name{ipol}
\alias{ipol}
\title{Interpolation wrapper}
\description{
  A wrapper for all the functions in the package
}

\usage{
ipol(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
           method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic'),
           ...) 
}

\arguments{
\item{val}{Array or function. Function values on a grid, or the function
  itself. If it is the values, the \code{dim}-attribute
  must be appropriately set.}
\item{dims}{integer vector. Integer vector, the number of grid points in
  each dimension. Not needed if \code{grid} is used.}
\item{intervals}{list of length 2 numeric vectors. The lower and upper
  bound in each dimension. Only used if \code{dims} is specified.}
\item{grid}{A list.  Each element is a vector of ordered grid-points for a
  dimension.  These need not be Chebyshev-knots, nor evenly spaced.}
\item{knots}{matrix. Each column is a point in an M-dimensional space.}
\item{k}{numeric. Additional value, used with some methods.}
\item{method}{character. The interpolation method to use.}
\item{...}{Further arguments to the function, if \code{is.function(val)}.}
}


\value{
A \code{function(x, ...)} defined on hypercube, an interpolant for the given
function.  The function yields values for arguments outside the
hypercube as well, though it will typically be a poor approximation.
}

\details{
 \code{ipol} is just a wrapper around the other functions in \pkg{chebpol}.
  Which arguments to specify depends on the method.

  The argument \code{dims} is used by the methods \code{"chebyshev"} and
  and \code{"uniform"}. These are the methods where the
  interpolation points can be inferred by how many there are. See
  \code{\link{chebappx}} and \code{\link{ucappx}}.
  The argument \code{intervals} generally goes with \code{dims} when
  something else than standard intervals \code{[-1, 1]} are used.

  The argument \code{grid} is used with methods \code{"multilinear"},
  \code{"fh"}, and \code{"general"}. These are the methods which can use
  arbitrary Cartesian grids.  See \code{\link{mlappx}},
  \code{\link{fhappx}}, and \code{\link{chebappxg}}.

  The arguments \code{knots} and \code{k} are used with method
  \code{"polyharmonic"}, see \code{\link{polyh}}. \code{k} with method
  \code{"fh"} is used as the \code{d} argument to \code{\link{fhappx}},
  the degree of the blending polynomials, and defaults to 4.
}

\examples{
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ml1 <- ipol(exp,grid=list(s), method='multilin')
fh1 <- ipol(exp, grid=list(s), method='fh')
## test it, since exp is convex, the linear approximation lies above
## the exp between the grid points
ml1(su) - exp(su)
fh1(su) - exp(su)

## multi linear approx
f <- function(x) exp(sum(x^2))
grid <- list(s,su)

ml2 <- ipol(f, grid=grid, method='multilin')
fh2 <- ipol(f, grid=grid, method='fh')
ch2 <- ipol(f, dims=c(10,10), intervals=list(0:1,0:1), method='cheb')
knots <- matrix(runif(2*100),2)
ph2 <- ipol(f, knots=knots, k=2, method='poly')
a <- runif(2)
f(a); ml2(a); fh2(a); ch2(a); ph2(a)

}
