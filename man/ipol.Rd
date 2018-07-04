\name{ipol}
\alias{ipol}
\title{Interpolation wrapper}
\description{
  A wrapper for all the functions in the package
}

\usage{
ipol(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
     method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic','crbf'),
     ...) 
}

\arguments{
\item{val}{array or function. Function values on a grid, or the function
  itself. If it is the values, the \code{dim}-attribute
  must be appropriately set.}
\item{dims}{integer vector. The number of grid points in
  each dimension. Not needed if \code{grid} is used.}
\item{intervals}{list of length 2 numeric vectors. The lower and upper
  bound in each dimension. Only used if \code{dims} is specified.}
\item{grid}{list. Each element is a vector of ordered grid-points for a
  dimension.  These need not be Chebyshev-knots, nor evenly spaced.}
\item{knots}{matrix. Each column is a point in an M-dimensional space.}
\item{k}{numeric. Additional value, used with some methods.}
\item{method}{character. The interpolation method to use.}
\item{...}{Further arguments to the function, if \code{is.function(val)}.}
}


\value{
  A \code{function(x, threads)} defined on a hypercube, an interpolant for the given
  function. The argument \code{x} can be a matrix of column vectors
  which are evaluated in parallel in a number of threads.
  The function yields values for arguments outside the
  hypercube as well, though it will typically be a poor approximation.
}

\details{
 \code{ipol} is just a wrapper around the other functions in \pkg{chebpol}.
  Which arguments to specify depends on the method.

  The method \code{"chebyshev"} needs only the number of Chebyshev knots
  to use in each dimension, and possibly which intervals if different
  from [-1, 1]. The number of knots is specified in the \code{dims}
  argument, the intervals in the \code{intervals} argument. See
  \code{\link{chebappx}}. 

  The method \code{"uniform"} is similar to the \code{"chebyshev"}, but
  uniformly spaced knots are created. See \code{\link{ucappx}}.
  The argument \code{intervals} generally goes with \code{dims} when
  something else than standard intervals \code{[-1, 1]} are used.

  The method \code{"multilinear"}, \code{"fh"} (Floater-Hormann), and
  \code{"general"} needs the argument \code{grid}.  These are the
  methods which can use arbitrary Cartesian grids.  See
  \code{\link{mlappx}}, \code{\link{fhappx}}, and
  \code{\link{chebappxg}}.

  The method \code{"polyharmonic"} needs the arguments \code{knots} and \code{k},
  see \code{\link{polyh}}.

  The Floater-Hormann method (\code{"fh"}) needs the \code{k} argument,
  which is passed to the \code{d} argument of \code{\link{fhappx}}, the
  degree of the blending polynomials. It defaults to 4.

  The \code{"crbf"} is the multilayer compact radial basis function
  interpolation from ALGLIB. It is only available if ALGLIB was
  available at compile time. In this case \code{k} must be vector of
  length 3, the rbase, the number of layers and the non-linear smoothing
  parameter lambda.
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

## multi dimensional approximation
f <- function(x) exp(-sum(x^2))
# a 3-dimensional 10x10x10 grid, first and third coordinate are non-uniform
grid <- list(s, su, sort(1-s))

# make multilinear, Floater-Hormann, Chebyshev and polyharmonic spline.
ml2 <- ipol(f, grid=grid, method='multilin')
fh2 <- ipol(f, grid=grid, method='fh')
ch2 <- ipol(f, dims=c(10,10,10), intervals=list(0:1,0:1,0:1), method='cheb')
knots <- matrix(runif(3*1000),3)
ph2 <- ipol(f, knots=knots, k=2, method='poly')
# make 5 points in R3 to test them on
m <- matrix(runif(3*5),3)
rbind(true=apply(m,2,f), ml=ml2(m), fh=fh2(m), cheb=ch2(m), poly=ph2(m))
}
