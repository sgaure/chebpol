\name{lagappx}
\alias{lagappx}
\title{Lagrange interpolation on a grid}
\concept{Lagrange interpolation}
\description{
  A modified Lagrange interpolation on an arbitrary Cartesian product.
}

\usage{
  lagappx(val, dims=NULL, intervals=NULL, grid=NULL, ...)
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
\item{...}{Further arguments to the function, if \code{is.function(val)}.}
}


\value{
A \code{function(x)} defined on the hypercube, approximating the given
function.  The interpolant function uses the barycenteric Lagrangian interpolation.
}

\details{
  A call \code{fun <- lagappx(val,grid)} creates an 
  interpolant function on the grid.  The value on the grid points will be
  exact, the value between the grid points is interpolated with a polynomial.

  If \code{val} is a function it will be evaluated on the grid.

  If only \code{dims} is specified, Chebyshev sets are generated as
  \code{grid=chebknots(dims,intervals)}. This yields a Chebyshev
  approximation.
  Other grids can be tried out, but remember that the Lagrangian
  interpolant is a polynomial, so uniform grids may suffer from the Runge phenomenon.
  The modification from standard Lagrange is that the weights are
  constant, no matter which grid, we use the weights \code{c(0.5,1,...,1,0.5)} from a Chebyshev
  grid of the form \code{cos((0:n)*pi/n)}. The author do not know why
  this works so well even for uniform grids, but it seems to do.
}

\examples{
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
lg1 <- lagappx(exp,grid=list(s))
## test it
lg1(su) - exp(su)

## two dimensional approximation
f <- function(x) exp(sum(x^2))
grid <- list(s,su)

lg2 <- lagappx(evalongrid(f,grid=grid),grid=grid)
# an equivalent would be lg2 <- lagappx(f,grid)

a <- runif(2); lg2(a); f(a)
}