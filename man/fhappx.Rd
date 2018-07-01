\name{fhappx}
\alias{fhappx}
\title{Floater-Hormann interpolation on a grid}
\concept{Floater-Hormann interpolation}
\description{
  The Floater-Hormann interpolation.
}

\usage{
  fhappx(val, grid=NULL, d=1, ...)
}

\arguments{
\item{val}{Array or function. Function values on a grid, or the function
  itself. If it is the values, the \code{dim}-attribute
  must be appropriately set.}
\item{grid}{A list.  Each element is a vector of ordered grid-points for a
  dimension.  These need not be Chebyshev-knots, nor evenly spaced.}
\item{d}{integer. The degree of the blending polynomial.}
\item{...}{Further arguments to the function, if \code{is.function(val)}.}
}


\value{
A \code{function(x)} defined on the hypercube, approximating the given
function.  The interpolant function uses the barycenteric
Floater-Hormann interpolation.
}

\details{
  A call \code{fun <- fhappx(val,grid)} creates a Floater-Hormann
  rational interpolant function on the grid. The degree of the blending
  polynomial,\code{d}, can be a vector, different for each grid-dimension.

  If \code{val} is a function it will be evaluated on the grid.
}

\examples{
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
fh1 <- fhappx(exp,grid=list(s))
## test it
fh1(su) - exp(su)

## two dimensional approximation
f <- function(x) exp(sum(x^2))
grid <- list(s,su)

fh2 <- fhappx(evalongrid(f,grid=grid),grid=grid)
# an equivalent would be fh2 <- fhappx(f,grid)

a <- runif(2); fh2(a); f(a)
}
