\name{chebpol-package}
\alias{chebpol-package}
\alias{chebpol}
\docType{package}
\title{Methods for creating multivariate Chebyshev interpolations on hypercubes}
\concept{Multivariate Chebyshev approximation}
\concept{Discrete Cosine Transform}
\description{
  The package contains methods for creating multivariate Chebyshev
  interpolations for real-valued functions on hypercubes.  Some methods for
  non-Chebyshev grids are also provided, multilinear and Floater-Hormann
  for Cartesian product grids, and a polyharmonic spline for scattered multi
  dimensional data.
}

\details{
  Given a real-valued function on a hypercube, or hyper-rectangle,
  it may be approximated by a (multivariate) Chebyshev polynomial.
  In the one-dimensional case, the Chebyshev approximation is
  constructed by evaluating the function in certain points, a
  Chebyshev grid, and fitting a polynomial to these values.
  Thus, one needs the function values on a set of
  prespecified points.
  The multivariate case is similar, the grid is the Cartesian product of
  one-dimensional grids.  I.e. all combinations of all grid-points.

  The Chebyshev coefficients for the interpolating polynomial in the
  one-dimensional case is a simple linear transform of the function
  values.  The Chebyshev-transform, or Discrete Cosine Transform, being
  a variant of the Fourier transform, factors over tensor products, thus
  the multivariate transform is just a tensor product of several
  one-dimensional transforms.  If \href{http://fftw.org}{\acronym{FFTW}} was
  available at compile time, \pkg{chebpol} uses it to generate
  the Chebyshev coefficients, otherwise a slower and more
  memory-demanding matrix method is used, and a warning message is issued when
  the package is attached.
  The Chebyshev-approximation is defined on the interval [-1,1], but it
  is straightforward to map any interval into [-1,1], thus making
  Chebyshev approximation on an interval of choice.  Or, a hypercube of choice.

  The primary method of the package is \code{\link{ipol}} which just
  dispatches to one of the following methods.
  For Chebyshev interpolation, there is \code{\link{chebappx}} which takes as
  input the function values on the grid, possibly together with a
  hypercube specification in the form of a list of intervals.  It
  produces a function which interpolates on the hypercube.
  There is also a wrapper called \code{\link{chebappxf}} which may be used if
  one has the function to be approximated rather than only its values in
  the grid-points.
  There is even an interpolation for uniform grids in
  \code{\link{ucappx}}, with a wrapper in \code{\link{ucappxf}}
  with some examples.  And a more general for arbitrary
  Cartesian-product grids in \code{\link{chebappxg}} with a wrapper in
  \code{\link{chebappxgf}}. These are based on transforms of Chebyshev-polynomials.
  A multilinear interpolation is available in \code{\link{mlappx}}. The
  Floater-Hormann rational interpolation is available in
  \code{\link{fhappx}}. For scattered multi dimensional data there is
  a polyharmonic spline interpolation in \code{\link{polyh}}.
  
  There are also functions for producing Chebyshev
  grids (\code{\link{chebknots}}) as well as a support function for evaluating a
  function on a grid (\code{\link{evalongrid}}), and a function for
  finding the Chebyshev coefficients (\code{\link{chebcoef}}).

  All the generated interpolant functions accept as an argument a matrix of column
  vectors. Except for the polyharmonic spline, the generated functions
  also accept an argument \code{threads=getOption('chebpol.threads')} to
  utilize more than one CPU if a matrix of column vectors is evaluated.
  The option \code{chebpol.threads} is initialized from the environment
  variable \code{CHEBPOL_THREADS} upon loading of the package.
  The polyharmonic spline will instead benefit from a multi-threaded
  BLAS, such as Intel's MKL.
}

\examples{
## make some function values on a 50x50x50 grid
dims <- c(x=50,y=50,z=50)
f <- function(x) sqrt(1+x[1])*exp(x[2])*sin(5*x[3])^2
value <- evalongrid(f , dims)
##fit a Chebyshev approximation to it. Note that the value-array contains the
##grid-structure. 
ch <- ipol(value,dims=dims,method='cheb')
## To see the full grid, use the chebknots function and expand.grid
\dontrun{
head(cbind(expand.grid(chebknots(dims)), value=as.numeric(value),
      appx=as.numeric(evalongrid(ch,dims))))
}
## Make a Floater-Hormann approximation on a uniform grid as well
fh <- ipol(f,grid=lapply(dims,function(n) seq(-1,1,length.out=n)),method='fh',k=5)
## evaluate in some random points in R3
m <- matrix(runif(15,-1,1),3)
data.frame(true=apply(m,2,f), cheb=ch(m), fh=fh(m))
}
\keyword{DCT}
\keyword{Discrete Cosine Transform}
\keyword{Chebyshev approximation}
