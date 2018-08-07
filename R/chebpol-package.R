#' Methods for creating multivariate interpolations on hypercubes
#' 
#' The package contains methods for creating multivariate
#' interpolations for real-valued functions on hypercubes. The
#' methods include classical Chebyshev interpolation, multilinear
#' and Floater-Hormann for arbitrary Cartesian product grids, and simplex linear
#' and polyharmonic spline for scattered multi dimensional data.
#' 
#' The primary method of the package is \code{\link{ipol}} which
#' dispatches to some other method.  All the generated
#' \link{interpolant}s accept as an argument a matrix of column
#' vectors. The generated functions also accept an argument
#' \code{threads=getOption('chebpol.threads')} to utilize more than
#' one CPU if a matrix of column vectors is evaluated.  The option
#' \code{chebpol.threads} is initialized from the environment variable
#' \code{CHEBPOL_THREADS} upon loading of the package. It defaults to \code{1}.
#'
#' The interpolants are ordinary R-objects and can be saved with \code{save()} and loaded
#' later with \code{load()} or serialized/unserialized with other tools, just like any R-object.
#' However, they contain calls to functions in the package, and while the author will make efforts
#' to ensure that generated interpolants are compatible with future versions of \pkg{chebpol},
#' I can issue no such absolute guarantee.
#'
#' @section Chebyshev:
#' If we are free to evaluate the function to interpolate in arbitrary points, we can use
#' a Chebyshev interpolation. The classical one is available with
#' \code{\link{ipol}(...,method='chebyshev')}. 
#' @section Uniform grids:
#' There are several options if your function must be evaluated in a uniform grid.
#' There is the Floater-Hormann rational interpolation available with \code{\link{ipol}(...,method='fh')}.
#' There is a transformed Chebyshev variant \code{\link{ipol}(..., method='uniform')}.
#' @section Arbitrary grids:
#' For grids which are not uniform, but still Cartesian products of one-dimensional grids,
#' there is the Floater-Hormann interpolation \code{\link{ipol}(...,method='fh')}, and a transformed
#' Chebyshev variant \code{\link{ipol}(...,method='general')}, as well as a multilinear
#' \code{\link{ipol}(...,method='multilinear')}. These methods work on uniform grids as well.
#' @section Scattered data:
#' For scattered data, not necessarily organised as a Cartesian product grid, there is
#' a simplex linear interpolation available with \code{\link{ipol}(...,method='simplexlinear')},
#' and a polyharmonic spline with \code{\link{ipol}(...,method='polyharmonic')}.
#' @section Support functions:
#' There are also functions for producing Chebyshev grids
#' (\code{\link{chebknots}}) as well as a function for evaluating a
#' function on a grid (\code{\link{evalongrid}}), and a function for finding
#' the Chebyshev coefficients (\code{\link{chebcoef}}).
#' 
#' @name chebpol-package
#' @aliases chebpol-package chebpol
#' @docType package
#' @concept DCT Floater-Hormann
#' @keywords DCT Chebyshev interpolation Floater-Hormann
#' @seealso \link{ipol}, \link{interpolant}
#' @examples
#' 
#' ## make some function values on a 50x50x50 grid
#' dims <- c(x=50,y=50,z=50)
#' f <- function(x) sqrt(1+x[1])*exp(x[2])*sin(5*x[3])^2
#' value <- evalongrid(f , dims)
#' ##fit a Chebyshev approximation to it. Note that the value-array contains the
#' ##grid-structure. 
#' ch <- ipol(value,method='cheb')
#' ## To see the full grid, use the chebknots function and expand.grid
#' \dontrun{
#' head(cbind(expand.grid(chebknots(dims)), value=as.numeric(value),
#'       appx=as.numeric(evalongrid(ch,dims))))
#' }
#' ## Make a Floater-Hormann approximation on a uniform grid as well
#' fh <- ipol(f,grid=lapply(dims,function(n) seq(-1,1,length.out=n)),method='fh',k=5)
#' ## evaluate in some random points in R3
#' m <- matrix(runif(15,-1,1),3)
#' rbind(true=apply(m,2,f), cheb=ch(m), fh=fh(m))
#' 
#' @useDynLib chebpol, .registration=TRUE, .fixes='C_'
NULL
