# General grids, Floater-Hormann


#' Floater-Hormann interpolation on a grid
#' 
#' The Floater-Hormann interpolation.
#' 
#' A call \code{fun <- fhappx(val,grid)} creates a Floater-Hormann rational
#' interpolant function \code{fun} on the grid. The degree of the blending
#' polynomial,\code{d}, can be a vector, different for each grid-dimension.  In
#' theory, \code{d} can be any integer between 0 and the grid-dimension, a
#' higher \code{d} yields a smoother fit. However, increasing \code{d} leads to
#' exponential growth in rounding errors, so there is a tradeoff somewhere
#' which depends on the analyticity of the function, and not in an obvious way.
#' Current recommendations is to start low, at 3 or 4, and increase if
#' necessary.
#' 
#' If \code{val} is a function it will be evaluated on the grid.
#' 
#' @param val array or function. Function values on a grid, or the function
#' itself. If it is the values, the \code{dim}-attribute must be appropriately
#' set.
#' @param grid list.  Each element is a vector of ordered grid-points for a
#' dimension.  These need not be Chebyshev-knots, nor evenly spaced.
#' @param d integer. The degree of the blending polynomial.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x)} defined on the hypercube, approximating the
#' given function.  The interpolant function uses the barycentric
#' Floater-Hormann interpolation.
#' @examples
#' 
#' ## evenly spaced grid-points
#' su <- seq(0,1,length.out=10)
#' ## irregularly spaced grid-points
#' s <- su^3
#' ## create approximation on the irregularly spaced grid
#' fh1 <- fhappx(exp,grid=list(s))
#' ## test it
#' fh1(su) - exp(su)
#' 
#' ## two dimensional approximation
#' f <- function(x) exp(sum(x^2))
#' grid <- list(s,su)
#' 
#' fh2 <- fhappx(evalongrid(f,grid=grid),grid=grid)
#' # an equivalent would be fh2 <- fhappx(f,grid)
#' 
#' a <- runif(2); fh2(a); f(a)
#' 
#' @export
fhappx <- function(val,grid=NULL, d=1, ...) {
  x <- threads <- NULL; rm(x,threads) # avoid warning about undefined vars
  if(is.null(grid)) 
    stop('Must specify grid')
  if(!is.list(grid)) grid <- list(grid)
  grid <- lapply(grid,as.numeric)
  if(is.function(val)) val <- evalongrid(val, grid=grid, ...)
  dd <- as.integer(d)
  dd <- rep(dd, length(grid) %/% length(d))
  # calculate weights, formula 18 in Floater & Hormann, in C, parallelized
  weights <- .Call(C_FHweights, grid, dd, getOption('chebpol.threads'))
  local(vectorfun(.Call(C_FH,x,val,grid,weights,threads,NULL), 
                  args=alist(x=,threads=getOption('chebpol.threads')),
                  arity=length(grid),
                  domain=lapply(grid,range)),
        list(val=val,grid=grid,weights=weights))
}
