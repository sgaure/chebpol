# we can actually find the grid-maps for uniform grids.
# the Chebyshev knots are cos(pi*(j+0.5)/n) for j=0..n-1 These should
# map into the n grid points. These have distance 2/(n-1), starting in -1, ending in 1
# so they are -1 + 2*j/(n-1). After some manipulation, the function is:

ugm <- function(x,n) sin(0.5*pi*x*(1-n)/n)



#' Interpolation on a uniform grid
#' 
#' A poor-man's approximation on uniform grids.  If you for some reason can't
#' evaluate your function on a Chebyshev-grid, but instead have a uniform grid,
#' you may use this function to create an interpolation.
#' 
#' This does about the same as \code{\link{chebappxg}} for unform grids, though
#' no grid map function is constructed, as a fixed such function is used.
#' 
#' A Chebyshev-interpolation \code{ch} is made for \code{val} with
#' \code{\link{chebappx}}. Upon evaluation the uniform grid in each dimension
#' is mapped differentiably to the Chebyshev-knots so that \code{ch} is
#' evaluated in \eqn{sin(\frac{\pi }{sin(0.5*pi*x*(1-n)/n)}\eqn{
#' x(1-n)}{2n})}{sin(0.5*pi*x*(1-n)/n)} where \code{n} is the number of knots
#' in the dimension, possibly after \code{x} has been remapped from the
#' hypercube interval to [-1,1].
#' 
#' Thus, the interpolation is not a polynomial.
#' 
#' For \code{ucappx} the function values are provided, the number of grid
#' points in each dimension is to be found in \code{dim(val)}. For
#' \code{ucappxf} the function to be interpolated is \code{fun}, and the number
#' of grid points is passed in \code{dims}.
#' 
#' As the example shows, this approximation is better than the Chebyshev
#' approximation for some functions.
#' 
#' @param val Array. Function values on a grid.
#' @param intervals List of vectors of length two. Specifying the hypercube
#' extent in each dimension
#' @return A \code{function(x)} defined on the hypercube, approximating the
#' given function.
#' @examples
#' 
#' \dontrun{
#' # Runge function
#' f <- function(x) 1/(1+25*x^2)
#' grid <- seq(-1,1,length.out=15)
#' val <- f(grid)
#' uc <- Vectorize(ucappx(val))
#' # and the Chebyshev
#' ch <- Vectorize(chebappxf(f,15))
#' # test it at 10 random points
#' t(replicate(10,{a<-runif(1,-1,1); c(arg=a, uc=uc(a), true=f(a), cheb=ch(a))}))
#' }
#' @export
#' @keywords internal
ucappx <- function(...) deprecated('ucappx',...)
ucappx.real <- function(val, intervals=NULL) {
  if(is.null(dim(val))) dim(val) <- length(val)
  dims <- dim(val)
  ch <- chebappx.real(val)
  if(is.null(intervals)) {
    gridmap <- compiler::cmpfun(function(x) mapply(function(xi,d) ugm(xi,d),x,dims))
  } else {
    # precompute interval mid points and inverse lengths
    md <- lapply(intervals,mean)
    ispan <- lapply(intervals, function(i) 2/diff(i))
    gridmap <- compiler::cmpfun(function(x) mapply(function(xi,mid,is,d) ugm(is*(xi-mid),d),x,md,ispan,dims))
  }
  gm <- compiler::cmpfun(function(x) {
    if(is.matrix(x)) apply(x,2,gridmap) else gridmap(x)
  })
  vectorfun(function(x, threads=getOption('chebpol.threads')) ch(gm(x), threads), 
            arity=length(dims), 
            domain=if(is.null(intervals)) lapply(rep(-1,length(dim(val))),c,1) else intervals)

}

#' @rdname ucappx
#' @param fun Function to be interpolated.
#' @param dims Integer. Number of grid points in each dimension.
#' @param ... Further arguments to \code{fun}.
#' @export
#' @keywords internal
ucappxf <- function(...) deprecated('ucappxf',...)
ucappxf.real <- function(fun, dims, intervals=NULL,...) {
  if(is.null(intervals))
    return(ucappx.real(evalongrid(fun,...,grid=lapply(dims,function(d) seq(-1,1,length.out=d)))))
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  return(
    ucappx.real(evalongrid(fun,...,
                      grid=mapply(function(d,i) seq(min(i),max(i),length.out=d),
                                  dims, intervals,SIMPLIFY=FALSE)),
           intervals))
}
