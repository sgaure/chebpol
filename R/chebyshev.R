# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots1 <- function(n, interval=NULL) {
  kn <- cos(pi*((1:n)-0.5)/n)
  if((n %% 2) == 1) kn[[(n+1)/2]] = 0
  if(is.null(interval)) return(kn)
  kn*diff(interval)/2 + mean(interval)
}



#' Create a Chebyshev-grid
#' 
#' Create a Chebyshev grid on a hypercube.
#' 
#' If \code{intervals} is not provided, it is assumed that the domain of the
#' function in each dimension is [-1,1].  Thus, standard Chebyshev knots are
#' produced.  If \code{dims} is of length 1, \code{intervals} may be a vector
#' of length 2 rather than a list with a vector of length 2.
#' 
#' @param dims The number of grid-points in each dimension.  For
#' Chebyshev-polynomial of degree \code{dims-1}.
#' @param intervals A list of vectors of length 2.  The lower and upper bounds
#' of the hypercube.
#' @return A array of dimension \code{dims}.  The Chebyshev grid-points.
#' @examples
#' 
#' ## Standard knots for degree 3
#' chebknots(4)
#' ## Knots in the interval [2,3] for degree 3
#' chebknots(4,interval=c(2,3))
#' ## Multivariate knots
#' chebknots(c(x=3,y=4,z=3))
#' ## Multivariate grid
#' \dontrun{
#' expand.grid(chebknots(c(x=3,y=4,z=5), list(c(1,3), c(4,6), c(800,900))))
#' }
#' 
#' @export
chebknots <- function(dims, intervals=NULL) {
  if(is.null(intervals)) {
    res <- lapply(dims,chebknots1)
  } else {
    if(length(dims) == 1 && is.numeric(intervals)) intervals=list(intervals)
    if(!is.list(intervals)) stop("intervals should be a list")
    res <- mapply(chebknots1,dims,intervals,SIMPLIFY=FALSE)
  }

  res
}

#' Compute Chebyshev-coefficients given values on a Chebyshev grid
#' 
#' Compute the multivariate Chebyshev-coefficients, given values on a Chebyshev
#' grid.
#' 
#' If \code{val} has no \code{dim}-attribute, it is assumed to be
#' one-dimensional of length the length of \code{val}.
#' 
#' If \pkg{chebpol} was compiled without \acronym{FFTW}, running
#' \code{chebcoef} on large grids may be slow and memory-demanding.
#' 
#' @param val An \code{array} of function values on a Chebyshev grid.  The
#' \code{dim}-attribute must be appropriately set.  If not set, it is assumed
#' to be one-dimensional.
#' @param dct Logical. Since the Chebyshev coefficients are closely related to
#' the DCT-II transform of \code{val}, the non-normalized real-even DCT-II
#' coefficients may be retrieved instead. I.e. those from FFTW_REDFT10 in each
#' dimension.  This is not used anywhere in the package, it is merely provided
#' as a convenience for those who might need it.
#' @return An array of Chebyshev-coefficients for an interpolating
#' Chebyshev-polynomial.
#' @seealso \code{\link{havefftw}}
#' @examples
#' 
#' ## Coefficients for a 2x3x4 grid
#' a <- array(rnorm(24),dim=c(2,3,4))
#' chebcoef(a)
#' 
#' @export
chebcoef <- function(val, dct=FALSE) {
  structure(.Call(C_chebcoef,as.array(val),dct,getOption('chebpol.threads')),dimnames=dimnames(val))
}



#' Evaluate a Chebyshev interpolation in a point
#' 
#' Given Chebyshev coefficients, evaluate the interpolation in a point.
#' 
#' 
#' @param x The point to evaluate.
#' @param coef The Chebyshev coefficients. Typically from a call to
#' \code{\link{chebcoef}}, possibly modified.
#' @param intervals A list of minimum and maximum values. One for each
#' dimension of the hypercube.
#' @param threads And integer. In case \code{x} is a matrix of column vectors,
#' use this number of threads in parallel to evaluate.
#' @return A numeric. The interpolated value.
#' @examples
#' 
#' # make a function which is known to be unsuitable for Chebyshev approximation
#' f <- function(x) sign(x)
#' # make a standard Chebyshev interpolation
#' ch <- ipol(f,dims=50,method='chebyshev')
#' # then do a truncated interpolation
#' val <- evalongrid(f,50)
#' coef <- chebcoef(val)
#' # truncate the high frequencies
#' coef[-(1:10)] <- 0
#' # make a truncated approximation
#' tch <- Vectorize(function(x) chebeval(x,coef))
#' # make a lower degree also
#' ch2 <- ipol(f,dims=10,method='chebyshev')
#' # plot the functions
#' \dontrun{
#' s <- seq(-1,1,length.out=400)
#' plot(s,ch(s),col='red',type='l')
#' lines(s,tch(s),col='blue')
#' lines(s,f(s))
#' lines(s,ch2(s),col='green')
#' }
#' 
#' @export
chebeval <- function(x,coef,intervals=NULL,threads=getOption('chebpol.threads')) {
  if(is.null(intervals)) return(.Call(C_evalcheb,coef,x,threads,NULL))
  # map into intervals
  .Call(C_evalcheb,coef,mapply(function(x,i) 2*(x[[1]]-mean(i))/diff(i),x,intervals),threads,NULL)
}

# return a function which is a Chebyshev interpolation


#' Chebyshev interpolation on a hypercube
#' 
#' Given function, or function values on a Chebyshev grid, create an
#' interpolatin function defined in the whole hypercube.
#' 
#' If \code{intervals} is not provided, it is assumed that the domain of the
#' function is the Cartesian product [-1,1] x [-1,1] x ... x [-1,1].  Where the
#' number of grid-points are given by \code{dim(val)}.
#' 
#' For \code{chebappxf}, the function is provided instead, and the number of
#' grid points in each dimension is in the vector \code{dims}. The function is
#' evaluated on the Chebyshev grid.
#' 
#' If \code{intervals} is provided, it should be a \code{list} with elements of
#' length 2, providing minimum and maximum for each dimension. Arguments to the
#' function will be transformed from these intervals into [-1,1] intervals.
#' 
#' The approximation function may be evaluated outside the hypercube, but be
#' aware that it may be highly erratic there, especially if of high degree.
#' 
#' @param val The function values on the Chebyshev grid. \code{val} should be
#' an array with appropriate dimension attribute.
#' @param intervals A list of minimum and maximum values. One for each
#' dimension of the hypercube. If NULL, assume [-1,1] in each dimension.
#' @return A function defined on the hypercube. A Chebyshev approximation to
#' the function \code{fun}, or the values provided in \code{val}.
#' @examples
#' \dontrun{
#' 
#' f <- function(x) exp(-sum(x^2))
#' ## we want 3 dimensions, i.e. something like
#' ## f(x,y,z) = exp(-(x^2 + y^2 + z^2))
#' ## 8 points in each dimension
#' gridsize <- list(8,8,8)
#' # get the function values on the Chebyshev grid
#' values <- evalongrid(f,gridsize)
#' # make an approximation
#' ch <- chebappx(values)
#' ## test it:
#' a <- runif(3,-1,1);ch(a)-f(a)
#' 
#' ## then one with domain [0.1,0.3] x [-1,-0.5] x [0.5,2]
#' intervals <- list(c(0.1,0.3),c(-1,-0.5),c(0.5,2))
#' # evaluate on the grid
#' values <- evalongrid(f,gridsize,intervals)
#' # make an approximation
#' ch2 <- chebappx(values,intervals)
#' a <- c(0.25,-0.68,1.43); ch2(a)-f(a)
#' # outside of domain:
#' a <- runif(3) ; ch2(a); f(a)
#' 
#' # Make a function on [0,2] x [0,1]
#' f <- function(y) uniroot(function(x) x-y[[1]]*cos(pi*x^2),lower=0,upper=1)$root*sum(y^2)
#' # approximate it
#' ch <- chebappxf(f,c(12,12),intervals=list(c(0,2),c(0,1)))
#' # test it:
#' a <- c(runif(1,0,2),runif(1,0,1)); ch(a); f(a)
#' 
#' # Lambert's W:
#' f <- function(y) uniroot(function(x) y - x*exp(x), lower=-1,upper=3)$root
#' W <- chebappxf(f,100,c(-exp(-1),3*exp(3)))
#' W(10*pi)*exp(W(10*pi))/pi
#' }
#' @export 
#' @keywords internal
chebappx <- function(...) deprecated('chebappx',...)
chebappx.real <- function(val,intervals=NULL) {
  if(is.null(dim(val))) {
    # allow for one-dimensional
    dim(val) <- length(val)
  }
  K <- length(dim(val))
   # allow for vector, e.g. intervals=c(0,1), put it inside list
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  cf <- chebcoef(val)

  if(is.null(intervals)) {
    # it's [-1,1] intervals, so drop transformation
    fun <- vectorfun(function(x,threads=getOption('chebpol.threads')) .Call(C_evalcheb,cf,x,threads,NULL),
                     arity=K,
                     domain=lapply(rep(-1,K),c,1))
    rm(val)
  } else {
    # it's intervals, create mapping into [-1,1]
    if(!is.list(intervals)) stop("intervals should be a list")
    if(any(sapply(intervals,length) != 2)) stop("interval elements should have length 2")
    if(length(intervals) != length(dim(val))) stop("values should have the same dimension as intervals ",
               length(intervals),' ',length(dim(val)))
    ispan <- sapply(intervals,function(x) 2/diff(x))
    mid <- sapply(intervals,function(x) mean(x))
    imap <- function(x) (x-mid)*ispan

    fun <- vectorfun(function(x,threads=getOption('chebpol.threads')) .Call(C_evalcheb,cf,imap(x), threads, NULL),
                     arity=K,
                     domain=intervals)
    rm(val)
  }
  fun
}

#' @rdname chebappx
#' @param fun The function to be approximated.
#' @param dims Integer. The number of Chebyshev points in each dimension.
#' @param ... Further arguments to \code{fun}.
#' @export
#' @keywords internal
chebappxf <- function(...) deprecated('chebappxf',...)

chebappxf.real <- function(fun,dims,intervals=NULL,...) {
  chebappx.real(evalongrid(fun,dims,intervals,...),intervals)
}

# interpolate on a non-Chebyshev grid. This is useful if you for some reason
# do not have the function values on a Chebyshev-grid, but on some other grid.  It comes at a cost,
# The interpolation may not be very good compared to the Chebyshev-one.

# val are the function values on a grid, an array of appropriate dimension
# in the order of expand.grid()

# if grid is unspecified, it is assumed that it is on a Chebyshev grid in [-1,1]
# If grid is specified, it is a list of vectors. The length of the list
# is the dimension of the grid. Each vector contains grid-points in increasing or decreasing order.
# val is assumed to be the function values on expand.grid(grid)
# The user-grid is mapped to Chebshev knots in [-1,1] by splinefun in pkg stats





#' Interpolation on a non-Chebyshev grid
#' 
#' A poor-man's approximation on non-Chebyshev grids.  If you for some reason
#' can't evaluate your function on a Chebyshev-grid, but instead have some
#' other grid which still is a Cartesian product of one-dimensional grids, you
#' may use this function to create an interpolation.
#' 
#' 
#' A call \code{fun <- chebappxg(val,grid)} does the following.  A Chebyshev
#' interpolation \code{ch} for \code{val} is created, on the [-1,1] hypercube.
#' For each dimension a grid-map function \code{gm} is created which maps the
#' grid-points monotonically into Chebyshev knots. For this, the function
#' \code{\link{splinefun}} with \code{method='hyman'} is used.  When
#' \code{fun(x)} is called, it translates to \code{ch(gm(x))}.  For uniform
#' grids, the function \code{\link{ucappx}} will produce a faster interpolation
#' in that a closed form \code{gm} is used.
#' 
#' \code{chebappxgf} is used if the function, rather than its values, is
#' available. The function will be evaluated on the grid.
#' 
#' Even though this approach works in simple cases it is not a panacea.  The
#' grid in each dimension should probably not be too irregularly spaced. I.e.
#' short and long gaps interspersed is likely to cause problems.
#' 
#' @param val Array. Function values on a grid.
#' @param grid A list. Each element is a sorted vector of grid-points for a
#' dimension. These need not be Chebyshev-knots, nor evenly spaced.
#' @param ... Further arguments to \code{fun}.
#' @return A \code{function(x)} defined on the hypercube, approximating the
#' given function.
#' @examples
#' \dontrun{
#' ## evenly spaced grid-points
#' su <- seq(0,1,length.out=10)
#' ## irregularly spaced grid-points
#' s <- su^3
#' ## create approximation on the irregularly spaced grid
#' ch <- Vectorize(chebappxg(exp(s),list(s)))
#' ## test it:
#' ch(su) - exp(su)
#' # try one with three variables
#' f <- function(x) exp(-sum(x^2))
#' grid <- list(s,su,su^2)
#' ch2 <- chebappxg(evalongrid(f,grid=grid),grid)
#' # test it at 10 random points
#' replicate(10,{a<-runif(3); ch2(a)-f(a)})
#' 
#' # Try Runge's function on a uniformly spaced grid.
#' # Ordinary polynomial fitting of high degree of Runge's function on a uniform grid
#' # creates large oscillations near the end of the interval. Not so with chebappxgf
#' f <- function(x) 1/(1+25*x^2)
#' chg <- Vectorize(chebappxgf(f,seq(-1,1,length.out=15)))
#' # also compare with Chebyshev interpolation
#' ch <- Vectorize(chebappxf(f,15))
#' \dontrun{
#'  # plot it
#'  s <- seq(-1,1,length.out=200)
#'  plot(s, f(s), type='l', col='black')
#'  lines(s, chg(s), col='blue')
#'  lines(s, ch(s), col='red')
#'  legend('topright',
#'         legend=c('Runge function','chebappxg on uniform grid','Chebyshev'),
#'         col=c('black','blue','red'), lty=1)
#' }
#' }
#' 
#' @export 
#' @keywords internal
chebappxg <- function(...) deprecated('chebappxg',...)
chebappxg.real <- function(val,grid=NULL,mapdim=NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  if(is.null(grid)) return(chebappx.real(val))
  if(is.null(dim(val))) dim(val) <- length(val)
  if(!is.list(grid) && length(grid) == length(val)) grid <- list(grid)
  if(prod(sapply(grid,length)) != length(val)) stop('grid size must match data length')
  dim(val) <- sapply(grid,length)

  # ok, grid is something like list(c(...),c(...),c(...))
  # create piecewise linear functions which maps grid-points to chebyshev grids


  intervals <- lapply(grid,range)
  # create monotone splines with splinefun, method hyman
  gridmaps <- mapply(stats::splinefun, grid, chebknots(dim(val)), MoreArgs=list(method='hyman'))
#  gridmaps <- mapply(polyh, chebknots(dim(val)), grid, MoreArgs=list(k=1))
  gridmap <- function(x) {
    if(!is.matrix(x)) return(mapply(function(gm,x) gm(x),gridmaps,x))
    apply(x,2,function(x) mapply(function(gm,x) gm(x),gridmaps,x))
  }
  ch <- chebappx.real(val)
  vectorfun(function(x,threads=getOption('chebpol.threads')) ch(gridmap(x),threads), 
            arity=length(grid), 
            domain=intervals)
}

#' @rdname chebappxg
#' @param fun The function to be approximated.
#' @param mapdim Deprecated.
#' @export
#' @keywords internal
chebappxgf <- function(...) deprecated('chebappxgf',...)
chebappxgf.real <- function(fun, grid, ..., mapdim=NULL) {

  if(!is.list(grid)) grid <- list(grid)
  chebappxg.real(evalongrid(fun, ..., grid=grid),grid,mapdim)
}
