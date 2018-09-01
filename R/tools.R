
vectorfun <- function(fun,arity,domain=NULL) {
  force(arity)
  f <- function(x) {
    mc <- match.call(expand.dots=TRUE)
    mc[[1L]] <- quote(list)
    arglist <- eval.parent(mc)
    x <- arglist[[1]]
    if(is.matrix(x) || length(x) == arity) return(do.call(fun,arglist))
    if(arity == 1) {arglist[[1]] <- matrix(x,1); return(do.call(fun,arglist))}
    # try to coerce argument to matrix
    if(length(x) %% arity == 0) {
      numvec <- length(x) / arity
      if(numvec > 1)
        warning(sprintf('Coercing vector of length %d into %dx%d matrix',length(x),arity,numvec))
      dim(x) <- c(arity,numvec)
      arglist[[1]] <- x
      do.call(fun,arglist)
    } else
      stop(sprintf('Function should take %d arguments, you supplied a vector of length %d',arity,length(x)))
  }
  formals(f) <- formals(fun)
  structure(f,arity=arity,domain=as.data.frame(domain),
            chebpol.version=utils::packageVersion('chebpol'))
}

defaultblend <- c('linear','cubic','sigmoid','parodic','square')
blenddef <- function(fun,blend=defaultblend) {
  if(is.null(blend)) return(fun)
  blend <- match.arg(blend)
  pos <- match(blend,defaultblend)
  f <- formals(fun)
  if(!('blend' %in% names(f))) return(f)
  f[['blend']] <- as.call(c(list(as.name('c')),c(defaultblend[pos],defaultblend[-pos])))
  formals(fun) <- f
  ffun <- get('fun',environment(fun))
  ff <- formals(ffun)
  if(!('blend' %in% names(ff))) return(fun)
  ff[['blend']] <- f[['blend']]
  formals(ffun) <- ff
  assign('fun',ffun,envir=environment(fun))
  fun
}

separgs <- function(fun,...) {
  fun <- match.fun(fun)
  f <- function() {
    mc <- match.call()
    mc[[1L]] <- quote(list)
    arg <- eval.parent(mc)
    nn <- sapply(arg,length)
    if(all(nn == 1)) return(fun(as.numeric(arg)))
    return(fun(do.call(rbind,arg)))
  }
  args <- as.list(sys.call())[-(1:2)]
  formals(f) <- args
  compiler::cmpfun(f)
}

# evaluate a function on a Chebyshev grid


#' Evaluate a function on a grid
#' 
#' Evaluate a function on a Chebyshev grid, or on a user-specified grid.
#' 
#' The function \code{fun} should be a \code{function(x,...)}, where
#' \code{length(x)} equals \code{length(dims)} (or \code{length(grid)}).
#' 
#' If \code{grid} is provided, \code{fun} is evaluated on each point in the
#' Cartesian product of the vectors in \code{grid}.
#' 
#' If \code{intervals} is not provided, it is assumed that the domain of the
#' function is the hypercube [-1,1] x [-1,1] x ... x [-1,1].  Thus, the
#' function is evaluated on a standard Chebyshev grid.
#' 
#' If \code{intervals} is provided, it should be a \code{list} with elements of
#' length 2, providing minimum and maximum for each dimension.
#' 
#' The grid itself may be produced by
#' \code{expand.grid(\link{chebknots}(dims,intervals))}, or
#' \code{expand.grid(grid)}.
#' 
#' This function does the same as \code{apply(expand.grid(grid),1,fun)}, but
#' it's faster and more memory-efficient for large grids because it does not
#' actually expand the grid.
#' 
#' The function \code{evalongridV} is for vectorized functions, i.e. those that
#' can take a matrix of column vectors as argument.  It's equivalent to
#' \code{fun(t(expand.grid(grid)))}.
#' 
#' @param fun Multivariate real-valued function to be evaluated. Must be
#' defined on the hypercube described by \code{intervals}.
#' @param dims A vector of integers. The number of grid-points in each
#' dimension.
#' @param intervals A list. Each entry is a vector of length 2 with the lower
#' and upper end of the interval in each dimension.
#' @param ... Further arguments to fun.
#' @param grid Rather than specifying dims and intervals to get a Chebyshev
#' grid, you may specify your own \code{grid} as a list of vectors whose
#' Cartesian product will be the grid, as in \code{expand.grid(grid)}.
#' @return An array with the value of \code{fun} on each grid point. The
#' \code{dim} attribute has been appropriately set for the grid.  If \code{fun}
#' returns a vector, this will be the first dimension of the returned array.
#' @examples
#' 
#' f <- function(x) {a <- sum(x^2); ifelse(a == 0,0,exp(-1/a))}
#' ## Standard Chebyshev grid
#' evalongrid(f,dims=c(3,5))
#' ## Then Chebyshev on [0,1] x [2,3]
#' evalongrid(f,dims=c(3,5),intervals=list(c(0,1),c(2,3)))
#' ## And on my own grid
#' grid <- list(sort(rnorm(3)),sort(rnorm(5)))
#' evalongrid(f,grid=grid)
#' g <- ipol(f,grid=grid,method='fh')
#' evalongridV(g, grid=grid, threads=2)
#' ## vector valued function
#' f <- function(x) c(prod(x),sum(x^2))
#' evalongrid(f,grid=grid)
#' 
#' @export
evalongrid <- function(fun,dims,intervals=NULL,...,grid=NULL) {
# do a call to stuff which doesn't really expand the grid
  if(is.numeric(grid)) grid <- list(grid)
  if(is.null(grid)) grid <- chebknots(dims,intervals)
  mf <- match.fun(fun)
  .Call(C_evalongrid,function(x) mf(x,...), grid)
}

#' @rdname evalongrid
#' @export
evalongridV <- function(fun, dims, intervals=NULL, ..., grid=NULL) {
  if(is.numeric(grid)) grid <- list(grid)
  if(is.null(grid)) grid <- chebknots(dims,intervals)
  fun <- match.fun(fun)
  structure(fun(t(expand.grid(grid)), ...), dim=if(length(grid) > 1) sapply(grid,length) else NULL)
}



#' Check whether chebpol uses FFTW
#' 
#' It is possible to compile chebpol without FFTW.  If this is done, it will
#' not be feasible with high-degree Chebyshev polynomials.  I.e. a 100x100x100
#' approximation will be possible, but not a one-dimensional 1000000.  This
#' function checks whether chebpol uses FFTW.
#' 
#' 
#' @return Returns TRUE if chebpol uses FFTW. Otherwise FALSE.
#'
#' @export havefftw
havefftw <- function() .Call(C_havefftw)

