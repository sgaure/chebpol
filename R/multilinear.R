#' Multilinear interpolation on a grid
#' 
#' Multilinear interpolation on an arbitrary Cartesian product.
#' 
#' A call \code{fun <- mlappx(val,grid)} creates a multilinear interpolation on
#' the grid.  The value on the grid points will be exact, the value between the
#' grid points is a convex combination of the values in the corners of the
#' hypercube surrounding it.
#' 
#' If \code{val} is a function it will be evaluated on the grid.
#' 
#' @aliases mlappx mlappxf
#' @param val Array or function. Function values on a grid, or the function
#' itself. If it is the values, the \code{dim}-attribute must be appropriately
#' set.
#' @param grid A list.  Each element is a vector of ordered grid-points for a
#' dimension.  These need not be Chebyshev-knots, nor evenly spaced.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x)} defined on the hypercube, approximating the
#' given function.  The function yields values for arguments outside the
#' hypercube as well, as a linear extension.
#' @examples
#' 
#' ## evenly spaced grid-points
#' su <- seq(0,1,length.out=10)
#' ## irregularly spaced grid-points
#' s <- su^3
#' ## create approximation on the irregularly spaced grid
#' ml1 <- Vectorize(mlappx(exp,list(s)))
#' ## test it, since exp is convex, the linear approximation lies above
#' ## the exp between the grid points
#' ml1(su) - exp(su)
#' 
#' ## multi linear approx
#' f <- function(x) exp(sum(x^2))
#' grid <- list(s,su)
#' 
#' ml2 <- mlappx(evalongrid(f,grid=grid),grid)
#' # an equivalent would be ml2 <- mlappx(f,grid)
#' 
#' a <- runif(2); ml2(a); f(a)
#' # we also get an approximation outside of the domain, of disputable quality
#' ml2(c(1,2)); f(c(1,2))
#' 
#' @export mlappx
mlappx <- function(val, grid, ...) {
  x <- threads <- NULL; rm(x,threads) # avoid cran check warning 
  if(is.numeric(grid)) grid <- list(grid)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)
  val <- as.numeric(val)
#  if(adjust!=0) {
#    val <- val + (val - .Call(C_predmlip,grid,as.numeric(val)))*adjust
#  }
  local(vectorfun(.Call(C_evalmlip,grid,val,x,threads), length(grid), 
                  args=alist(x=,threads=getOption('chebpol.threads')),
                  domain=lapply(grid,range)),
        list(grid=grid,val=val))
}
