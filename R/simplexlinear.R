#' Simplex Linear approximation on scattered data
#'
#' A call \code{fun <- slappx(val, knots)} creates a piecewise multilinear interpolation
#' on the Delaunay triangulation of the \code{knots}.
#'
#' @param val Array or function. Function values in the knots, or the function
#' itself.
#' @param knots matrix. Each column is a point in M-dimensional space.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x)} interpolating the values.
#' @note
#' By default, the interpolant will yield NaN for points outside the convex hull of the knots,
#' but some extrapolation will be returned instead if the interpolant is called with the
#' argument \code{epol=TRUE}.
#' @examples
#' \dontrun{
#' knots <- matrix(runif(3*1000), 3)
#' f <- function(x) exp(-sum(x^2))
#' g <- slappx(f, knots)
#' a <- matrix(runif(3*6), 3)
#' rbind(true=apply(a,2,f), sl=g(a))
#'}
#' 
#' @export
#' @keywords internal
slappx <- function(...) deprecated('slappx',...)
slappx.real <- function(val, knots, ...) {
  if(is.function(val)) val <- apply(knots,2,val)
  dtri <- t(geometry::delaunayn(t(knots),options="Qt Pp"))
  # For fast evaluation: For each simplex, we precompute the LU-factorization
  # needed for transforming to barycentric coordinates. These are used for finding the right simplex.
  adata <- .Call(C_analyzesimplex, dtri, knots, getOption('chebpol.threads'))
  vectorfun(function(x,threads=getOption('chebpol.threads'), epol=FALSE,
                     blend=c('cubic','linear','sigmoid','parodic','square','mean')) {
    blend <- switch(match.arg(blend),linear=0L,sigmoid=1L,parodic=2L,cubic=3L,square=4L,mean=5L)
    .Call(C_evalsl, x, knots, dtri, adata, val, epol, threads, as.integer(blend))},
    arity=nrow(knots))
}

