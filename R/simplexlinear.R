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
#' @note The simplex linear interpolation uses package \pkg{geometry} heavily, it is not 
#' optimized for dimensions larger than 2, so it is quite slow, in particular evaluation of the
#' interpolant. It does not yield values for points outside the convex hull of the knots.
#' @examples
#' knots <- matrix(runif(3*100), 3)
#' f <- function(x) exp(-sum(x^2))
#' g <- slappx(f, knots)
#' a <- matrix(runif(3*5), 3)
#' rbind(true=apply(a,2,f), sl=g(a))
#'
#' @export
slappx <- function(val, knots, ...) {
  x <- threads <- NULL; rm(x,threads)
  if(is.function(val)) val <- apply(knots,2,val)
  tknots <- t(knots)
  tess <- geometry::delaunayn(tknots)
  local(vectorfun(sleval(x, val, tess, tknots, threads, NULL),
                  arity=nrow(knots), args=alist(x=, threads=getOption('chebpol.threads'))),
        list(val=val, tknots=tknots, tess=tess))
}

sleval <- function(x, val, tess, tknots, threads, spare) {
  idx <- geometry::tsearchn(tknots, tess, t(x))
  sapply(seq_along(idx$idx), function(i) sum(val[tess[idx$idx[i],]]*idx$p[i,]))
}
