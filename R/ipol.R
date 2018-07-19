#' Interpolation wrapper
#' 
#' A wrapper for all the functions in the package
#' 
#' \code{ipol} is just a wrapper around the other functions in \pkg{chebpol}.
#' Which arguments to specify depends on the method.
#' 
#' The method \code{"chebyshev"} needs only the number of Chebyshev knots to in
#' each dimension. This is either the \code{"dim"} attribute of the array
#' \code{val}, or the \code{dims} argument if \code{val} is a function.  Also
#' the intervals can be specified if different from [-1, 1]. See
#' \code{\link{chebappx}}.
#' 
#' The method \code{"uniform"} is similar to the \code{"chebyshev"}, but
#' uniformly spaced knots are created. See \code{\link{ucappx}}.  The argument
#' \code{intervals} generally goes with \code{dims} when something else than
#' standard intervals \code{[-1, 1]} are used.
#' 
#' The method \code{"multilinear"}, \code{"fh"} (Floater-Hormann), and
#' \code{"general"} needs the argument \code{grid}.  These are the methods
#' which can use arbitrary Cartesian grids.  See \code{\link{mlappx}},
#' \code{\link{fhappx}}, and \code{\link{chebappxg}}.  The Floater-Hormann
#' method (\code{"fh"}) also needs the \code{k} argument, which is passed to
#' the \code{d} argument of \code{\link{fhappx}}, the degree of the blending
#' polynomials. It defaults to 4.
#' 
#' The method \code{"polyharmonic"} needs the arguments \code{knots} and
#' \code{k}, see \code{\link{polyh}}.
#'
#' The method \code{"simplexlinear"} needs the arguments \code{knots}. It creates a
#' Delaunay triangulation from the knots, and does linear interpolation in each simplex
#' by weighting the vertex values with the barycentric coordinates, see also \code{\link{slappx}}.
#' 
#' The \code{"crbf"} is the multilayer compact radial basis function
#' interpolation in ALGLIB (\url{http://www.alglib.net/interpolation/fastrbf.php}).
#' It is only available if ALGLIB was available at
#' compile time. In this case \code{k} must be a vector of length 3, the rbase,
#' the number of layers and the non-linear smoothing parameter lambda.
#' 
#' There are also some usage examples and more in \code{vignette("chebpol")} and \code{vignette('chebusage')}.
#' 
#' @param val array or function. Function values on a grid, or the function
#' itself. If it is the values, the \code{"dim"}-attribute must be
#' appropriately set. If it is a function, it will be evaluated in the grid points.
#' @param dims integer vector. The number of grid points in each dimension. Not
#' needed if \code{val} is an array or \code{grid} is used.
#' @param intervals list of length 2 numeric vectors. The lower and upper bound
#' in each dimension. Not used if \code{grid} is specified.
#' @param grid list. Each element is a vector of ordered grid-points for a
#' dimension.  These need not be Chebyshev-knots, nor evenly spaced.
#' @param knots matrix. Each column is a point in an M-dimensional space.
#' @param k numeric. Additional value, used with some methods.
#' @param method character. The interpolation method to use.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x, threads=getOption('chebpol.threads'))} defined on a hypercube, an \link{interpolant}
#' for the given function. The argument \code{x} can be a matrix of column
#' vectors which are evaluated in parallel in a number of threads.  The
#' function yields values for arguments outside the hypercube as well, though
#' it will typically be a poor approximation.  \code{threads} is an integer
#' specifying the number of parallel threads which should be used when
#' evaluating a matrix of column vectors.
#' @examples
#' 
#' ## evenly spaced grid-points
#' su <- seq(0,1,length.out=10)
#' ## irregularly spaced grid-points
#' s <- su^3
#' ## create approximation on the irregularly spaced grid
#' ml1 <- ipol(exp, grid=list(s), method='multilin')
#' fh1 <- ipol(exp, grid=list(s), method='fh')
#' ## test it, since exp is convex, the linear approximation lies above
#' ## the exp between the grid points
#' ml1(su) - exp(su)
#' fh1(su) - exp(su)
#' 
#' ## multi dimensional approximation
#' f <- function(x) 10/(1+25*mean(x^2))
#' # a 3-dimensional 10x10x10 grid, first and third coordinate are non-uniform
#' grid <- list(s, su, sort(1-s))
#' 
#' # make multilinear, Floater-Hormann, Chebyshev and polyharmonic spline.
#' ml2 <- ipol(f, grid=grid, method='multilin')
#' fh2 <- ipol(f, grid=grid, method='fh')
#' ch2 <- ipol(f, dims=c(10,10,10), intervals=list(0:1,0:1,0:1), method='cheb')
#' knots <- matrix(runif(3*1000),3)
#' ph2 <- ipol(f, knots=knots, k=2, method='poly')
#' sl2 <- ipol(f, knots=knots, method='simplexlinear')
#' # my alglib is a bit slow, so stick to 100 knots
#' if(havealglib()) crb <- ipol(f, knots=knots[,1:100], k=c(2,5,0), method='crbf') 
#' # make 7 points in R3 to test them on
#' m <- matrix(runif(3*7),3)
#' rbind(true=apply(m,2,f), ml=ml2(m), fh=fh2(m), cheb=ch2(m), poly=ph2(m), sl=sl2(m),
#' crbf=if(havealglib()) crb(m) else NULL )
#' 
#' @export
ipol <- function(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
                 method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic',
                          'simplexlinear', 'crbf'),
                 ...) {
  method <- match.arg(method)
  switch(method,
         chebyshev={
           if(is.function(val) && is.null(dims)) stop('dims must be specified')
           if(is.function(val)) return(chebappxf(val,dims,intervals,...))
           if(is.null(dim(val)) && !is.null(dims)) dim(val) <- dims
           return(chebappx(val,intervals))
         },
         multilinear={
           if(is.null(grid)) stop('grid must be specified for multi linear interpolation')
           if(!is.list(grid)) grid <- list(grid)
           grid <- lapply(grid,as.numeric)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           return(mlappx(val,grid,...))
         },
         simplexlinear={
           if(is.null(knots)) stop('knots must be specified for simplex linear interpolation')
           return(slappx(val,knots,...))
         },
         fh={
           if(is.null(grid)) stop('grid must be specified for Floater-Hormann interpolation')
           if(!is.list(grid)) grid <- list(grid)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           grid <- lapply(grid,as.numeric)
           if(is.null(k)) k <- pmin(4,sapply(grid,length)-1)
           return(fhappx(val,grid,d=k))
         },
         uniform={
           if(is.function(val) && is.null(dims)) stop('Must specify dims for uniform intervals')
           if(is.function(val)) return(ucappxf(val,dims,intervals,...))
           if(is.null(dim(val)) && !is.null(dims)) dim(val) <- dims
           return(ucappx(val,intervals))
         },
         general={
           if(is.null(grid)) stop('grid must be specified for general interpolation')
           if(!is.list(grid)) grid <- list(grid)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           grid <- lapply(grid,as.numeric)
           if(is.function(val)) return(chebappxgf(val,grid,...))
           return(chebappxg(val,grid))
         },
         polyharmonic={
           if(is.null(knots)) stop('Must specify knots for polyharmonic splines.')
           if(is.null(k)) k <- 3
           return(polyh(val,knots,k,...))
         },
         crbf={
           if(is.null(knots)) stop('Must specify knots for radial basis functions.')
           if(is.null(k)) k <- c(2,5,0)
           rbase <- k[1]; layers <- k[2]; lambda <- k[3]
           return(rbf.alglib(val,knots,rbase,layers,lambda))
         },
         stop('Unknown interpolation method: ', method)
         )
}

unsortedgrid <- function(g) {
  any(sapply(g,function(s) is.unsorted(s,strictly=TRUE) && is.unsorted(-s,strictly=TRUE)))
}
