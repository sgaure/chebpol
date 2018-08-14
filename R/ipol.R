#' Create interpolating function.
#' 
#' Create an interpolating function from given values. Several interpolation methods are
#' supported.
#' 
#' \code{ipol} is a wrapper around various interpolation methods in package \pkg{chebpol}.
#' Which arguments to specify depends on the method. The interpolation methods are described
#' in \code{vignette("chebpol",package="chebpol")}.
#' 
#' The method \code{"chebyshev"} needs only the number of Chebyshev knots in
#' each dimension. This is either the \code{"dim"} attribute of the array
#' \code{val}, or the \code{dims} argument if \code{val} is a function.  Also
#' the intervals can be specified if different from [-1, 1]. 
#' 
#' The method \code{"uniform"} is similar to the \code{"chebyshev"}, but
#' uniformly spaced knots are created.  The argument
#' \code{intervals} generally goes with \code{dims} when something else than
#' standard intervals \code{[-1, 1]} are used. 
#' 
#' The methods \code{"multilinear"}, \code{"fh"} (Floater-Hormann), \code{"stalker"}, and
#' \code{"general"} needs the argument \code{grid}.  These are the methods
#' which can use arbitrary Cartesian grids.
#' The stalker spline
#' is described in \code{vignette("stalker",package="chebpol")}. The Floater-Hormann
#' method (\code{"fh"}) also needs the \code{k} argument, the degree of the blending
#' polynomials. It defaults to 4.
#' 
#' The method \code{"polyharmonic"} needs the arguments \code{knots} and
#' \code{k}. In addition it can take the logical argument \code{normalize}
#' for normalizing the knots to the unit hypercube. The default is \code{NA}, which
#' uses normalization if any of the knots are outside the unit hypercube. Also, a logical
#' \code{nowarn} is accepted, it is used to suppress a warning in case the system can't
#' be solved exactly and a least squares fallback method is used.
#'
#' The method \code{"simplexlinear"} needs the argument \code{knots}. It creates a
#' Delaunay triangulation from the knots, and does linear interpolation in each simplex
#' by weighting the vertex values with the barycentric coordinates.
#'
#' If knots are required, but the grid argument is given, knots are constructed as
#' \code{t(expand.grid(grid))}
#' 
#' The \code{"crbf"} is the multilayer compact radial basis function
#' interpolation in ALGLIB (\url{http://www.alglib.net/interpolation/fastrbf.php}).
#' It is only available if ALGLIB was available at
#' compile time. It takes the extra arguments \code{"rbase"}, \code{"layers"}, and
#' \code{"lambda"}. These are discussed in the ALGLIB documentation.
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
#' @param ... Further arguments to the function, if \code{is.function(val)}. And some
#' extra arguments for interpolant creation described in section Details.
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
#' if(havealglib()) crb <- ipol(f, knots=knots[,1:100], method='crbf',
#'   rbase=2, layers=5, lambda=0) 
#' # make 7 points in R3 to test them on
#' m <- matrix(runif(3*7),3)
#' rbind(true=apply(m,2,f), ml=ml2(m), fh=fh2(m), cheb=ch2(m), poly=ph2(m), sl=sl2(m),
#' crbf=if(havealglib()) crb(m) else NULL )
#' 
#' @export
ipol <- function(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
                 method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic',
                          'simplexlinear', 'stalker', 'crbf'),
                 ...) {
  method <- match.arg(method)
  args <- list(...)
  switch(method,
         chebyshev={
           if(is.function(val) && is.null(dims)) stop('dims must be specified')
           if(is.function(val)) return(chebappxf.real(val,dims,intervals,...))
           if(is.null(dim(val)) && !is.null(dims)) dim(val) <- dims
           return(chebappx.real(val,intervals))
         },
         multilinear={
           if(is.null(grid)) stop('grid must be specified for multi linear interpolation')
           if(!is.list(grid)) grid <- list(grid)
           grid <- lapply(grid,as.numeric)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           return(mlappx.real(val,grid,...))
         },
         simplexlinear={
           if(is.null(knots)) {
             if(!is.null(grid)) 
               knots <- t(expand.grid(grid))
             else
               stop('knots must be specified for simplex linear interpolation')
           }
           return(slappx.real(val,knots,...))
         },
         fh={
           if(is.null(grid)) stop('grid must be specified for Floater-Hormann interpolation')
           if(!is.list(grid)) grid <- list(grid)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           grid <- lapply(grid,as.numeric)
           if(is.null(k)) k <- pmin(4,sapply(grid,length)-1)
           return(fhappx.real(val,grid,d=k))
         },
         uniform={
           if(is.function(val) && is.null(dims)) stop('Must specify dims for uniform intervals')
           if(is.null(intervals) && !is.null(grid)) {
             intervals <- lapply(grid,range)
             warning("intervals constructed from ranges of grid-argument")
           }
           if(is.function(val)) return(ucappxf.real(val,dims,intervals,...))
           if(is.null(dim(val)) && !is.null(dims)) dim(val) <- dims
           return(ucappx.real(val,intervals))
         },
         general={
           if(is.null(grid)) stop('grid must be specified for general interpolation')
           if(!is.list(grid)) grid <- list(grid)
           if(unsortedgrid(grid)) stop("grid must be distinct ordered values")
           grid <- lapply(grid,as.numeric)
           if(is.function(val)) return(chebappxgf.real(val,grid,...))
           return(chebappxg.real(val,grid))
         },
         polyharmonic={
           if(is.null(knots)) {
             if(!is.null(grid)) 
               knots <- t(expand.grid(grid))
             else
               stop('Must specify knots for polyharmonic splines.')
           }
           if(is.null(k)) k <- 3
           # Now, ... may contain normalize and nowarn in addition to
           # arguments to val if it is a function
           # peel off normalize and nowarn
           normalize <- args[['normalize']]
           if(is.null(normalize)) normalize <- NA
           nowarn <- args[['nowarn']]
           if(is.null(nowarn)) nowarn <- FALSE
           newarg <- args[-match(c('normalize','nowarn'),names(args), nomatch=0)]
           return(do.call(polyh.real,c(list(val,knots,k,normalize,nowarn),newarg)))
         },
         stalker={
           if(is.null(grid)) stop('grid must be specified for stalker interpolation')
           if(!is.list(grid)) grid <- list(grid)
           if(unsortedgrid(grid)) stop('grid must be distinct ordered values')
           grid <- lapply(grid,as.numeric)
           if(is.null(k)) k = 2
           hyman <- args[['hyman']]
           if(is.null(hyman)) hyman <- FALSE
           newarg <- args[-match('hyman',names(args),nomatch=0)]
           return(do.call(stalkerappx,c(list(val,grid,r=k,hyman=hyman),newarg)))
#           return(stalkerappx(val,grid,r=k,hyman=hyman))
         },
         crbf={
           if(is.null(knots)) stop('Must specify knots for radial basis functions.')
           if(length(args) == 0) {
               # old interface
               if(is.null(k)) {
                   rbase <- 2; layers <- 5; lambda = 0
               } else {
                   rbase <- k[1]; layers <- k[2]; lambda <- k[3]
               }
               return(rbf.alglib(val,knots,rbase,layers,lambda))
           }
           rbase <- args[['rbase']]; if(is.null(rbase)) rbase <- 2
           layers <- args[['layers']]; if(is.null(layers)) layers <- 5L
           lambda <- args[['lambda']]; if(is.null(lambda)) lambda <- 0
           newarg <- args[-match(c('rbase','layers','lambda'),names(args), nomatch=0)]
           return(do.call(rbf.alglib,c(list(val,knots,rbase,layers,lambda),newarg)))
         },
         stop('Unknown interpolation method: ', method)
         )
}

unsortedgrid <- function(g) {
  any(sapply(g,function(s) is.unsorted(s,strictly=TRUE) && is.unsorted(-s,strictly=TRUE)))
}
