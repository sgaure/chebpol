.onAttach <- function(libname,pkgname) {
  if(!havefftw()) {
    packageStartupMessage("*** ",pkgname,": FFTW not used.\n*** You should install it from http://fftw.org\n*** or check if your OS-distribution provides it, and recompile.")
  }
}

.onLoad <- function(libname,pkgname) {
  if(is.na(thr <- as.integer(Sys.getenv('CHEBPOL_THREADS')))) thr <- 1L
  options(chebpol.threads=thr)
}

# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots1 <- function(n, interval=NULL) {
  kn <- cos(pi*((1:n)-0.5)/n)
  if((n %% 2) == 1) kn[[(n+1)/2]] = 0
  if(is.null(interval)) return(kn)
  kn*diff(interval)/2 + mean(interval)
}

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

# evaluate a function on a Chebyshev grid
evalongrid <- function(fun,dims,intervals=NULL,...,grid=NULL) {
# do a call to stuff which doesn't really expand the grid
  if(is.numeric(grid)) grid <- list(grid)
  if(is.null(grid)) grid <- chebknots(dims,intervals)
  mf <- match.fun(fun)
  .Call(C_evalongrid,function(x) mf(x,...), grid)
#  structure(apply(expand.grid(chebknots(dims,intervals)),1,fun,...), dim=dims)
}


# Chebyshev coefficients for x, which may be an array
chebcoef <- function(val, dct=FALSE) {
  structure(.Call(C_chebcoef,as.array(val),dct),dimnames=dimnames(val))
}

chebeval <- function(x,coef,intervals=NULL,threads=getOption('chebpol.threads')) {
  if(is.null(intervals)) return(.Call(C_evalcheb,coef,x,threads))
  # map into intervals
  .Call(C_evalcheb,coef,mapply(function(x,i) 2*(x[[1]]-mean(i))/diff(i),x,intervals),threads)
}

# return a function which is a Chebyshev interpolation
chebappx <- function(val,intervals=NULL) {
  if(is.null(dim(val))) {
    # allow for one-dimensional
    dim(val) <- length(val)
  }
  K <- length(dim(val))
   # allow for vector, e.g. intervals=c(0,1), put it inside list
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  cf <- chebcoef(val)

  x <- threads <- NULL; rm(x,threads) # avoid cran check warning 
  if(is.null(intervals)) {
    # it's [-1,1] intervals, so drop transformation
    fun <- structure(vectorfun(.Call(C_evalcheb,cf,x,threads), K,
                               args=alist(x=,threads=getOption('chebpol.threads'))),
                     arity=K)
    rm(val)
  } else {
    # it's intervals, create mapping into [-1,1]
    if(!is.list(intervals)) stop("intervals should be a list")
    if(any(sapply(intervals,length) != 2)) stop("interval elements should have length 2")
    if(length(intervals) != length(dim(val))) stop("values should have the same dimension as intervals ",
               length(intervals),' ',length(dim(val)))
    ispan <- sapply(intervals,function(x) 2/diff(x))
    mid <- sapply(intervals,function(x) mean(x))
    imap <- compiler::cmpfun(function(x) (x-mid)*ispan)

    fun <- structure(vectorfun(.Call(C_evalcheb,cf,imap(x), threads), K,
                               args=alist(x=,threads=getOption('chebpol.threads'))),
                     arity=length(dim(val)),domain=intervals)
    rm(val)
  }
  fun
}

# interpolate a function
chebappxf <- function(fun,dims,intervals=NULL,...) {
  chebappx(evalongrid(fun,dims,intervals,...),intervals)
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



chebappxg <- function(val,grid=NULL,mapdim=NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  x <- threads <- NULL; rm(x,threads) # avoid cran check warning 
  if(is.null(grid)) return(chebappx(val))
  if(is.null(dim(val))) dim(val) <- length(val)
  if(!is.list(grid) && length(grid) == length(val)) grid <- list(grid)
  if(prod(sapply(grid,length)) != length(val)) stop('grid size must match data length')
  dim(val) <- sapply(grid,length)

  # ok, grid is something like list(c(...),c(...),c(...))
  # create piecewise linear functions which maps grid-points to chebyshev grids


  intervals <- lapply(grid,function(x) c(min(x),max(x)))
  # create monotone splines with splinefun, method monoH.FC
  gridmaps <- mapply(splinefun, grid, chebknots(dim(val)), MoreArgs=list(method='monoH.FC'))
#  gridmaps <- mapply(polyh, chebknots(dim(val)), grid, MoreArgs=list(k=1))
  gridmap <- function(x) {
    if(!is.matrix(x)) return(mapply(function(gm,x) gm(x),gridmaps,x))
    apply(x,2,function(x) mapply(function(gm,x) gm(x),gridmaps,x))
  }
  ch <- chebappx(val)
  structure(vectorfun(ch(gridmap(x),threads), length(grid), args=alist(x=,threads=getOption('chebpol.threads'))),
            arity=length(grid),domain=intervals,grid=grid)
}

chebappxgf <- function(fun, grid, ..., mapdim=NULL) {

  if(!is.list(grid)) grid <- list(grid)
  chebappxg(evalongrid(fun, ..., grid=grid),grid,mapdim)
}


# General grids, Floater-Hormann
fhappx <- function(val,grid=NULL, d=1, ...) {
  x <- threads <- NULL; rm(x,threads) # avoid warning about undefined vars
  if(is.null(grid)) 
    stop('Must specify grid')
  if(!is.list(grid)) grid <- list(grid)
  if(is.function(val)) val <- evalongrid(val, grid=grid, ...)
  dd <- as.integer(d)
  dd <- rep(dd, length(grid) %/% length(d))
  # calculate weights, formula 18 in Floater & Hormann
  weights <- lapply(seq_along(grid), function(gg) {
    g <- grid[[gg]]
    d <- dd[gg]
    n = length(g)-1
    if(d > n+1) stop(sprintf('d (%d) must be less or equal to dimension (%d)',d,n+1))
    sapply(seq_along(g)-1L, function(k) {
      gk <- g[k+1]
      sum(sapply(intersect(seq(k-d, k), 0:(n-d)), function(i) {
        sign = if(i %% 2 == 0) 1 else -1
        pidx <- setdiff(seq(i,i+d),k)
        sign/prod(gk - g[pidx+1])
      }))
    })
  })
#  weights <- lapply(grid, function(g) 1/sapply(seq_along(g), function(i) prod((g[i] - g[-i]))))
#  weights <- NULL
  vectorfun(.Call(C_FH,x,val,grid,weights,threads), 
            args=alist(x=,threads=getOption('chebpol.threads')),
            arity=length(grid))
}

# we can actually find the grid-maps for uniform grids.
# the Chebyshev knots are cos(pi*(j+0.5)/n) for j=0..n-1 These should
# map into the n grid points. These have distance 2/(n-1), starting in -1, ending in 1
# so they are -1 + 2*j/(n-1). After some manipulation, the function is:

ugm <- function(x,n) sin(0.5*pi*x*(1-n)/n)

ucappx <- function(val, intervals=NULL) {
  x <- threads <- NULL; rm(x,threads) # avoid cran check warning 
  if(is.null(dim(val))) dim(val) <- length(val)
  dims <- dim(val)
  ch <- chebappx(val)
  if(is.null(intervals)) {
    gridmap <- function(x) mapply(function(xi,d) ugm(xi,d),x,dims)
  } else {
    # precompute interval mid points and inverse lengths
    md <- lapply(intervals,mean)
    ispan <- lapply(intervals, function(i) 2/diff(i))
    gridmap <- function(x) mapply(function(xi,mid,is,d) ugm(is*(xi-mid),d),x,md,ispan,dims)
  }
  gm <- function(x) {
    if(is.matrix(x)) apply(x,2,gridmap) else gridmap(x)
  }
  return(structure(vectorfun(ch(gm(x), threads), length(dims), 
                             args=alist(x=,threads=getOption('chebpol.threads'))), 
                   arity=length(dims)))
}

ucappxf <- function(fun, dims, intervals=NULL,...) {
  if(is.null(intervals))
    return(ucappx(evalongrid(fun,...,grid=lapply(dims,function(d) seq(-1,1,length.out=d)))))
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  return(
    ucappx(evalongrid(fun,...,
                      grid=mapply(function(d,i) seq(min(i),max(i),length.out=d),
                                  dims, intervals,SIMPLIFY=FALSE)),
           intervals))
}

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
  vectorfun(.Call(C_evalmlip,grid,val,x,threads), length(grid), 
            args=alist(x=,threads=getOption('chebpol.threads')))
}

havefftw <- function() .Call(C_havefftw)

polyh <- function(val, knots, k=2, normalize=NA, nowarn=FALSE, ...) {
# Linear polyharmonic splines. Centres are columns in matrix knots. Function values in val.
# Quite slow for ncol(knots) > 3000 or so.  k=2 yields thin-plate splines.
# There exist faster evaluation methods for dimensions <= 4
# k < 0 yields Gaussian kernel splines, with sigma^2 = -1/k
# we compute r^2, so powers etc. are adjusted for that case
# Hmm, I should take a look at http://dx.doi.org/10.1016/j.jat.2012.11.008 for the unit ball
# perhaps some normalization should be added?
# I also need to look at fast summation with NFFT
# https://www-user.tu-chemnitz.de/~potts/nfft/fastsum.php
# Likewise its interpolation applications which I don't understand yet.
  if(is.null(dim(knots))) dim(knots) <- c(1,length(knots))
  if(is.function(val)) val <- apply(knots,2,val,...)
  N <- ncol(knots)
  M <- nrow(knots)
  if(k > 0) if(abs(round(k)-k) > sqrt(.Machine$double.eps))
              stop(sprintf('k is positive, but not an integer: %.17f',k))
            else 
              k <- as.integer(round(k))

  if(is.na(normalize)) normalize <- (min(knots) < 0) || (max(knots) > 1)
  else if(length(normalize) == 2L) dim(normalize) <- 1:2
  if(is.matrix(normalize)) {
    if(M %% nrow(normalize) != 0 || ncol(normalize) != 2) 
      stop('normalize must but be a n x 2 matrix with ', M, ' divisible by n(', nrow(normalize), ')')
    wa <- normalize[,1]
    wb <- normalize[,2]
    normfun <- function(x) wa*(x-wb)
    normalize <- TRUE
  } else if(normalize) {
    wa <- as.vector(1/diff(apply(knots,1,range)))
    wb <- apply(knots,1,min)
    normfun <- function(x) wa*(x-wb)
    knots <- normfun(knots)
    normalize <- TRUE
  } else {
    normfun <- identity
    normalize <- FALSE
  }

  ki <- k/2
  if(k < 0) {
    phi <- local(compiler::cmpfun(function(r2) exp(k*r2)),list(k=k))
  } else if(k %% 2L == 1L) {
    phi <- local(compiler::cmpfun(function(r2) r2^ki), list(ki=ki))
  } else {
    ki <- as.integer(ki)-1L # trick to handle r2=0. Works because 0^0 = 1 in R
    phi <- local(compiler::cmpfun(function(r2) r2^ki * log(r2^r2)),list(ki=ki))
  }

  # would it be faster to apply phi only on lower tri?  Without crossprod?
  # and either fill in the upper tri, or tailor a solver? I think
  # evaluation of phi on the matrix is fast compared to creating the matrix,
  # though I haven't measured it. Could've done it in C, probably somewhat faster.
  # use abs when we call phi. The argument can't be negative, but numerically it can
  sqnm <- colSums(knots^2)
  #A <- phi(apply(knots,2, function(ck) colSums((ck-knots)^2)))
#  A <- phi(.Call(C_sqdiffs,knots,knots))  
  # faster:
  A <- phi(abs(-2*crossprod(knots) + sqnm + rep(sqnm,each=N)))
  diag(A) <- phi(0) 

  B <- rbind(1,knots)
  mat <- cbind(rbind(A,B),rbind(t(B),matrix(0,M+1,M+1)))
  rhs <- c(val,rep(0,M+1))
  wv <- try(solve(mat, rhs, tol=10*.Machine$double.eps), silent=TRUE)
  if(inherits(wv,'try-error')) {
    if(!nowarn)
      warning('Failed to fit exactly, fallback to least squares fit.',
              if(!normalize) ' You could try normalize=TRUE.' else '')
    wv <- lm.fit(mat,rhs)$coefficients
    wv[is.na(wv)] <- 0
  }
  w <- wv[1:N]
  v <- wv[(N+1):length(wv)]

  # There's an alternative here: https://mathematica.stackexchange.com/questions/65763/understanding-polyharmonic-splines
  # with W = B'
  # solve W' A^{-1} W v = W' A^{-1} val  for v
  # compute w as A^{-1} val - A^{-1} W v
  # It's slower than the above, at least with MKL.
  # W <- t(rbind(1,knots))
  # Ai <- solve(A)
  # AiW <- crossprod(Ai,W)
  # v <- as.numeric(solve(crossprod(W,AiW),t(AiW) %*% val))
  # w <- as.numeric(Ai %*% val - AiW %*% v)

  local(function(x, threads) {
    if(!missing(threads)) warning('polyh does not support parallel threads')
    if(!is.matrix(x) && length(x) == M) {
      nx <- normfun(x)
      sum(w*phi(abs(-2*crossprod(nx,knots) + sqnm + sum(nx^2)))) + sum(v*c(1,nx))
    } else {
      if(is.null(dim(x))) dim(x) <- c(M,length(x)/M)
      if(nrow(x) != M) stop('spline was built for dimension ',M,' not ',nrow(x))
      nx <- normfun(x)
      colSums(w*phi(abs(-2*crossprod(knots,nx) + sqnm + rep(colSums(nx^2),each=N)))) + colSums(v*rbind(1,nx))
    }
  }, list(w=w,v=v,knots=knots,phi=phi,sqnm=sqnm,M=M,normfun=normfun))
}


vectorfun <- function(e,arity,args=alist(x=)) {
  fun <- function() {}
  formals(fun) <- args
  body(fun) <- substitute(e)
  environment(fun) <- parent.frame()
  force(arity)
  f <- local(function(x) {
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
  }, list(fun=compiler::cmpfun(fun)))
  formals(f) <- args
  f
}
