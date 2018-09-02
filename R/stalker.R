# Precompute some values used in the stalker evaluation
stalkercontext <- function(dmin,dplus,r,grid) {
  rank <- seq_along(dmin)
  powmin <- lapply(rank, function(i) dmin[[i]]^r[[i]])
  powplus <- lapply(rank, function(i) dplus[[i]]^r[[i]])
  det <- lapply(rank, function(i) {
    1/(dmin[[i]]*powplus[[i]] + dplus[[i]]*powmin[[i]])
  })
  list(det=det,pmin=powmin,pplus=powplus)
}
stalkerappx <- function(val, grid, r=2, ...) {
  if(!is.list(grid)) grid <- list(grid)
  grid <- lapply(grid,as.numeric)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)

  r <- rep(as.numeric(r),length.out=length(grid))
  uniform <- sapply(grid,function(g) all(abs(diff(g,differences=2L)) < 1e-8))
  # precompute distance between grid points
  dmin <- lapply(grid, function(gr) c(NaN,diff(gr)))
  dplus <- lapply(dmin, function(dd) c(dd[-1],NaN))
  val <- as.numeric(val)
  stalker <- c(list(val=val, grid=grid), stalkercontext(dmin,dplus,r,grid),uniform=list(uniform))

  if(!havegsl() && any(is.na(r) & !uniform)) {
    stop('Non-uniform grid and varying degree needs GSL. Recompile with GSL.')
  }

  vectorfun(function(x,threads=getOption('chebpol.threads'),degree=r,
                     blend=c('linear','cubic','sigmoid','parodic','square')) {
    blend <- switch(match.arg(blend),linear=0L,sigmoid=1L,parodic=2L,cubic=3L,square=4L)
    if(!identical(degree,r)) {
      degree <- rep(as.numeric(degree),length.out=length(grid))
      if(!havegsl() && any(is.na(degree) & !uniform)) {
        stop('Non-uniform grid and varying degree needs GSL. Recompile with GSL.')
      }
      stalker <- c(list(val=val, grid=grid), stalkercontext(dmin,dplus,degree,grid),uniform=list(uniform))
    }
    .Call(C_evalstalker,x,stalker,degree, as.integer(blend), as.integer(threads))
  }, 
  arity=length(grid), 
  domain=lapply(grid,range))
}

hstalkerappx <- function(val, grid, ...) {
  if(!is.list(grid)) grid <- list(grid)
  grid <- lapply(grid,as.numeric)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  stalker <- .Call(C_makehyp,val,grid)
  vectorfun(function(x,threads=getOption('chebpol.threads'),
                     blend=c('cubic','linear','sigmoid','parodic','square')) {
    blend <- switch(match.arg(blend),linear=0L,sigmoid=1L,parodic=2L,cubic=3L,square=4L)
    .Call(C_evalhyp,x,stalker,as.integer(blend),as.integer(threads))
  },
  arity=length(grid),
  domain=lapply(grid,range))
}

newstalkerappx <- function(val, grid, ...) {
  if(!is.list(grid)) grid <- list(grid)
  grid <- lapply(grid,as.numeric)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  stalker <- .Call(C_makestalk,val,grid)
  vectorfun(function(x,threads=getOption('chebpol.threads'),degree=r,
                     blend=c('cubic','linear','sigmoid','parodic','square')) {
    blend <- switch(match.arg(blend),linear=0L,sigmoid=1L,parodic=2L,cubic=3L,square=4L)
    .Call(C_evalstalk,x,stalker,as.integer(blend),as.integer(threads))
  },
  arity=length(grid),
  domain=lapply(grid,range))
}

havegsl <- function() .Call(C_havegsl)
