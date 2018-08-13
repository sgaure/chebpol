# Precompute some values used in the stalker evaluation
stalkercontext <- function(dmin,dplus,r) {
  rank <- seq_along(dmin)
  powmin <- lapply(rank, function(i) dmin[[i]]^r[[i]])
  powplus <- lapply(rank, function(i) dplus[[i]]^r[[i]])
  det <- lapply(rank, function(i) {
    1/(dmin[[i]]*powplus[[i]] + dplus[[i]]*powmin[[i]])
  })
  list(det=det,pmin=powmin,pplus=powplus)
}
stalkerappx <- function(val, grid, r=2, ...) {
  if(is.numeric(grid)) grid <- list(grid)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)

  r <- rep(as.numeric(r),length.out=length(grid))
  # are any of the dimensions non-uniform?
  dmin <- lapply(grid, function(gr) c(NaN,diff(gr)))
  dplus <- lapply(dmin, function(dd) c(dd[-1],NaN))
  val <- as.numeric(val)
  stalker <- c(list(val=val, grid=grid), stalkercontext(dmin,dplus,r))

  vectorfun(function(x,threads=getOption('chebpol.threads'),degree=r) {
    degree <- rep(as.numeric(degree),length.out=length(grid))
    if(!identical(degree,r))
        stalker <- c(list(val=val, grid=grid), stalkercontext(dmin,dplus,degree))
    .Call(C_evalstalker,x,stalker,degree, as.integer(threads))
  }, 
  arity=length(grid), 
  domain=lapply(grid,range))
  
}
