stalkerappx <- function(val, grid, ...) {
  if(is.numeric(grid)) grid <- list(grid)
  if(any(sapply(grid,is.unsorted))) {
    if(!is.function(val)) stop('Grid points must be ordered in increasing order')
    grid <- lapply(grid,sort)
  }
  if(is.function(val)) val <- evalongrid(val,grid=grid,...)
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)

  # are any of the dimensions non-uniform?


  ## unimap <- identity
  ## if(any(abs(unlist(lapply(grid,diff,differences=2L))) > 1e-12)) {
  ##   # create monotonic maps into a uniform grid
  ##   unimaps <- lapply(grid, function(gr) {
  ##     if(all(abs(diff(gr,differences=2L) < 1e-12))) return(identity)
  ##     return(splinefun(gr,seq_along(gr), method='hyman'))
  ##   })
  ##   grid <- lapply(seq_along(grid), function(i) unimaps[[i]](grid[[i]]));
  ##   unimap <- function(x) {
  ##     if(!is.matrix(x)) x <- matrix(x,length(grid))
  ##     t(sapply(seq_along(grid), function(i) sapply(x[i,],unimaps[[i]])))
  ##   }
  ## } 

  stalker <- list(val=as.numeric(val), grid=grid)
  vectorfun(function(x,threads=getOption('chebpol.threads'),mindeg=1,maxdeg=2) {
    mindeg <- rep(mindeg,length.out=length(grid))
    maxdeg <- rep(maxdeg,length.out=length(grid))
    .Call(C_evalstalker,x,stalker,as.numeric(mindeg),as.numeric(maxdeg),
          as.integer(threads))
  }, 
  length(grid), 
  domain=lapply(grid,range))
  
}
