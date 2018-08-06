stalkerappx <- function(val, grid, ...) {
  x <- threads <- maxdeg <- mindeg <- smooth <- NULL; rm(x,threads,maxdeg,mindeg,smooth) # avoid cran check warning 
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
  stalker <- .Call(C_makestalker, val, grid, getOption('chebpol.threads'))
  local(vectorfun({
    .Call(C_evalstalker,x,stalker,mindeg,maxdeg,as.logical(smooth),threads)
  }, 
  length(grid), 
  args=alist(x=,threads=getOption('chebpol.threads'),mindeg=1,maxdeg=2,smooth=FALSE),
  domain=lapply(grid,range)),
  list(stalker=stalker))
  
}
