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
      if(is.na(smooth)) smooth <- 2-mindeg
      .Call(C_evalstalker,x,stalker,as.numeric(mindeg),as.numeric(maxdeg),
            as.numeric(smooth),threads)
  }, 
  length(grid), 
  args=alist(x=,threads=getOption('chebpol.threads'),mindeg=1,maxdeg=2,smooth=0),
  domain=lapply(grid,range)),
  list(stalker=stalker))
  
}
