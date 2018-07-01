ipol <- function(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
                 method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic'),
                 ...) {
  method <- match.arg(method)
  switch(method,
         chebyshev={
           if(is.function(val)) return(chebappxf(val,dims,intervals,...))
           return(chebappx(val,intervals))
         },
         multilinear={
           return(mlappx(val,grid,...))
         },
         fh={
           if(is.null(k)) k <- min(4,sapply(grid,length))
           return(fhappx(val,grid,d=k))
         },
         uniform={
           if(is.function(val) && is.null(dims)) stop('Must specify dims for uniform intervals')
           if(is.function(val)) return(ucappxf(val,dims,intervals,...))
           return(ucappx(val,intervals))
         },
         general={
           if(is.function(val)) return(chebappxgf(val,grid,...))
           return(chebappxg(val,grid))
         },
         polyharmonic={
           if(is.null(knots)) stop('Must specify knots for polyharmonic splines.')
           return(polyh(val,knots,k,...))
         },
         stop('Unknown interpolation method: ', method)
         )
}
