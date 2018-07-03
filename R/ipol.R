ipol <- function(val,dims=NULL,intervals=NULL,grid=NULL,knots=NULL,k=NULL,
                 method=c('chebyshev','multilinear','fh','uniform','general','polyharmonic','rbf'),
                 ...) {
  method <- match.arg(method)
  switch(method,
         chebyshev={
           if(is.function(val) && is.null(dims)) stop('dims must be specified')
           if(is.function(val)) return(chebappxf(val,dims,intervals,...))
           return(chebappx(val,intervals))
         },
         multilinear={
           if(is.null(grid)) stop('grid must be specified for multi linear interpolation')
           return(mlappx(val,grid,...))
         },
         fh={
           if(is.null(grid)) stop('grid must be specified for Floater-Hormann interpolation')
           if(is.null(k)) k <- min(4,sapply(grid,length))
           return(fhappx(val,grid,d=k))
         },
         uniform={
           if(is.function(val) && is.null(dims)) stop('Must specify dims for uniform intervals')
           if(is.function(val)) return(ucappxf(val,dims,intervals,...))
           return(ucappx(val,intervals))
         },
         general={
           if(is.null(grid)) stop('grid must be specified for general interpolation')
           if(is.function(val)) return(chebappxgf(val,grid,...))
           return(chebappxg(val,grid))
         },
         polyharmonic={
           if(is.null(knots)) stop('Must specify knots for polyharmonic splines.')
           if(is.null(k)) k <- 3
           return(polyh(val,knots,k,...))
         },
         rbf={
           if(is.null(knots)) stop('Must specify knots for radial basis functions.')
           if(is.null(k)) k <- c(2,5,0)
           rbase <- k[1]; layers <- k[2]; lambda <- k[3]
           return(rbf.alglib(val,knots,rbase,layers,lambda))
         },
         stop('Unknown interpolation method: ', method)
         )
}
