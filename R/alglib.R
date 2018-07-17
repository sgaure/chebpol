
rbf.alglib <- function(val, knots, rbase=2,  layers=5, lambda=0, ...) {
  x <- threads <- NULL; rm(x,threads)
  if(is.null(dim(knots))) dim(knots) <- c(1,length(knots))
  if(is.function(val)) val <- apply(knots,2,val,...)
  model <- .Call(C_makerbf, rbind(knots,val), layers, rbase, lambda)
  local(vectorfun(.Call(C_evalrbf, model, x, threads), args=alist(x=,threads=1L), arity=nrow(knots),
                  domain=data.frame(apply(knots,1,range))),
        list(model=model))
}
#' Check whether chebpol has the ALGLIB library
#'
#' If ALGLIB was available at compile time, it can be used for compact
#' support radial basis function interpolations on scattered data.
#' This function checks whether ALGLIB is available.
#'
#' @return Returns TRUE if ALGLIB is available. Otherwise FALSE.
#' @export
havealglib <- function() .Call(C_havealglib)
