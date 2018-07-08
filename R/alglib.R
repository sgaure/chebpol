
rbf.alglib <- function(val, knots, rbase=2,  layers=5, lambda=0, ...) {
  x <- threads <- NULL; rm(x,threads)
  if(is.null(dim(knots))) dim(knots) <- c(1,length(knots))
  if(is.function(val)) val <- apply(knots,2,val,...)
  model <- .Call(C_makerbf, rbind(knots,val), layers, rbase, lambda)
  local(vectorfun(.Call(C_evalrbf, model, x, threads), args=alist(x=,threads=1L), arity=nrow(knots),
                  domain=data.frame(apply(knots,1,range))),
        list(model=model))
}
havealglib <- function() .Call(C_havealglib)
