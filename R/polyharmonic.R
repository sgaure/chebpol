#' Polyharmonic splines on scattered data
#' 
#' Polyharmonic splines on scattered data.
#' 
#' \code{polyh} fits a polyharmonic spline with radial basis function
#' \code{x^k} for odd \code{k}, and \code{x^k log(x)} for even \code{k}. If
#' \code{k < 0}, the basis \code{exp(k x^2)} is used. There are more details in
#' a vignette.
#' 
#' If \code{val} is a function it will be evaluated on the knots.
#' 
#' \code{normalize} can be used to change the scaling on the space. Set
#' \code{normalize=TRUE} to do an affine transformation on the knots into the
#' unit hybercube. The default is to transform if any of the knot-coordinates
#' are outside the interval \eqn{[0,1]}. You may also specify \code{normalize}
#' as an \eqn{M x 2} matrix, where the columns \eqn{a} and \eqn{b} are used for
#' the normalization: \eqn{x -> a*(x-b)}.  \code{normalize} can be set to a
#' vector of length 2 \code{c(a, b)} if the same normalization should apply in
#' each dimension.  Set \code{normalize=FALSE} if you do not want any scaling.
#' 
#' @param val array or function. Function values on scattered data, or the
#' function itself.
#' @param knots matrix. Each column is a point in an M-dimensional space.
#' @param k positive integer or negative numeric. The degree of the
#' polyharmonic spline.
#' @param normalize logical, vector or matrix. Should coordinates be
#' normalized?
#' @param nowarn logical. Avoid warning about fallback to least squares fit.
#' @param ... Further arguments to the function, if \code{is.function(val)}.
#' @return A \code{function(x)} defined on the multidimensional space,
#' approximating the given function.
#' @examples
#' \dontrun{
#' # a function on a 20-dimensional space
#' r <- runif(20)
#' r <- r/sum(r)
#' f <- function(x) 1/mean(log1p(r*x))
#' # 1000 random knots 
#' knots <- matrix(runif(20000), 20)
#' phs <- polyh(f, knots, 3)
#' # test it in a random point
#' s <- runif(20)
#' c(true=f(s), phs(s))
#' }
#' @export
#' @keywords internal
polyh <- function(...) deprecated('polyh',...)
polyh.real <- function(val, knots, k=2, normalize=NA, nowarn=FALSE, ...) {
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

  # trickery to get it in place
  phi <- function(x) {
    eval.parent(as.call(list(quote(.Call), C_phifunc, substitute(x), k, getOption('chebpol.threads'))))
  }

  #A <- phi(apply(knots,2, function(ck) colSums((ck-knots)^2)))
  # one day I will look into a faster solver, 
  A <- phi(.Call(C_sqdiffs,knots,knots,getOption('chebpol.threads')))  
  B <- rbind(1,knots)
  if(FALSE) {
    # Outta https://mathematica.stackexchange.com/questions/65763/understanding-polyharmonic-splines
    B <- t(B)
    Ai <- solve(A)
    v <- solve(crossprod(B,Ai %*% B), crossprod(B,Ai %*% val))
    w <- Ai %*% (val - B %*% v)
#
#    v <- solve(crossprod(B,solve(A,B)), crossprod(B,solve(A, val)))
#    w <- solve(A,val - B %*% v)
  } else {
    mat <- matrix(0,nrow(A)+nrow(B),nrow(A)+nrow(B))
    mat[1:nrow(A),1:ncol(A)] <- A
    mat[(nrow(A)+1):nrow(mat), 1:ncol(B)] <- B
    mat[1:ncol(B), (ncol(A)+1):ncol(mat)] <- t(B)
#    mat <- cbind(rbind(A,B),rbind(t(B),matrix(0,M+1,M+1)))
    rhs <- c(val,rep(0,M+1))
    wv <- try(solve(mat, rhs), silent=TRUE)
    if(inherits(wv,'try-error')) {
      if(!nowarn)
        warning('Failed to fit exactly, fallback to least squares fit.',
                if(!normalize) ' You could try normalize=TRUE.' else '')
      wv <- stats::lm.fit(mat,rhs)$coefficients
      wv[is.na(wv)] <- 0
    }
    w <- wv[1:N]
    v <- wv[(N+1):length(wv)]
    rm(mat,rhs,wv)
  }
  rm(A,B)
  vectorfun(function(x,threads=getOption('chebpol.threads')) .Call(C_evalpolyh, normfun(x), knots, w, v, k, threads, NULL),
            arity=M,
            domain=data.frame(apply(knots,1,range)))
}
