\name{chebcoef}
\alias{chebcoef}
\title{Compute Chebyshev-coefficients given values on a Chebyshev grid}
\concept{Chebyshev polynomial}
\description{
  Compute the multivariate Chebyshev-coefficients, given values on a Chebyshev grid.
}

\usage{
chebcoef(val, dct=FALSE)
}

\arguments{
  \item{val}{An \code{array} of function values on a Chebyshev grid.
    The \code{dim}-attribute must be appropriately set.
    If not set, it is assumed to be one-dimensional.
  }
  \item{dct}{Logical. Since the Chebyshev coefficients are closely related
    to the DCT-II transform of \code{val}, the non-normalized real-even
    DCT-II coefficients may be retrieved instead. I.e. those from
    FFTW_REDFT10 in each dimension.
    This is not used anywhere in the package, it is merely provided as a convenience
    for those who might need it.
  }
}

\value{
  An array of Chebyshev-coefficients for an interpolating Chebyshev-polynomial.
}

\details{
  If \code{val} has no \code{dim}-attribute, it is assumed to be
  one-dimensional of length the length of \code{val}.

  If \pkg{chebpol} was compiled without \acronym{FFTW}, running
  \code{chebcoef} on large grids may be slow and memory-demanding.
}
\seealso{\code{\link{havefftw}}}
\examples{
## Coefficients for a 2x3x4 grid
a <- array(rnorm(24),dim=c(2,3,4))
chebcoef(a)
}