#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/Visibility.h>

SEXP R_makerbf(SEXP, SEXP, SEXP, SEXP);
SEXP R_evalrbf(SEXP, SEXP, SEXP);
SEXP R_evalstalker(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_makestalker(SEXP, SEXP, SEXP);
SEXP R_havealglib();
