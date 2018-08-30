#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/Visibility.h>
#include "config.h"

SEXP R_makerbf(SEXP, SEXP, SEXP, SEXP);
SEXP R_evalrbf(SEXP, SEXP, SEXP);
SEXP R_evalstalker(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_computehyp(SEXP, SEXP);
SEXP R_evalhyp(SEXP, SEXP, SEXP, SEXP);
SEXP R_havealglib();
SEXP havegsl();
