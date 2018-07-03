#ifdef HAVE_ALGLIB
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Visibility.h>
#include "alglib/src/stdafx.h"
#include "alglib/src/interpolation.h"
#include "config.h"


/* http://www.alglib.net/translator/man/manual.cpp.html#sub_rbfcreate
   The sequence is:
   rbfcreate()
   rbfsetpoints()
   Optional: rbfsetalgohierarchical()/rbfsetalgomultilayer()
   rbfbuildmodel()

   rbfcalc()
*/

// alglib::xparams attribute_visible xdefault = {0};

SEXP R_makerbf(SEXP Sknots, SEXP Svalues) {
  if(!isMatrix(Sknots)) error("knots must be a matrix");
  if(ncols(Sknots) != length(Svalues)) error("columns of knot-matrix must match number of values");
  alglib::ae_int_t nx = nrows(Sknots);
  int numvecs = ncols(Sknots);
  alglib::ae_int_t ny = 1;  // not vector valued yet
  alglib::ae_int_t N = length(Svalues);
  alglib::ae_int_t nlayers = 10;
  double lambdaN = 0.3;
  double rbase = 2.0;
  alglib::rbfmodel s;
  alglib::rbfcreate(nx, ny, s);
  // bless it, it wants the knots and values cbind'ed, and in its private matrix format real_2d_array
  double *knotvalues = (double *) R_alloc( (nx+ny)*N, sizeof(double));
  double *knots = REAL(Sknots);
  double *values = REAL(Svalues);

  // fill in the knotvalues from knots and values

  alglib::real_2d_array dataset;
  dataset.attach_to_ptr(N, nx+ny, knotvalues);
  alglib::rbfsetpoints(s, dataset, N);

  alglib::rbfsetalgohierarchical(s, rbase, nlayers, lambdaN);

  // return some structure which agrees with R's gc
}

SEXP R_evalrbf(SEXP vecs) {
}
#endif
