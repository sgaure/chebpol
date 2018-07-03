#include "config.h"
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
#include <exception>

using namespace alglib;
using namespace std;

extern "C" {
  #include "chebpol.h"
}
/* http://www.alglib.net/translator/man/manual.cpp.html#sub_rbfcreate
   The sequence is:
   rbfcreate()
   rbfsetpoints()
   Optional: rbfsetalgohierarchical()/rbfsetalgomultilayer()
   rbfbuildmodel()

   rbfcalc()
*/

extern "C" {

  static void mfinal(SEXP model) {
    rbfmodel *s = (rbfmodel*) R_ExternalPtrAddr(model);
    delete s;
  }
  
  SEXP R_makerbf(SEXP Sknotvalues, SEXP Slay,  SEXP Srbase, SEXP lambda) {
    if(!isMatrix(Sknotvalues)) error("knots must be a matrix");
    ae_int_t ny = 1;  // not vector valued yet
    ae_int_t nx = nrows(Sknotvalues)-ny;
    int N = ncols(Sknotvalues);
    ae_int_t nlayers = INTEGER(AS_INTEGER(Slay))[0];
    double rbase = REAL(AS_NUMERIC(Srbase))[0];
    double lambdaN = REAL(AS_NUMERIC(lambda))[0];
    try {
      rbfmodel *s = new rbfmodel;
      rbfcreate(nx, ny, *s);
      real_2d_array dataset;
      dataset.attach_to_ptr(N, nx+ny, REAL(Sknotvalues));
      rbfsetpoints(*s, dataset, N);
      rbfsetalgohierarchical(*s, rbase, nlayers, lambdaN);
      rbfreport rep;
      rbfbuildmodel(*s, rep);
      SEXP result = R_MakeExternalPtr(s, NULL, Sknotvalues);
      
      R_RegisterCFinalizerEx(result, mfinal, FALSE);
      // return some structure which agrees with R's gc
      return result;
    } catch(...) {error("exception from alglib");}
    }
    
    SEXP R_evalrbf(SEXP model, SEXP vectors, SEXP Sthreads) {
      rbfmodel *s = (rbfmodel*) R_ExternalPtrAddr(model);
    const int M = nrows(vectors);
    const int N = ncols(vectors);
    SEXP res = PROTECT(NEW_NUMERIC(N));
    double *out = REAL(res);
    double *vec = REAL(vectors);
    int threads = INTEGER(AS_INTEGER(Sthreads))[0];
    // No OpenMP yet, awkward interface in alglib
    try {
      for(int i = 0; i < N; i++) {
	real_1d_array x, y;
	x.attach_to_ptr(M, vec+i*M);
	y.attach_to_ptr(M, out+i);
	rbfcalcbuf(*s, x, y);
      }
    } catch(...) {error("exception from alglib");}
    UNPROTECT(1);
    return res;
  }
}
#endif
