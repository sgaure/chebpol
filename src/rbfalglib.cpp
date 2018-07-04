#include "config.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#ifdef _WIN32
#undef HAVE_ALGLIB
#endif
#ifdef HAVE_ALGLIB
#include "stdafx.h"
#include "interpolation.h"
#include <exception>
#ifdef _OPENMP
#include <omp.h>
int omp_get_thread_num();
#endif
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
    bool useomp = false;
#ifdef _OPENMP
    useomp = N > 1 && threads > 1;
#endif
    if(useomp) {
      bool *init = new bool[threads];
      rbfcalcbuffer *bufs = new rbfcalcbuffer[threads];
      for(int t = 0; t < threads; t++) init[t] = false;
      try {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
	for(int i = 0; i < N; i++) {
	  int thr = omp_get_thread_num();
	  real_1d_array x,y;
	  if(!init[thr]) {
	    init[thr] = true;
	    rbfcreatecalcbuffer(*s, bufs[thr]);
	  }
	  x.attach_to_ptr(M, vec+i*M);
	  y.attach_to_ptr(M, out+i);
	  rbftscalcbuf(*s, bufs[thr], x, y);
	}
      } catch(...) {error("exception from alglib/omp");}
      delete [] init;
      delete [] bufs;
    } else {
      try {
	for(int i = 0; i < N; i++) {
	  real_1d_array x, y;
	  x.attach_to_ptr(M, vec+i*M);
	  y.attach_to_ptr(M, out+i);
	  rbfcalcbuf(*s, x, y);
	}
      } catch(...) {error("exception from alglib");}
    }
    UNPROTECT(1);
    return res;
  }
  SEXP R_havealglib() {
    return ScalarLogical(1);
  }
}
#else
extern "C" {
  SEXP R_makerbf(SEXP a, SEXP b,  SEXP c, SEXP d) {
    if(a == NULL || b == NULL || c == NULL || d == NULL) {};
    error("alglib not supported");
  }
  SEXP R_evalrbf(SEXP a, SEXP b,  SEXP c) {
    if(a == NULL || b == NULL || c == NULL) {};
    error("alglib not supported");
  }
  SEXP R_havealglib() {
    return ScalarLogical(0);
  }

}
#endif
