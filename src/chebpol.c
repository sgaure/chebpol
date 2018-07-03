#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Visibility.h>
#include "config.h"
#include "chebpol.h"
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static void chebcoef(double *x, int *dims, int *dimlen, double *F, int dct) {
  int siz = 1;
  const int rank = *dimlen;
  for(int i = 0; i < rank; i++) siz *= dims[i];

#ifdef HAVE_FFTW
  double isiz = 1.0/siz;
  int rdims[rank];
  /* Create a plan */
  fftw_r2r_kind kind[rank];
  for(int i = 0; i < rank; i++) kind[i] = FFTW_REDFT10;  // type II DCT

  // reverse the dimensions. fftw uses row-major order. 
  for(int i = 0; i < rank; i++) rdims[i] = dims[rank-i-1];

  // Plan and execute
  
  fftw_plan plan = fftw_plan_r2r(rank, rdims, x, F, kind, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
  if(plan == NULL) error("FFTW can't create execution plan for transform");
  fftw_execute(plan);
  // A new plan with the same parameters is fast to create it says, so we destroy this to
  // clean up memory
  fftw_destroy_plan(plan);
  // adjust scale to fit our setup
  if(!dct) for(int i = 0; i < siz; i++) F[i] *= isiz;
#else
  double *src, *dest, *buf;
  double **mat;
  double beta=0;

  // Some work space
  // 
  buf = (double*) R_alloc(siz,sizeof(double));
  mat = (double**) R_alloc(rank,sizeof(double*));
  // Create the needed transform matrices
  // reuse same dimension matrices
  for(int j = 0; j < rank; j++) {
    // Is it there already?
    int N = dims[j];
    double *jmat = NULL;
    for(int k = 0; k < j; k++) {
      if(dims[k] == N) {
	jmat = mat[k];
	break;
      }
    }
    if(jmat != NULL) {
      mat[j] = jmat;
      continue;
    }
    jmat = mat[j] = (double*) R_alloc(N*N,sizeof(double));
    for(int k = 0; k < N; k++) {
      double *jkvec = &jmat[k*N];
      for(int i = 0; i < N; i++) {
	jkvec[i] = cospi(k*(i+0.5)/N);
      }
    }
  }
  // We will switch src and dest between F and buf
  // We should end up with dest=F, so if 
  // rank is odd we should start with src = F, dest=buf
  if((rank & 1) == 1) {
    src = F;
    dest = buf;
  } else {
    src = buf;
    dest = F;
  }
  memcpy(dest,x,siz*sizeof(double));
  for(int i = rank-1; i >= 0; i--) {
    // transformation of dimension i, put in front
    int N = dims[i];  // length of transform
    int stride = siz/N;
    double alpha = dct ? 2.0 : 2.0/N;  // Constant for transform
    // swap src and dest
    double *p = src;
    src = dest;
    dest = p;

    F77_CALL(dgemm)("TRANS","TRANS",&N,&stride,&N,&alpha,mat[i],&N,src,&stride,&beta,dest,&N);
    // Fix the first element in each vector, nah do it at the end
    //    for(int k = 0; k < siz; k+= N) dest[k] *= 0.5;
  }
#endif
  if(!dct) {
    // We need to adjust the first element in each dimension by 0.5
    // The DCT-II isn't exactly equal to the Chebyshev transform
    int blocklen = 1;  // length of block
    int stride = 1; // distance between blocks
    for(int i = 0; i < rank; i++) {
      stride *= dims[i];
      for(int j = 0; j < siz; j += stride) {
	for(int s = 0; s < blocklen; s++) {
	  F[j+s] *= 0.5;
	}
      }
      blocklen = stride;
    }
  }
}



static SEXP R_chebcoef(SEXP x, SEXP sdct) {

  SEXP dim;
  int rank;
  int *dims;
  int siz;
  SEXP resvec;
  int sdims;
  int dct;
  if(!isLogical(sdct) || LENGTH(sdct) < 1) error("dct must be a logical");
  dct = LOGICAL(sdct)[0];
  dim = getAttrib(x,R_DimSymbol);
  // If no dim-attribute, assume one-dimensional
  if(isNull(dim)) {
    sdims = LENGTH(x);
    dims = &sdims;
    rank = 1;

  } else {
    rank = LENGTH(dim);
    if(rank <= 0) error("Rank must be positive");
    dims = INTEGER(dim);
  }
  siz = 1;
  for(int i = 0; i < rank; i++) siz *= dims[i];
  if(siz == 0) error("array size must be positive");

  PROTECT(resvec = NEW_NUMERIC(siz));
  chebcoef(REAL(x),dims,&rank,REAL(resvec),dct);
  setAttrib(resvec,R_DimSymbol,dim);
  setAttrib(resvec,R_DimNamesSymbol,getAttrib(x,R_DimNamesSymbol));
  UNPROTECT(1);

  return resvec;
}

static double C_FH(double *fv, double *x, double **knots, int *dims, const int rank, double **weights) {
  // Use Floater-Hormann method
  if(rank == 0) return fv[0];
  int siz = 1;
  const int newrank = rank-1;
  const int N = dims[newrank];
  const double xx = x[newrank];
  double *kn = knots[newrank];
  double *w = weights[newrank];
  for(int i = 0; i < newrank; i++) siz *= dims[i];
  double num=0, denom=0;

  // Special case:
  for(int i = 0; i < N; i++) {
    if(xx == kn[i]) return C_FH(&fv[i*siz], x, knots, dims, newrank, weights);
  }

#if 1
  // save a recursion step
  if(newrank == 0) {
    for(int i = 0; i < N; i++) {
      const double val = fv[i];
      double pole = w[i] / (xx-kn[i]);  // Should have used the weights instead of 1.0 and the sign
      //      if( (i&1) == 1) pole = -pole;
      //      if(i == 0 || i == N-1) pole = 0.5*pole;
      num += pole * val;
      denom += pole;
    }
    return num/denom;
  }
#endif 
  for(int i = 0,j=0; i < N; i++,j+=siz) {
    const double val = C_FH(&fv[j], x, knots, dims, newrank, weights);
    double pole = w[i] / (xx-kn[i]); // Should have used the weights instead of 1.0 and the sign
    //    if( (i&1) == 1) pole = -pole;
    //    if(i == 0 || i == N-1) pole = 0.5*pole;
    num += pole * val;
    denom += pole;
  }
  return num/denom;
}

static SEXP R_FH(SEXP inx, SEXP vals, SEXP grid, SEXP Sweights, SEXP Rthreads) {
  int *dims;
  int siz = 1;
  double *val = REAL(vals);
  const int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  SEXP dim;
  int rank;

  // Create some pointers and stuff. 
  dim = getAttrib(vals,R_DimSymbol);
  dims = INTEGER(dim);
  rank = LENGTH(dim);
  if(rank <= 0) error("rank must be positive");
  if(isMatrix(inx)) {
    if(rank != nrows(inx))
      error("coeffcient rank(%d) must match number of rows(%d)",LENGTH(dim),nrows(inx));
  } else {
    if(rank != LENGTH(inx))
      error("coefficient rank(%d) does not match argument length(%d)",
	    LENGTH(dim),LENGTH(inx));
  }
  /* This shouldn't happen, but we check it anyway since we'll bomb if it's wrong */
  for(int i = 0; i < rank; i++)  siz *= dims[i];
  if(LENGTH(vals) != siz)
    error("coefficient length(%d) does not match data length(%d)",
	  LENGTH(vals),siz);

  if(LENGTH(grid) != rank)
    error("There must be one value for each knot");
  double **knots = (double **) R_alloc(rank,sizeof(double*));
  double **weights = (double **) R_alloc(rank, sizeof(double*));
  for(int i = 0; i < rank; i++) {
    knots[i] = REAL(VECTOR_ELT(grid,i));
    weights[i] = REAL(VECTOR_ELT(Sweights,i));
  }
  double *xp = REAL(inx);
  const int numvec = isMatrix(inx) ? ncols(inx) : 1;
#ifdef RETMAT
  SEXP resvec = PROTECT(allocMatrix(REALSXP, numvec, 1));
#else
  SEXP resvec = PROTECT(NEW_NUMERIC(numvec));
#endif
  double *out = REAL(resvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for(int i = 0; i < numvec; i++) {
    out[i] = C_FH(val, xp+i*rank, knots, dims, rank, weights);
  }
  UNPROTECT(1);
  return resvec;
}

static double C_evalcheb(double *cf, double *x, int *dims, const int rank) {
  if(rank == 0) return cf[0];
  // Otherwise, use the Clenshaw algorithm
  int siz = 1;
  const int newrank = rank-1;
  const int N = dims[newrank];
  double x2 = 2.0*x[newrank], bn1=0, bn2=0, bn=0;

#if 1
  // semantically unnecessary, just save a recursion level, it speeds up things.
  if(newrank == 0) {
    for(int i = N-1; i > 0; i--) {
      bn2 = bn1; bn1 = bn;
      bn = x2*bn1 - bn2 + cf[i];
    }
    return x[0]*bn - bn1 + cf[0];
  }
#endif
  for(int i= 0; i < newrank; i++) siz *= dims[i];
  for(int i = N-1,j=siz*(N-1); i >= 0; i--,j-=siz) {
    bn2 = bn1; bn1 = bn;
    bn = x2*bn1 - bn2 + C_evalcheb(&cf[j],x,dims,newrank);
  }
  return bn - x[newrank]*bn1;
}

static SEXP R_evalcheb(SEXP coef, SEXP inx, SEXP Rthreads) {
  int *dims;
  int siz = 1;
  double *cf = REAL(coef);
  const int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  SEXP dim;
  int rank;

  // Create some pointers and stuff. 
  dim = getAttrib(coef,R_DimSymbol);
  dims = INTEGER(dim);
  rank = LENGTH(dim);
  if(rank <= 0) error("rank must be positive");
  if(isMatrix(inx)) {
    if(rank != nrows(inx))
      error("coeffcient rank(%d) must match number of rows(%d)",LENGTH(dim),nrows(inx));
  } else {
    if(rank != LENGTH(inx))
      error("coefficient rank(%d) does not match argument length(%d)",
	    LENGTH(dim),LENGTH(inx));
  }
  /* This shouldn't happen, but we check it anyway since we'll bomb if it's wrong */
  for(int i = 0; i < rank; i++)  siz *= dims[i];
  if(LENGTH(coef) != siz)
    error("coefficient length(%d) does not match data length(%d)",
	  LENGTH(coef),siz);

  double *xp = REAL(inx);
  const int numvec = isMatrix(inx) ? ncols(inx) : 1;
#ifdef RETMAT
  SEXP resvec = PROTECT(allocMatrix(REALSXP, numvec, 1));
#else
  SEXP resvec = PROTECT(NEW_NUMERIC(numvec));
#endif
  double *out = REAL(resvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalcheb(cf, xp+i*rank, dims, rank);
  }
  UNPROTECT(1);
  return resvec;
}

// an R-free evalongrid. Recursive
static void C_evalongrid(void (*fun)(double *x, double *y, int valuedim, void *ud),
		double *arg, double **grid,
		const int *dims, const int rank, const int valuedim, double *result, void *userdata) {
  int mrank = rank-1;
  int stride = valuedim;

  if(rank == 0) {
    fun(arg,&result[0],valuedim,userdata);
    return;
  }

  for(int i = 0; i < mrank; i++) stride *= dims[i];
  for(int i = 0,j = 0; i < dims[mrank]; i++, j += stride) {
    arg[mrank] = grid[mrank][i];
    C_evalongrid(fun, arg, grid, dims, mrank, valuedim, &result[j], userdata);
  }
}

static void C_call(double *x, double *y, const int valuedim, void *userdata) {
    // don't need x, because the arg-pointer which is used is set in the
    // R_fcall structure
  if(x == NULL) {}; // avoid warning
  double *fv = REAL(eval( *(SEXP*)userdata, R_BaseEnv));
  for(int i = 0; i < valuedim; i++) y[i] = fv[i];
}


static SEXP R_evalongrid(SEXP fun, SEXP sgrid) {
  int rank = LENGTH(sgrid);
  double *grid[rank];
  int *dims;
  R_xlen_t len=1.0;
  SEXP resvec;
  SEXP R_fcall;
  SEXP R_arg;
  SEXP R_dim;
  int valuedim;

  if(!isFunction(fun)) error("'fun' must be a function");
  PROTECT(R_dim = NEW_INTEGER(rank));
  dims = INTEGER(R_dim);
  for(int i = 0; i < rank; i++) {
    grid[i] = REAL(VECTOR_ELT(sgrid,i));
    dims[i] = LENGTH(VECTOR_ELT(sgrid,i));
    len *= dims[i];
  }

  // Create a call to fun
  PROTECT(R_fcall = lang2(fun, R_NilValue));
  // and an argument list to it
  PROTECT(R_arg = NEW_NUMERIC(rank));
  SETCADR(R_fcall, R_arg);

  // find the dimension of the result
  for(int i = 0; i < rank; i++) REAL(R_arg)[i] = grid[i][0];
  SEXP fv = eval(R_fcall,R_BaseEnv);
  if(!IS_NUMERIC(fv)) error("fun must return a real value");
  valuedim = LENGTH(fv);
  
  // Now, it's just to fill in arg, and 
  // do REAL(eval(R_fcall,NULL))[0] 
  PROTECT(resvec = NEW_NUMERIC(len*valuedim));
  C_evalongrid(C_call,REAL(R_arg),grid,dims,rank,valuedim,REAL(resvec), (void*) &R_fcall);

  // Set the dim-attribute
  if(valuedim == 1) {
    setAttrib(resvec,R_DimSymbol,R_dim);
  } else {
    SEXP snewdim;
    int *newdim;
    PROTECT(snewdim = NEW_INTEGER(rank+1));
    newdim = INTEGER(snewdim);
    newdim[0] = valuedim;
    for(int i = 1; i < rank+1; i++)
      newdim[i] = dims[i-1];
    setAttrib(resvec,R_DimSymbol,snewdim);
    UNPROTECT(1);
  }
  UNPROTECT(4);
  return resvec;
}


// return index of the smallest element larger or equal to val
// binary search is slower than linear for small n, but we don't optimize that here
// small grids anyway run faster
// we never return 0, if val < arr[0] we return 1
// if val > arr[n-1], we return n-1
// this is suited for our particular application
static inline int binsearch(const int n,double *arr, const double val) {
  int toosmall, toolarge, test;
  toosmall = 0;
  toolarge = n-1;
  while(toosmall+1 < toolarge) {
    test = (toosmall + toolarge) >> 1;
    if(arr[test] >= val) {
      toolarge = test;
    } else if(arr[test] < val) {
      toosmall = test;
    } 
  }
  return toolarge;
}

static void C_predmlip(const int rank, double **grid, int *dims, double *values, double *output) {
  int stride[rank+1];
  stride[0] = 1;
  for(int i = 1; i <= rank; ++i) stride[i] = dims[i-1]*stride[i-1];
  int N = stride[rank];
  int dim[rank];
  for(int i = 0; i < rank; ++i) {
    dim[i] = 0;
  }
  /* Loop over entire grid. For each point, predict it as the mean of all its immediate neighbours. */
  for(int j = 0; j < N; j++,dim[0]++) {
    for(int i = 0; i < rank-1; i++) {
      if(dim[i] >= dims[i]) {
	dim[i] -= dims[i];
	dim[i+1]++;
      } else 
	break;
    }
    //    printf("j %d: ",j);for(int i = 0; i < rank; i++) printf("%d ",dim[i]); printf("\n");
    double pred = 0;
    int w = 0;
    for(int i = 0; i < rank; i++) {
      if(dim[i] <= 0 || dim[i]+1 >= dims[i]) continue;
      double *gr = grid[i];
      double r = (gr[dim[i]] - gr[dim[i]-1])/(gr[dim[i]+1] - gr[dim[i]-1]);
      w++;
      pred += r*values[j+stride[i]] + (1-r)*values[j-stride[i]];
    }
    if(w == 0)
      output[j] = values[j];
    else
      output[j] = pred/w;
  }
}

static SEXP R_mlippred(SEXP sgrid, SEXP values) {
  int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  double *grid[rank];
  SEXP resvec;
  if(!IS_NUMERIC(values)) error("values must be numeric");

  /* Make some pointers into the grid data */
  for(int i = 0; i < rank; i++) {
    dims[i] = LENGTH(VECTOR_ELT(sgrid,i));
    grid[i] = REAL(VECTOR_ELT(sgrid,i));
    gridsize *= dims[i];
  }

  if(LENGTH(values) != gridsize) error("grid has size %d, you supplied %d values",
				       gridsize,LENGTH(values));

  PROTECT(resvec = NEW_NUMERIC(gridsize));
  (void) C_predmlip(rank,grid,dims,REAL(values),REAL(resvec));
  UNPROTECT(1);
  return resvec;
}

static double C_evalmlip(const int rank, double *x, double **grid, int *dims, double *values) {

  double weight[rank];
  int valpos = 0;
  double ipval = 0;
  /*
   Find first grid point which is larger or equal, compute
   the weight, store it. As well as the linear position valpos in the grid
   A binary search is faster if there are more grid points
   but up to 10-20 linear search is usually equally fast or faster due to simplicity
   Note that we never check whether x is within the grid. If it's outside we assume
   the function continues linearly outwards.
   Remaining optimizations:  We could have precomputed the stuff we divide by, if we should
   do more than one vector in a row, this would save some time.  Likewise the stride (which is
   integer multiplication, which is slow.)
  */

  int stride = 1;
  for(int i = 0; i < rank; i++) {
    double *gdvec = grid[i];
    // We could use Rs findInterval, but how many ways are there to do this?
    int gp = binsearch(dims[i],gdvec,x[i]);
    weight[i] = (x[i]-gdvec[gp-1])/(gdvec[gp]-gdvec[gp-1]);
    valpos += stride*gp;
    stride *= dims[i];
  }

  // loop over the corners of the box, sum values with weights
  for(int i = 0; i < (1<<rank); i++) {
    // i represents a corner. bit=1 if upper corner, 0 if lower corner.
    // We should find its weight
    // it's a product of the weights from each dimension
    int vpos = valpos;
    int stride = 1;
    double cw = 1;
    for(int g = 0; g < rank; g++) {
      if( (1<<g) & i) {
	cw *= weight[g];
      } else {
	cw *= 1.0-weight[g];
	vpos -= stride;
      }
      stride *= dims[g];
    }
    ipval += cw*values[vpos];
  }
  return ipval;
}

/* Then a multilinear approximation */
static SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x, SEXP Rthreads) {
  const int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  double *grid[rank];

  if(!IS_NUMERIC(values)) error("values must be numeric");
  if(!IS_NUMERIC(x)) error("argument x must be numeric");  
  if(isMatrix(x) ? (nrows(x) != rank) : (LENGTH(x) != rank))
    error("grid has dimension %d, you supplied a length %d vector", 
	  rank, isMatrix(x) ? nrows(x) : LENGTH(x));

  /* Make some pointers into the grid data */
  for(int i = 0; i < rank; i++) {
    dims[i] = LENGTH(VECTOR_ELT(sgrid,i));
    grid[i] = REAL(VECTOR_ELT(sgrid,i));
    gridsize *= dims[i];
  }

  if(LENGTH(values) != gridsize) error("grid has size %d, you supplied %d values",
				       gridsize,LENGTH(values));
  const int numvec = isMatrix(x) ? ncols(x) : 1;
  double *xp = REAL(x);
#ifdef RETMAT
  SEXP resvec = PROTECT(allocMatrix(REALSXP, numvec, 1));
#else
  SEXP resvec = PROTECT(NEW_NUMERIC(numvec));
#endif
  double *out = REAL(resvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalmlip(rank,xp+i*rank,grid,dims,REAL(values));
  }
  UNPROTECT(1);
  return resvec;
}

static SEXP R_sqdiffs(SEXP x1, SEXP x2) {
  // each column in x1 should be subtracted from each column in x2,
  // the squared column sums should be returned.
  int r1 = nrows(x1), c1 = ncols(x1), r2 = nrows(x2), c2 = ncols(x2);
  int N = c1*c2;
  SEXP res = PROTECT(NEW_NUMERIC(N));
  double *dres = REAL(res);
  double *np = dres;
  for(int i = 0; i < c1; i++) {
    double *x1p = REAL(x1) + i*r1;
    for(int j = 0; j < c2; j++) {
      double *x2p = REAL(x2) + j*r2;
      double csum = 0.0;
      for(int k = 0; k < r1; k++) {
	csum += (x1p[k] - x2p[k])*(x1p[k] - x2p[k]);
      }
      *np++ = csum;
    }
  }
  SEXP snewdim = PROTECT(NEW_INTEGER(2));
  INTEGER(snewdim)[0] = c1;
  INTEGER(snewdim)[1] = c2;
  setAttrib(res,R_DimSymbol,snewdim);
  UNPROTECT(2);
  return res;
}

static SEXP R_havefftw() {
  SEXP res;
  PROTECT(res = NEW_LOGICAL(1));
#ifdef HAVE_FFTW
  LOGICAL(res)[0] = TRUE;
#else
  LOGICAL(res)[0] = FALSE;
#endif
  UNPROTECT(1);
  return res;
}

R_CallMethodDef callMethods[] = {
  {"evalcheb", (DL_FUNC) &R_evalcheb, 3},
  {"chebcoef", (DL_FUNC) &R_chebcoef, 2},
  {"FH", (DL_FUNC) &R_FH, 5},
  {"evalmlip", (DL_FUNC) &R_evalmlip, 4},
  //  {"predmlip", (DL_FUNC) &R_mlippred, 2},
  {"evalongrid", (DL_FUNC) &R_evalongrid, 2},
  {"havefftw", (DL_FUNC) &R_havefftw, 0},
  {"sqdiffs", (DL_FUNC) &R_sqdiffs, 2},
#ifdef HAVE_ALGLIB
  {"makerbf", (DL_FUNC) &R_makerbf, 4},
  {"evalrbf", (DL_FUNC) &R_evalrbf, 3},
#endif
  {NULL, NULL, 0}
};


void attribute_visible R_init_chebpol(DllInfo *info) {
  if(info != NULL) {}; // avoid warning about unused parameter
  /* register our routines */
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info, FALSE);
  R_RegisterCCallable("chebpol", "evalcheb", (DL_FUNC) C_evalcheb);
  R_RegisterCCallable("chebpol", "evalmlip", (DL_FUNC) C_evalmlip);
  R_RegisterCCallable("chebpol", "evalongrid", (DL_FUNC) C_evalongrid);
}
#ifdef HAVE_FFTW
void attribute_visible R_unload_chebpol(DllInfo *info) {
  if(info == NULL) {};
  // Clean up fftw
  fftw_cleanup();
}
#endif
