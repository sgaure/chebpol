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
#include "chebpol.h"
#define UNUSED(x) (void)(x)
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static void chebcoef(double *x, int *dims, const int rank, double *F, int dct, int threads) {
  R_xlen_t siz = 1;
  for(int i = 0; i < rank; i++) siz *= dims[i];

#ifdef HAVE_FFTW
  int brokenfftw = 0;
  double isiz = 1.0/siz;
  int rdims[rank];
  /* Create a plan */
  fftw_r2r_kind kind[rank];
  for(int i = 0; i < rank; i++) kind[i] = FFTW_REDFT10;  // type II DCT
  
  // reverse the dimensions. fftw uses row-major order. 
  for(int i = 0; i < rank; i++) rdims[i] = dims[rank-i-1];
  
  // Plan and execute
  
  fftw_plan plan = fftw_plan_r2r(rank, rdims, x, F, kind, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
  if(plan != NULL) {
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  } else {
    // That's funny, perhaps we've been given MKL with its amputated FFTW?
    // Let's try a series of strided 1d-transforms
    for(R_xlen_t i = 0; i < siz; i++) F[i] = x[i];
    int stride = 1;
    for(int i = 0; i < rank; i++) {
      fftw_iodim dim;
      fftw_plan plan;
      fftw_r2r_kind kind = FFTW_REDFT10;
      const int tlen = dims[i];
      const int ntrans = siz/tlen;
      dim.n = tlen;
      dim.is = stride;
      dim.os = stride;
      plan = fftw_plan_guru_r2r(1, &dim, 0, NULL, F, F, &kind, FFTW_ESTIMATE);
      if(plan == NULL) {
	if(tlen > 10000) warning("FFTW fails on strided long (%d) vector, is it MKL?",tlen);
	brokenfftw = 1;
	break;
      }
      // Do them
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
      for(R_xlen_t j = 0; j < ntrans; j++) {
	div_t dv = div(j, stride);
	int offset = dv.rem + stride*tlen * dv.quot;
	fftw_execute_r2r(plan, F+offset, F+offset);
      }
      fftw_destroy_plan(plan);
      stride *= tlen;
    }
  }
  // adjust scale to fit our setup
  if(!dct) for(int i = 0; i < siz; i++) F[i] *= isiz;

  if(brokenfftw) {
#endif
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
      if(N > 10000) warning("Long vector (%d), no FFTW. This may go wrong.",N);
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
#ifdef HAVE_FFTW
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



static SEXP R_chebcoef(SEXP x, SEXP sdct, SEXP Sthreads) {

  SEXP dim;
  int rank;
  int *dims;
  int siz;
  SEXP resvec;
  int sdims;
  int dct;
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
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
  chebcoef(REAL(x),dims,rank,REAL(resvec),dct,threads);
  setAttrib(resvec,R_DimSymbol,dim);
  setAttrib(resvec,R_DimNamesSymbol,getAttrib(x,R_DimNamesSymbol));
  UNPROTECT(1);

  return resvec;
}

static double FH(double *fv, double *x, double **knots, int *dims, const int rank, double **weights) {
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
    if(fabs(xx - kn[i]) < 10.0*DOUBLE_EPS) return FH(&fv[i*siz], x, knots, dims, newrank, weights);
  }

#if 1
  // save a recursion step
  if(newrank == 0) {
    for(int i = 0; i < N; i++) {
      const double val = fv[i];
      double pole = w[i] / (xx-kn[i]);  
      num += pole * val;
      denom += pole;
    }
    return num/denom;
  }
#endif 
  for(int i = 0,j=0; i < N; i++,j+=siz) {
    const double val = FH(&fv[j], x, knots, dims, newrank, weights);
    double pole = w[i] / (xx-kn[i]);
    num += pole * val;
    denom += pole;
  }
  return num/denom;
}

static SEXP R_FH(SEXP inx, SEXP vals, SEXP grid, SEXP Sweights, SEXP Rthreads, SEXP spare) {
  int *dims;
  int siz = 1;
  double *val = REAL(vals);
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  SEXP dim;
  int rank;
  UNUSED(spare);
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
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = FH(val, xp+i*rank, knots, dims, rank, weights);
  }
  UNPROTECT(1);
  return resvec;
}

// Compute the weights for Floater-Hormann. Formula (18) of the FH-paper
static SEXP R_fhweights(SEXP Sgrid, SEXP Sd, SEXP Sthreads) {
  if(LENGTH(Sgrid) != LENGTH(Sd)) 
    error("Length of grid (%d) should equal length of d (%d)",LENGTH(Sgrid),LENGTH(Sd));

  const int rank = LENGTH(Sgrid);
  double *grid[rank];
  int *dims = (int*) R_alloc(rank,sizeof(int));
  for(int i = 0; i < rank; i++) {
    grid[i] = REAL(VECTOR_ELT(Sgrid,i));
    dims[i] = LENGTH(VECTOR_ELT(Sgrid,i));
  }
  double *wlist[rank];
  SEXP ret = PROTECT(NEW_LIST(rank));
  for(int i = 0; i < rank; i++) {
    SET_VECTOR_ELT(ret,i,NEW_NUMERIC(dims[i]));
    wlist[i] = REAL(VECTOR_ELT(ret,i));
  }

  const int *dd = INTEGER(AS_INTEGER(Sd));
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];

  for(int r = 0; r < rank; r++) {
    const double *gr = grid[r];
    double *w = wlist[r];
    const int d = dd[r];
    const int n = dims[r]-1;
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
    for(int k = 0; k <= n; k++) {
      const int start = (k < d) ? 0 : k-d, end = (k < n-d) ? k : n-d;
      double sum = 0.0;
      for(int i = start; i <= end; i++) {
	double prod = 1.0;
	for(int j = i; j <= i+d; j++) {
	  if(j==k) continue;
	  prod *= gr[k]-gr[j];
	}
	sum += 1.0/fabs(prod);
      }
      w[k] = sum * ( (abs(k-d) % 2 == 1) ? -1.0 : 1.0);
    }
  }
  UNPROTECT(1);
  return ret;
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

static SEXP R_evalcheb(SEXP coef, SEXP inx, SEXP Rthreads, SEXP spare) {
  int *dims;
  int siz = 1;
  double *cf = REAL(coef);
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  SEXP dim;
  int rank;
  UNUSED(spare);
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
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalcheb(cf, xp+i*rank, dims, rank);
  }
  UNPROTECT(1);
  return resvec;
}

double C_evalpolyh(const double *x, const double *knots, const double *weights, 
		   const double *lweights, const int rank, const int nknots, const double k) {
  // Find squared distance to each knot
  int ki = (int) k;
  double wsum = 0.0;
  if(k < 0) {
    for(int i = 0; i < nknots; i++) {
      const double *kn = &knots[i*rank];
      double sqdist = 0.0;
      for(int j = 0; j < rank; j++) sqdist += (x[j]-kn[j])*(x[j]-kn[j]);
      wsum += weights[i]*exp(k*sqdist);
    }
  } else if(ki % 2 == 1) {
    for(int i = 0; i < nknots; i++) {
      const double *kn = &knots[i*rank];
      double sqdist = 0.0;
      for(int j = 0; j < rank; j++) sqdist += (x[j]-kn[j])*(x[j]-kn[j]);
      if(sqdist == 0) continue;
      wsum += weights[i] * R_pow_di(sqrt(sqdist), ki);
    }
  } else {
    for(int i = 0; i < nknots; i++) {
      const double *kn = &knots[i*rank];
      double sqdist = 0.0;
      for(int j = 0; j < rank; j++) sqdist += (x[j]-kn[j])*(x[j]-kn[j]);
      if(sqdist == 0) continue;
      wsum += weights[i] * log(sqrt(sqdist)) * R_pow_di(sqrt(sqdist), ki);
    }
  }
  // then the linear part, constant is the first element
  wsum += lweights[0];
  for(int i = 0; i < rank; i++) wsum += lweights[i+1] * x[i];
  return wsum;
}
SEXP R_evalpolyh(SEXP inx, SEXP Sknots, SEXP weights, SEXP lweights, SEXP Sk, SEXP Sthreads, SEXP spare) {
  double *x = REAL(inx);
  double *knots = REAL(Sknots);
  double *w = REAL(weights);
  double *lw = REAL(lweights);
  double k = REAL(AS_NUMERIC(Sk))[0];
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  int numvec;
  const int rank = nrows(Sknots);
  UNUSED(spare);
  if(LENGTH(lweights) != rank+1) 
    error("linear weights (%d) should be one longer than rank(%d)",
	  LENGTH(lweights),rank);
  if(LENGTH(weights) != ncols(Sknots))
    error("number of weights (%d) should equal number of knots (%d)",
	  LENGTH(weights), ncols(Sknots));
  if(isMatrix(inx)) {
    if(nrows(inx) != nrows(Sknots))
      error("Dimension of input (%d) must match dimension of interpolant (%d)",nrows(inx),nrows(Sknots));
    numvec = ncols(inx);
  } else {
    if(length(inx) != nrows(Sknots))
      error("Dimension of input (%d) must match dimension of interpolant (%d)",length(inx),nrows(Sknots));
    numvec = 1;
  }
  SEXP res = PROTECT(NEW_NUMERIC(numvec));
  double *out = REAL(res);

#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalpolyh(x + i*rank, knots, w, lw, rank, ncols(Sknots), k);
  }
  UNPROTECT(1);
  return res;
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
static SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x, SEXP Rthreads, SEXP spare) {
  const int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  double *grid[rank];
  UNUSED(spare);
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
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalmlip(rank,xp+i*rank,grid,dims,REAL(values));
  }
  UNPROTECT(1);
  return resvec;
}

static SEXP R_sqdiffs(SEXP x1, SEXP x2, SEXP Sthreads) {
  // each column in x1 should be subtracted from each column in x2,
  // the squared column sums should be returned.
  const int r1 = nrows(x1), c1 = ncols(x1), r2 = nrows(x2), c2 = ncols(x2);
  int N = c1*c2;
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  SEXP res = PROTECT(NEW_NUMERIC(N));
  double *dres = REAL(res);
  double *np = dres;
#pragma omp parallel for num_threads(threads) schedule(static) if (c1 > 1)
  for(int i = 0; i < c1; i++) {
    double *x1p = REAL(x1) + i*r1;
    for(int j = 0; j < c2; j++) {
      double *x2p = REAL(x2) + j*r2;
      double csum = 0.0;
      for(int k = 0; k < r1; k++) {
	csum += (x1p[k] - x2p[k])*(x1p[k] - x2p[k]);
      }
      np[j + i*c2] = csum;
    }
  }
  SEXP snewdim = PROTECT(NEW_INTEGER(2));
  INTEGER(snewdim)[0] = c1;
  INTEGER(snewdim)[1] = c2;
  setAttrib(res,R_DimSymbol,snewdim);
  UNPROTECT(2);
  return res;
}

static double findsimplex(double *x, double *knots, int *dtri, SEXP Sort, double *bbox, int epol,
			    const int dim, const int numsimplex) {
  int retval = INT_MIN;

  for(int simplex = 0; simplex < numsimplex; simplex++) {
    // x must be in the bounding box
    double *box = bbox + simplex*2*dim;
    int bad = 0;
    for(int i = 0; i < dim; i++) {
      if(x[i] < box[0] || x[i] > box[1]) {bad = 1; break;}
      box += 2;
    }
    if(bad) continue;
    //    Rprintf("Check simplex %d\n",simplex+1);
    // Ok, it's in the bounding box. Is it in the simplex?
    // for each of the points in the simplex, check the
    // orthogonal vector pointing towards it. Does it have
    // positive inner product with our point?
    // Loop over the points in the simplex
    const int *tri = dtri + simplex * (dim+1); //nrows(Sdtri);
    const double *ort = REAL(VECTOR_ELT(Sort, simplex));
    for(int d = 0; d <= dim; d++) {
      // Remember tri is 1-based
      const double *ref = (d==0) ? knots+(tri[1]-1)*dim : knots+(tri[0]-1)*dim;
      const double *or = ort + d*dim;
      double ip = 0.0;
      for(int i = 0; i < dim; i++) ip += (x[i] - ref[i])*or[i];
      if(ip < -1e-9) {bad = 1; break;}
    }
    if(bad) continue;
    // We found it
    retval = simplex;
    break;
  }
  
  if(retval == INT_MIN && epol) {
    // Find the closest simplex
    double mindist = 1e99;
    int nearest = 0;
    for(int simplex = 0; simplex < numsimplex; simplex++) {
      double dist = 1e99;
      const int *tri = dtri + simplex * (dim+1);
      for(int d = 0; d <= dim; d++) {
	double *pt = knots + (tri[d]-1)*dim;
	double ddist = 0.0;
	for(int i = 0; i < dim; i++) ddist += (pt[i]-x[i])*(pt[i]-x[i]);
	if(ddist < dist) dist = ddist;
      }
      if(dist < mindist) {mindist = dist; nearest = simplex;}
    }
    retval = nearest;
  }
  return retval;
}

static SEXP R_evalsl(SEXP Sx, SEXP Sknots, SEXP Sdtri, SEXP Sort, SEXP Sbbox, SEXP Sval,
			  SEXP extrapolate, SEXP Sthreads, SEXP spare) {
  UNUSED(spare); 
  // Loop over the triangulation
  int *dtri = INTEGER(Sdtri);
  const int numsimplex = ncols(Sdtri);
  double *x = REAL(Sx);
  const int dim = nrows(Sknots);
  const int numvec = ncols(Sx);
  double *bbox = REAL(Sbbox);
  double *knots = REAL(Sknots);
  double *val = REAL(Sval);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  //  SEXP ret = PROTECT(NEW_INTEGER(ncols(Sx)));
  //  int *pret = INTEGER(ret);
  SEXP ret = PROTECT(NEW_NUMERIC(ncols(Sx)));
  double *resvec = REAL(ret);
  int epol = LOGICAL(AS_LOGICAL(extrapolate))[0];
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1 && numvec > 1)
  for(int i = 0; i < numvec; i++) {
    double *xx = x + i*dim;
    const int simplex = findsimplex(xx, knots, dtri, Sort, bbox, epol, dim, numsimplex);
    if(simplex == INT_MIN) {resvec[i] = NA_REAL; continue;}
    // Collect the vertices to solve for the barycentric coordinates of xx
    // Hmm, that should be possible to do directly with the orthogonal vectors in Sort
    // Have a look at it later.
    double mat[(dim+1)*(dim+1)], vec[dim+1];
    int *tri = dtri + simplex*(dim+1);
    for(int j = 0; j <= dim; j++) {
      double *knot = knots + (tri[j]-1)*dim;
      for(int k=0; k < dim; k++) mat[j*(dim+1)+k] = knot[k];
      mat[j*(dim+1)+dim] = 1.0;
    }
    for(int j = 0; j < dim; j++) vec[j] = xx[j];
    vec[dim] = 1.0;
    // now, solve mat X = vec
    int N = dim+1, one=1, info;
    int ipiv[N];
    F77_CALL(dgesv)(&N, &one, mat, &N, ipiv, vec, &N, &info);
    // Now, the barycentric coordinates are in vec
    double sum = 0.0;
    for(int j = 0; j <= dim; j++) {
      // which knot?
      int knot = tri[j]-1;
      sum += val[knot]*vec[j];
    }
    resvec[i] = sum;
  }
  UNPROTECT(1);
  return ret;
}


static void findortho(int *tri, double *knots, const int dim, double *ortmat) {
  // Collect the knots in a matrix of size dim x (dim+1)
  double mat[dim*dim];
  double vec[dim];
  for(int vertex = 0; vertex <= dim; vertex++) {

    // Copy all vertices except this into mat, this one into vec
    int matcol = 0;
    double *ref = (vertex == 0) ? knots + (tri[1]-1)*dim : knots + (tri[0]-1)*dim;
    for(int v = 0; v <= dim; v++) {
      double *knot = knots + (tri[v]-1)*dim;
      if(v == vertex) {
	for(int j = 0; j < dim; j++) vec[j] = knot[j] - ref[j];
      } else {
	for(int j = 0; j < dim; j++) mat[matcol + j] = knot[j] - ref[j];
	matcol += dim;
      }
    }
    // Solve the system
    // Use R's least squares solver
    int one=1, N=dim, rank;
    double work[2*dim], tol=1e-10, coef[dim], eff[dim], qraux[dim];
    int jpvt[dim];
    for(int i=0; i < dim; i++) jpvt[i] = i+1;
    double *ort = ortmat + dim*vertex;
    F77_CALL(dqrls)(mat, &N, &N, vec, &one, &tol, coef, ort, eff, &rank,
		    jpvt, qraux, work);
    // Nnormalize it:
    double sum = 0.0;
    for(int i = 0; i < dim; i++) sum += ort[i]*ort[i];
    sum = 1/sqrt(sum);
    for(int i = 0; i < dim; i++) ort[i] *= sum;
  }
}
static SEXP R_findortho(SEXP Sdtri, SEXP Sknots, SEXP Sthreads) {
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  int *dtri = INTEGER(Sdtri);
  double *knots = REAL(Sknots);
  const int dim = nrows(Sknots);
  const int numsimplex = ncols(Sdtri);
  SEXP retlist = PROTECT(NEW_LIST(numsimplex));
  for(int i = 0; i < numsimplex; i++) {
    SET_VECTOR_ELT(retlist,i,allocMatrix(REALSXP,dim,dim+1));
  }
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1 && numsimplex > 1)
  for(int simplex = 0; simplex < numsimplex; simplex++) {
    findortho(dtri + simplex*(dim+1), knots, dim, REAL(VECTOR_ELT(retlist,simplex)));
  }
  UNPROTECT(1);
  return retlist;
}

static SEXP R_findbbox(SEXP Sdtri, SEXP Sknots, SEXP Sthreads) {
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  int *dtri = INTEGER(Sdtri);
  double *knots = REAL(Sknots);
  const int dim = nrows(Sknots);
  const int numsimplex = ncols(Sdtri);
  SEXP retmat = PROTECT(allocMatrix(REALSXP, 2*dim, numsimplex));
  double *mat = REAL(retmat);
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1 && numsimplex > 1)
  for(int simplex = 0; simplex < numsimplex; simplex++) {
    int *tri = dtri + simplex*(dim+1);
    double *m = mat + simplex*2*dim;
    for(int j = 0; j < dim; j++) {m[2*j] = 1e100; m[2*j+1] = -1e100;}
    for(int i = 0; i <= dim; i++) {
      double *kn = knots + (tri[i]-1)*dim;
      for(int j = 0; j < dim; j++) {
	if(kn[j] < m[2*j]) m[2*j] = kn[j];
	if(kn[j] > m[2*j+1]) m[2*j+1] = kn[j];
      }
    }
  }
  UNPROTECT(1);
  return retmat;
}

static SEXP R_havefftw() {
#ifdef HAVE_FFTW
  return ScalarLogical(TRUE);
#else
  return ScalarLogical(FALSE);
#endif
}

// inplace to save memory, do repated multiplication instead of pow(). It's faster.
static SEXP R_phifunc(SEXP Sx, SEXP Sk, SEXP Sthreads) {
  double k = REAL(AS_NUMERIC(Sk))[0];
  double *x = REAL(Sx);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  double *y;
  SEXP res;
  if(XLENGTH(Sx) == 1 && x[0] == 0.0) return ScalarReal((k<0)?1:0);

  if(MAYBE_REFERENCED(Sx)) {
    res = PROTECT(NEW_NUMERIC(XLENGTH(Sx)));
    y = REAL(res);
  } else {
    res = Sx;
    y = x;
  }

  if(k < 0) {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
    for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) y[i] = exp(k*x[i]);
  } else {
    int ki = INTEGER(AS_INTEGER(Sk))[0];
    if(ki % 2 == 1) {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
      for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) {
	// it's the sqrt(x) to ki'th power
	if(x[i] <= 0.0) {y[i]=0.0;continue;}
	y[i] = R_pow_di(sqrt(x[i]), ki);
      }
    } else {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
      for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) {
	// it's sqrt(x) to ki'th power, multiplied by 0.5 log(x)
	if(x[i] <= 0.0) {y[i]=0.0;continue;}
	y[i] = 0.5*log(x[i]) * R_pow_di(sqrt(x[i]),ki);
      }
    }
  }
  if(y != x) UNPROTECT(1);
  return res;
}

R_CallMethodDef callMethods[] = {
  {"evalcheb", (DL_FUNC) &R_evalcheb, 4},
  {"chebcoef", (DL_FUNC) &R_chebcoef, 3},
  {"FH", (DL_FUNC) &R_FH, 6},
  {"FHweights", (DL_FUNC) &R_fhweights, 3},
  {"evalmlip", (DL_FUNC) &R_evalmlip, 5},
  {"predmlip", (DL_FUNC) &R_mlippred, 2},
  {"evalongrid", (DL_FUNC) &R_evalongrid, 2},
  {"havefftw", (DL_FUNC) &R_havefftw, 0},
  {"sqdiffs", (DL_FUNC) &R_sqdiffs, 3},
  {"phifunc", (DL_FUNC) &R_phifunc, 3},
  {"makerbf", (DL_FUNC) &R_makerbf, 4},
  {"evalrbf", (DL_FUNC) &R_evalrbf, 3},
  {"evalpolyh", (DL_FUNC) &R_evalpolyh, 7},
  {"evalsl", (DL_FUNC) &R_evalsl, 9},
  {"findortho", (DL_FUNC) &R_findortho, 3},
  {"findbbox", (DL_FUNC) &R_findbbox, 3},
  {"havealglib", (DL_FUNC) &R_havealglib, 0},
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
