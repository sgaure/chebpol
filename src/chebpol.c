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
      
      F77_CALL(dgemm)("T", "T", &N, &stride, &N, &alpha, mat[i], &N, src, &stride, &beta, dest, &N FCONE FCONE);
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
    if(fabs(xx - kn[i]) < 10.0*DBL_EPSILON) return FH(&fv[i*siz], x, knots, dims, newrank, weights);
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
  SEXP resvec = PROTECT(NEW_NUMERIC(numvec));
  double *out = REAL(resvec);
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = FH(val, xp+i*rank, knots, dims, rank, weights);
  }
  UNPROTECT(1);
  return resvec;
}

R_INLINE static void fhweights(const int n, const double *grid, const int d, double *w, int threads) {
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
  for(int k = 0; k <= n; k++) {
    const int start = (k < d) ? 0 : k-d, end = (k < n-d) ? k : n-d;
    double sum = 0.0;
    for(int i = start; i <= end; i++) {
      double prod = 1.0;
      for(int j = i; j <= i+d; j++) {
	if(j==k) continue;
	prod *= grid[k]-grid[j];
      }
      sum += 1.0/fabs(prod);
    }
    w[k] = sum * ( ( k % 2 != d % 2) ? -1.0 : 1.0);
  }
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
    // The paper use index from 0 to n, i.e. n+1 knots in each dimension. Therefore dims[r]-1
    fhweights(dims[r]-1, grid[r], dd[r], wlist[r],threads);
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

static double C_evalmlip(const int rank, double *x, double **grid, int *dims, 
			 double *values, int blend) {

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
  double wsum = 0.0;
  for(int i = 0; i < (1<<rank); i++) {
    // i represents a corner. bit=1 if upper corner, 0 if lower corner.
    // We should find its weight
    // it's a product of the weights from each dimension
    int vpos = valpos;
    int stride = 1;
    double cw = 1;
    for(int g = 0; g < rank; g++) {
      if( (1<<g) & i) {
	cw *= blendfun(weight[g],blend);
      } else {
	cw *= blendfun(1-weight[g],blend);
	vpos -= stride;
      }
      stride *= dims[g];
    }
    wsum += cw;
    ipval += cw*values[vpos];
  }
  return ipval/wsum;
}

/* Then a multilinear approximation */
static SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x, SEXP Rthreads, SEXP Sblend) {
  const int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  double *grid[rank];
  int blend = 0;
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
  if(!isNull(Sblend)) blend = INTEGER(AS_INTEGER(Sblend))[0];
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
    out[i] = C_evalmlip(rank,xp+i*rank,grid,dims,REAL(values),blend);
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

static double findsimplex(double *x, double *knots, int *dtri, double *lumats, int *ipivs, double *bbox, int epol,
			  double *val, const int dim, const int numsimplex, int nknots, int blend) {

  double vec[dim+1];
  int N=dim+1, one=1, info;
  for(int simplex = 0; simplex < numsimplex; simplex++) {
    // x must be in the bounding box. Perhaps we should do binary search, or organize them hierarchically?
    double *box = bbox + simplex*2*dim;
    int bad = 0;
    for(int i = 0; i < dim; i++) {
      if(x[i] < box[0] || x[i] > box[1]) {bad = 1; break;}
      box += 2;
    }
    if(bad) continue;

    // bounding box matches. Transform to barycentric coordinates, check that they are positive.
    double *lumat = lumats + simplex*N*N;
    int *ipiv = ipivs + simplex*N;
    for(int d = 0; d < dim; d++) vec[d] = x[d];
    vec[dim] = 1.0;
    F77_CALL(dgetrs)("N", &N, &one, lumat, &N, ipiv, vec, &N, &info FCONE);
    for(int d = 0; d < N; d++) if(vec[d] < -1e-10) {bad = 1; break;}
    if(bad) continue;

    // We found it
    const int *tri = dtri + simplex * N;
    double sum = 0;
    double wsum = 0;
    for(int d = 0; d < N; d++) {
      const double w = blendfun(vec[d],blend);
      wsum += w;
      sum += w*val[tri[d]-1];
    }
    return sum/wsum;
  }
  
  if(epol) {
    // Find the closest knot
    double mindist = 1e99;
    int nearknot = 0;
    for(int k = 0; k < nknots; k++) {
      double dist = 0.0;
      for(int i = 0; i < dim; i++) dist += (x[i] - knots[k*dim + i])*(x[i] - knots[k*dim + i]);
      if(dist < mindist) {mindist = dist; nearknot = k;}
    }
    
    // Now, find a simplex with this knot
    int nearest = -1;
    for(int simplex = 0; simplex < numsimplex; simplex++) {
      for(int j = 0; j < N; j++) if(dtri[simplex*N + j] == nearknot+1) {nearest = simplex; break;}
      if(nearest != -1) break;
    }

    if(nearest == -1) return NA_REAL;
    double *lumat = lumats + nearest*N*N;
    int *ipiv = ipivs + nearest*N;
    const int *tri = dtri + nearest * (dim+1);
    for(int d = 0; d < dim; d++) vec[d] = x[d];
    vec[dim] = 1.0;
    F77_CALL(dgetrs)("N", &N, &one, lumat, &N, ipiv, vec, &N, &info FCONE);
    double sum = 0;
    for(int d = 0; d <= dim; d++) sum += vec[d]*val[tri[d]-1];
    return sum;
  }
  return NA_REAL;
}

static SEXP R_evalsl(SEXP Sx, SEXP Sknots, SEXP Sdtri, SEXP adata, SEXP Sval,
		     SEXP extrapolate, SEXP Sthreads, SEXP Sblend) {
  // Loop over the triangulation
  int *dtri = INTEGER(Sdtri);
  const int numsimplex = ncols(Sdtri);
  double *lumats = REAL(VECTOR_ELT(adata,0));
  int *pivots = INTEGER(VECTOR_ELT(adata,1));
  double *bboxmat = REAL(VECTOR_ELT(adata,2));
  double *x = REAL(Sx);
  const int dim = nrows(Sknots);
  const int numvec = ncols(Sx);
  double *knots = REAL(Sknots);
  double *val = REAL(Sval);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  SEXP ret = PROTECT(NEW_NUMERIC(ncols(Sx)));
  double *resvec = REAL(ret);
  int epol = LOGICAL(AS_LOGICAL(extrapolate))[0];
  int blend=0;
  if(!isNull(Sblend)) blend = INTEGER(AS_INTEGER(Sblend))[0];

#pragma omp parallel for num_threads(threads) schedule(guided) if(threads > 1 && numvec > 1)
  for(int i = 0; i < numvec; i++) {
    resvec[i] = findsimplex(x + i*dim, knots, dtri, lumats, pivots, bboxmat,
			    epol, val, dim, numsimplex, ncols(Sknots),blend);
  }
  UNPROTECT(1);
  return ret;
}

static SEXP R_analyzesimplex(SEXP Sdtri, SEXP Sknots, SEXP Sthreads) {
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  int *dtri = INTEGER(Sdtri);
  double *knots = REAL(Sknots);
  const int dim = nrows(Sknots);
  const int numsimplex = ncols(Sdtri);
  int N = dim+1;
  SEXP retlist = PROTECT(NEW_LIST(3));
  // Allocate everything before parallel for
  SET_VECTOR_ELT(retlist, 0, NEW_NUMERIC(numsimplex*N*N));
  SET_VECTOR_ELT(retlist, 1, NEW_INTEGER(numsimplex*N));
  SET_VECTOR_ELT(retlist, 2, allocMatrix(REALSXP, 2*dim, numsimplex));
  double *lumats = REAL(VECTOR_ELT(retlist,0));
  int *pivots = INTEGER(VECTOR_ELT(retlist,1));
  double *bboxmat = REAL(VECTOR_ELT(retlist,2));

#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1 && numsimplex > 1)
  for(int simplex = 0; simplex < numsimplex; simplex++) {
    int *tri = dtri + simplex*(dim+1);
    int *ipiv = pivots + simplex*N;
    for(int j = 0; j < N; j++) ipiv[j] = j+1;
    double *lumat = lumats + simplex*N*N;
    int info;
    double *m = bboxmat + simplex*2*dim;
    for(int j = 0; j < dim; j++) {m[2*j] = 1e100; m[2*j+1] = -1e100;}
    // Copy all vertices into mat, with a final row of ones to LU-factorize
    // Preparation for transforming into barycentric coordinates
    for(int v = 0; v <= dim; v++) {
      const double *knot = knots + (tri[v]-1)*dim;
      for(int j = 0; j < dim; j++) {
	lumat[v*(dim+1) + j] = knot[j];
	if(knot[j] < m[2*j]) m[2*j] = knot[j];
	if(knot[j] > m[2*j+1]) m[2*j+1] = knot[j];
      }
      lumat[v*(dim+1) + dim] = 1.0;
    }
    F77_CALL(dgetrf)(&N, &N, lumat, &N, ipiv, &info);
  }

  UNPROTECT(1);
  return retlist;
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
#ifndef HAVE_ALGLIB
SEXP R_makerbf(SEXP a, SEXP b, SEXP c, SEXP d) {
  UNUSED(a); UNUSED(b); UNUSED(c); UNUSED(d);
  error("ALGLIB not available");
  return R_NilValue;
}
SEXP R_evalrbf(SEXP a, SEXP b, SEXP c) {
  UNUSED(a); UNUSED(b); UNUSED(c);
  error("ALGLIB not available");
  return R_NilValue;
}
SEXP R_havealglib() {
  return ScalarLogical(FALSE);
}

#endif

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
  {"havealglib", (DL_FUNC) &R_havealglib, 0},
  {"evalpolyh", (DL_FUNC) &R_evalpolyh, 7},
  {"evalsl", (DL_FUNC) &R_evalsl, 8},
  {"evalstalker", (DL_FUNC) &R_evalstalker, 5},
  {"makehyp", (DL_FUNC) &R_makehyp, 2},
  {"evalhyp", (DL_FUNC) &R_evalhyp, 4},
  {"makestalk", (DL_FUNC) &R_makestalk, 2},
  {"evalstalk", (DL_FUNC) &R_evalstalk, 4},
  {"analyzesimplex", (DL_FUNC) &R_analyzesimplex, 3},
  {"havegsl", (DL_FUNC) &havegsl, 0},  
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
