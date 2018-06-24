#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Visibility.h>
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
	jkvec[i] = cos(M_PI*k*(i+0.5)/N);
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


/*
  The xvec is an array of pointers to vectors
  Vector i is of length d=dims[i] and contains
  the 1-dimensional Chebyshev vector on coordinate i of x
  I.e. the d-vector cos((0:(d-1))*acos((x[i]))) where
  x has been interval transformed into [-1,1]
*/

static double evalcheb(double *cf, double *xvec[], int *dims, const int rank) {
  double res = 0.0;
  double *xl;
  int siz = 1;
  const int newrank = rank-1;
  const int N = dims[newrank];

  if(newrank == 0) {
    for(int i = 0; i < N; i++) res += cf[i]*xvec[0][i];
    return res;
  }
  // loop over last dimension, recurse
  // integer multiplication is slow, but it seems we have to
  for(int i = 0; i < newrank; i++) siz *= dims[i];
  xl = xvec[newrank];
  for(int i = 0,j=0; i < N; i++,j+=siz) {
    res += xl[i]*evalcheb(&cf[j],xvec,dims,newrank);
  }
  return res;
}

static double C_evalcheb(double *cf, double *x, int *dims, const int rank) {
  double *xvec[rank];
  double res;

  if(rank == 1) {
    double x0 = x[0];
    // special case for rank=1. Clenshaw.
    // Could we generalize this so as not to allocate extra memory for rank > 1 too?
    double bn1=0, bn2 = 0, bn=0;
    for(int i = dims[0]-1; i >= 0; i--) {
      bn2 = bn1; bn1 = bn;
      bn = 2*x0*bn1 - bn2 + cf[i];
    }
    return bn - x0*bn1;
  }

  // make cos(j*acos(x)) for each dimension. 
  // Use recurrence relation (cos(nx) = 2cos(x)cos((n-1)x) - cos((n-2)x))
  // I.e. the Chebyshev identity

  for(int i = 0; i < rank; i++) {
    double *xv = Calloc(dims[i],double);
    xvec[i] = xv;
    xv[0] = 1.0;
    if(dims[i] == 1) break;
    xv[1] = x[i];
    for(int j = 2; j < dims[i]; j++) xv[j] = 2*x[i]*xv[j-1] - xv[j-2];
  }

  res = evalcheb(cf,xvec,dims,rank);
  for(int i = 0; i < rank; i++) Free(xvec[i]);
  return res;
}

static SEXP R_evalcheb(SEXP coef, SEXP inx) {
  int *dims;
  int siz = 1;
  double *cf = REAL(coef);
  SEXP resvec;
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
  /* This shouldn't happen, but we check it anyway since we're bomb if it's wrong */
  for(int i = 0; i < rank; i++)  siz *= dims[i];
  if(LENGTH(coef) != siz)
    error("coefficient length(%d) does not match data length(%d)",
	  LENGTH(coef),siz);

  double *xp = REAL(inx);
  int numvec = isMatrix(inx) ? ncols(inx) : 1;
  PROTECT(resvec = NEW_NUMERIC(numvec));
  for(int i = 0; i < ncols(inx); i++) {
    REAL(resvec)[i] = C_evalcheb(cf, xp, dims, rank);
    xp += rank;
  }
  UNPROTECT(1);
  return resvec;
}

// an R-free evalongrid. Recursive
void C_evalongrid(void (*fun)(double *x, double *y, int valuedim, void *ud),
		double *arg, double **grid,
		const int *dims, const int rank, const int valuedim, double *result, void *userdata) {
  int mrank = rank-1;
  int stride = valuedim;

  if(mrank == 0) {
    /* we could have terminated on rank==0 with a single result[0] = fun(arg),
       but we stop earlier to save some function calls */
    for(int i = 0, j=0; i < dims[0]; i++, j+=valuedim) {
      arg[0] = grid[0][i];
      fun(arg,&result[j],valuedim,userdata);
      //      result[i] = fun(arg);
    }
    return;
  }
  for(int i = 0; i < mrank; i++) stride *= dims[i];

  for(int i = 0,j = 0; i < dims[mrank]; i++, j += stride) {
    arg[mrank] = grid[mrank][i];
    C_evalongrid(fun, arg, grid, dims, mrank, valuedim, &result[j], userdata);
  }
}

void C_call(double *x, double *y, const int valuedim, void *userdata) {
    // don't need x, because the arg-pointer which is used is set in the
    // R_fcall structure

  double *fv = REAL(eval( *(SEXP*)userdata, R_BaseEnv));
  for(int i = 0; i < valuedim; i++) y[i] = fv[i];
}


static SEXP R_evalongrid(SEXP fun, SEXP sgrid) {
  int rank = LENGTH(sgrid);
  double *grid[rank];
  int *dims;
  double len=1.0;
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
static SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x) {
  const int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  double *grid[rank];
  SEXP resvec;
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
  PROTECT(resvec = NEW_NUMERIC(numvec));
  for(int i = 0; i < numvec; i++) {
    REAL(resvec)[i] = C_evalmlip(rank,xp,grid,dims,REAL(values));
    xp += rank;
  }
  UNPROTECT(1);
  return resvec;
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
  {"evalcheb", (DL_FUNC) &R_evalcheb, 2},
  {"chebcoef", (DL_FUNC) &R_chebcoef, 2},
  {"evalmlip", (DL_FUNC) &R_evalmlip, 3},
  {"predmlip", (DL_FUNC) &R_mlippred, 2},
  {"evalongrid", (DL_FUNC) &R_evalongrid, 2},
  {"havefftw", (DL_FUNC) &R_havefftw, 0},
  {NULL, NULL, 0}
};


void attribute_visible R_init_chebpol(DllInfo *info) {
  /* register our routines */
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info, FALSE);
  R_RegisterCCallable("chebpol", "evalcheb", (DL_FUNC) C_evalcheb);
  R_RegisterCCallable("chebpol", "evalmlip", (DL_FUNC) C_evalmlip);
  R_RegisterCCallable("chebpol", "evalongrid", (DL_FUNC) C_evalongrid);

}
#ifdef HAVE_FFTW
void attribute_visible R_unload_chebpol(DllInfo *info) {
  // Clean up fftw
  fftw_cleanup();
}
#endif
