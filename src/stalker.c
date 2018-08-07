#include "chebpol.h"

static double evalstalker(const double *x, const int nrank, const int *rank, double **grid,
			  const double *val, const double *b, const double *c, const double *d, 
			  const double mindeg, const double maxdeg, const double smooth) {
  int strides[nrank];
  double nx[nrank];
  const double irank = 1.0/nrank;
  int mflag;
  int lowlow = 0;
  int stride = 1;
  for(int r = 0; r < nrank; r++) {
    double *gr=grid[r];
    const int i = findInterval(gr, rank[r], x[r], TRUE, TRUE, 1, &mflag)-1;
    nx[r] = (x[r]-gr[i])/(gr[i+1]-gr[i]);
    lowlow += i*stride;
    strides[r] = stride;
    stride *= rank[r];
  }
  const double a3 = 0.5*(smooth - 1);
  const double a1 = 1-a3;

  // Now, lowlow is the index of the lower, left corner 
  double Val = 0.0;
  for(int crn = 0; crn < (1<<nrank); crn++) {
    // c is a corner. 
    int index = lowlow;
    for(int dir = 0; dir < nrank; dir++) {
      const int above = (1<<dir) & crn;
      if(above) index += strides[dir];
    }
    // index is the index into the val-array, but the b,c,d arrays
    // have nrank times as many, one for each direction
    // There are nrank directions from it
    // For each direction we find the value of the basis function
    // How do we weigh the directions?
    int idx = nrank*index;
    double value = 0.0;
    double ww = 1.0;
    for(int r = 0; r < nrank; r++) {
      const int above = (1<<r) & crn;
      const double p = above ? nx[r]-1 : nx[r];
      const double w = above ? nx[r] : 1-nx[r];
      double e = d[idx];
      if(e < mindeg) e = mindeg;
      if(e > maxdeg) e = maxdeg;
      if(smooth != 1.0) {
	double sw = 2*w - 1;
	ww *= 0.5*(1 + a1*sw + a3*sw*sw*sw);
      } else {
	ww *= w;
      }
      //      printf("x=%.2f, r=%d, a=%.2f, b=%.2f, c=%.2f, e=%.2f\n",x[r],r,val[index],b[idx],c[idx],e);
      double v;
      if(e == 2.0)
	v = val[index] + b[idx]*p + c[idx]*p*p;
      else if(e > 1.0)
	v = val[index] + b[idx]*p + c[idx]*pow(fabs(p),e);
      else
	v = val[index] + b[idx]*p + c[idx]*fabs(p);
      value += v;
      idx++;
    }
    Val += value*ww;
  }
  return Val*irank;
}

SEXP R_evalstalker(SEXP Sx, SEXP stalker, SEXP Smindeg, SEXP Smaxdeg, SEXP Ssmooth, SEXP Sthreads) {
  const double *x = REAL(Sx);
  const int K = nrows(Sx);
  const int N = ncols(Sx);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  const double *val = REAL(VECTOR_ELT(stalker,0));
  const double *b = REAL(VECTOR_ELT(stalker,1));
  const double *c = REAL(VECTOR_ELT(stalker,2));
  const double *d = REAL(VECTOR_ELT(stalker,3));
  SEXP Srank = VECTOR_ELT(stalker,4);
  const int *rank = INTEGER(Srank);
  const int nrank = LENGTH(Srank);
  double mindeg = REAL(AS_NUMERIC(Smindeg))[0];
  double maxdeg = REAL(AS_NUMERIC(Smaxdeg))[0];
  if(mindeg < 1) mindeg = 1.0;
  if(mindeg > 2) mindeg = 2.0;
  if(K != nrank) error("Rank of input(%d) does not match rank of spline",N,nrank);
  SEXP Sgrid = VECTOR_ELT(stalker,5);
  if(LENGTH(Sgrid) != nrank) error("Bad grid length %d , should be %d+",
				   LENGTH(Sgrid),nrank);
  double *grid[nrank];

  double smooth = REAL(Ssmooth)[0];
  smooth = 1-smooth;
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));
  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalstalker(x+i*nrank, nrank, rank, grid, val, b, c, d, mindeg, maxdeg, smooth);
  }
  UNPROTECT(1);
  return ret;
}

static void makestalker(int nrank, int *rank, double *val,
			double *b, double *c, double *d, int threads) {
  int N = 1;
  for(int r = 0; r < nrank; r++) N *= rank[r];
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 10 && threads > 1)
  for(int index = 0; index < N; index++) {
    int stride = 1;
    int idx = nrank*index;
    for(int r = 0; r < nrank; r++) {
      const int ridx = (index / stride) % rank[r];
      const int hi = index + stride;
      const int lo = index - stride;
      //    a[index] = val[index];
      d[idx] = 1.0;
      if(ridx == rank[r]-1) {
	// linear. Or should there be a double root at the end?
	b[idx] = val[index]-val[lo];
	c[idx] = 0.0;
      } else if(ridx == 0) {
	b[idx] = val[hi]-val[index];
	c[idx] = 0.0;
      } else {
	b[idx] = 0.5*(val[hi]-val[lo]);
	c[idx] = 0.5*(val[hi]+val[lo]) - val[index];
	const double ac = fabs(c[idx]), ab = fabs(b[idx]);
	if(ac != 0.0) {
	  if(ac <= ab && ab < 2.0*ac)
	    d[idx] = ab/ac;
	  else if(ab < ac && ac < 2.0*ab)
	    d[idx] = ac/ab;
	  else
	    d[idx] = 2.0;
	}
      }
      idx++;
      stride *= rank[r];
    }
  }
}

SEXP R_makestalker(SEXP Sval, SEXP Sgrid, SEXP Sthreads) {
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  int numprotect = 0;
  SEXP stalker = PROTECT(allocVector(VECSXP,6)); numprotect++;
  const int nrank = LENGTH(Sgrid);
  SEXP Srank = PROTECT(NEW_INTEGER(nrank)); numprotect++;
  int *rank = INTEGER(Srank);
  int N = 1;
  for(int i = 0; i < nrank; i++) {
    rank[i] = LENGTH(VECTOR_ELT(Sgrid,i));
    N *= rank[i];
  }
  if(N != LENGTH(Sval)) {
    error("Length of values (%d) does not match size of grid (%d)",LENGTH(Sval),N);
  }
  SEXP b = PROTECT(NEW_NUMERIC(N*nrank));numprotect++;
  SEXP c = PROTECT(NEW_NUMERIC(N*nrank));numprotect++;
  SEXP d = PROTECT(NEW_NUMERIC(N*nrank));numprotect++;
  SET_VECTOR_ELT(stalker, 0, Sval);
  SET_VECTOR_ELT(stalker, 1, b);
  SET_VECTOR_ELT(stalker, 2, c);
  SET_VECTOR_ELT(stalker, 3, d);
  SET_VECTOR_ELT(stalker, 4, Srank);
  SET_VECTOR_ELT(stalker, 5, Sgrid);
  makestalker(nrank, rank, REAL(Sval), REAL(b), REAL(c), REAL(d), threads);

  UNPROTECT(numprotect);
  return stalker;
}
