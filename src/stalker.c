#include "chebpol.h"

double evalstalker(const double *x, const int nrank, const int *dims, 
		 double **grid, const double *val,
		 const double mindeg, const double maxdeg, const double smooth) {
  if(nrank == 0) return val[0];
  double *gr = grid[nrank-1];
  const int N = dims[nrank-1];
  const double xx = x[nrank-1];
  int mflag;
  int stride = 1;
  for(int r = 0; r < nrank-1; r++) stride *= dims[r];
  // Find the interval of the last coordinate
  const int imin = findInterval(gr, N, xx,TRUE,TRUE,1,&mflag)-1;
  // normalized coordinates
  const double nx = (xx-gr[imin])/(gr[imin+1]-gr[imin]);
  // smoothing coordinates
  const double a3 = 0.5*(smooth - 1);
  const double a1 = 1-a3;

  // Now, find the function values on the two grid points on each side.
  // i.e. imin, imin-1, and imin+1 and imin+2
  // use these to create the stalker splines

  // values of the lower dimensional spline
  double vmin,vmin2,vplus,vplus2;
  vmin = evalstalker(x,nrank-1,dims,grid,val+imin*stride,mindeg,maxdeg,smooth);
  if(imin > 0)
    vmin2 = evalstalker(x,nrank-1,dims,grid,val+(imin-1)*stride,mindeg,maxdeg,smooth);
  vplus = evalstalker(x,nrank-1,dims,grid,val+(imin+1)*stride,mindeg,maxdeg,smooth);
  if(imin < N-1)
    vplus2 = evalstalker(x,nrank-1,dims,grid,val+(imin+2)*stride,mindeg,maxdeg,smooth);

  // Find the coefficients of the bases
  double bmin,cmin,dmin,bplus,cplus,dplus;
  if(imin == 0) {
    // linear left basis
    bmin = vplus - vmin;
    cmin = 0.0;
    // ordinary right basis
    bplus = 0.5*(vplus2-vmin);
    cplus = 0.5*(vplus2+vmin) - vplus;
  } else if(imin == N-2) {
    // ordinary left basis
    bmin = 0.5*(vplus-vmin2);
    cmin = 0.5*(vplus+vmin2) - vmin;
    // linear right basis
    bplus = vplus - vmin;
    cplus = 0.0;
  } else {
    // both ordinary
    bmin = 0.5*(vplus-vmin2);
    cmin = 0.5*(vplus+vmin2) - vmin;
    bplus = 0.5*(vplus2-vmin);
    cplus = 0.5*(vplus2+vmin) - vplus;
  }


  // Find the degrees of the bases
  double ac = fabs(cmin), ab = fabs(bmin);
  if(ac != 0.0) {
    if(ac <= ab && ab < 2.0*ac)
      dmin = ab/ac;
    else if(ab < ac && ac < 2.0*ab)
      dmin = ac/ab;
    else
      dmin = 2.0;
  }
  if(dmin < mindeg) dmin = mindeg;
  if(dmin > maxdeg) dmin = maxdeg;
  ac = fabs(cplus); ab = fabs(bplus);
  if(ac != 0.0) {
    if(ac <= ab && ab < 2.0*ac)
      dplus = ab/ac;
    else if(ab < ac && ac < 2.0*ab)
      dplus = ac/ab;
    else
      dplus = 2.0;
  }
  if(dplus < mindeg) dplus = mindeg;
  if(dplus > maxdeg) dplus = maxdeg;

  // evaluate the basis functions
  double low,high;
  if(dmin == 2.0) {
    low = vmin + bmin*nx + cmin*nx*nx;
  } else if(dmin > 1.0) {
    low = vmin + bmin*nx + cmin*pow(fabs(nx),dmin);
  } else {
    low = vmin + bmin*nx + cmin*fabs(nx);
  }
  if(dplus == 2.0) {
    high = vplus + bplus*(nx-1) + cplus*(nx-1)*(nx-1);
  } else if(dplus > 1.0) {
    high = vplus + bplus*(nx-1) + cplus*pow(fabs(nx-1),dplus);
  } else {
    high = vplus + bplus*(nx-1) + cplus*fabs(nx-1);
  }

  // combine the basis functions
  double w = 1-nx;
  if(smooth != 1.0) {
    double sw = 2*w - 1;
    w = 0.5*(1 + a1*sw + a3*sw*sw*sw);
  } 
  return w*low + (1-w)*high;
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
    out[i] = evalstalker(x+i*nrank, nrank, rank, grid, val, mindeg,maxdeg,smooth);
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
