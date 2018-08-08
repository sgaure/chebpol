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
  double a1=1, a3 = 0;
  if(smooth > 0.0) {
    a3 = 0.5*(smooth-1.0);
    a1 = 1.0-a3;
  }

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
  ac = fabs(cplus); ab = fabs(bplus);
  if(ac != 0.0) {
    if(ac <= ab && ab < 2.0*ac)
      dplus = ab/ac;
    else if(ab < ac && ac < 2.0*ab)
      dplus = ac/ab;
    else
      dplus = 2.0;
  }
  // map degree [1,2] into [mindeg,maxdeg]
  dmin = mindeg + (maxdeg - mindeg)*(dmin - 1);
  dplus = mindeg + (maxdeg - mindeg)*(dplus - 1);

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
  if(smooth != 0.0) {
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
  SEXP Sgrid = VECTOR_ELT(stalker,1);
  const int nrank = LENGTH(Sgrid);
  if(K != nrank) error("Rank of input(%d) does not match rank of spline",N,nrank);

  int dims[nrank];
  int len = 1;
  for(int r = 0; r < nrank; r++) {
    dims[r] = LENGTH(VECTOR_ELT(Sgrid,r));
    len *= dims[r];
  }
  if(len != LENGTH(VECTOR_ELT(stalker,0))) {
    error("number of values (%d) does not match grid size (%d)\n",
	  LENGTH(VECTOR_ELT(stalker,0)), len);
  }
  double mindeg = REAL(AS_NUMERIC(Smindeg))[0];
  double maxdeg = REAL(AS_NUMERIC(Smaxdeg))[0];
  double *grid[nrank];
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));
  
  double smooth = REAL(AS_NUMERIC(Ssmooth))[0];
  smooth = 1-smooth;
  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalstalker(x+i*nrank, nrank, dims, grid, val, mindeg,maxdeg,smooth);
  }
  UNPROTECT(1);
  return ret;
}
