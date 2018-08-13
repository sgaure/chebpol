#include "chebpol.h"

// Stalker spline with uniform grid, normalized coordinates
static R_INLINE double stalk1(double x, double vmin, double vplus,double dmin,
			      double dplus,double r, double iD,double powmin,double powplus) {
  double b,c;
  //  double iD = 1.0/(dmin*pow(dplus,r) + pow(dmin,r)*dplus);
  //  double b = (vplus*pow(dmin,r) - vmin*pow(dplus,r))*iD;
  //  double c = (vplus*dmin + vmin*dplus)*iD;
  if(!isnan(r)) {
    b = (vplus*powmin - vmin*powplus)*iD;
    c = (vplus*dmin + vmin*dplus)*iD;
    return b*x + c*pow(fabs(x),r);
  } else {

    // This works when dplus==dmin, otherwise we must solve an equation for r:
    // vmin*dplus^r - vplus*dmin^r = \pm r(vmin*dplus + vplus*dmin)
    // or (dplus^r - \pm r*dplus)/(dmin^r + \pm r*dmin) = vplus/vmin
    // Instead we map the interval [-dmin,0,dplus] -> [-1,0,1] with a rational function
    x = (dmin+dplus)*x/(2*dplus*dmin + (dplus-dmin)*x);
    //    x /= x>=0 ? dplus : dmin;
    b = 0.5*(vplus-vmin);
    c = 0.5*(vplus+vmin);
    
    if(c == 0.0) return b*x;
    
    // Find the degree. For uniform grids in normalized coordinates, this is quite easy
    double ac = fabs(c), ab = fabs(b);
    if(ac <= ab && ab < 2.0*ac) {
      r = ab/ac;
    } else if(ab < ac && ac < 2.0*ab) {
      r = ac/ab;
    } else {
      r = 2.0;
    }

    if(r == 2.0) {
      return b*x + c*x*x;
    } else if(r != 1.0) {
      return b*x + c*pow(fabs(x),r);
    } else {
      return b*x + c*fabs(x);
    }
  }
}

double evalstalker(const double *x, const int nrank, const int *dims, 
		   double **grid, const double *val,
		   const double *degree, double **det, double **pmin, double **pplus) {
  if(nrank == 0) return val[0];
  const int newrank = nrank-1;
  double *gr = grid[newrank];
  const int N = dims[newrank];
  const double xx = x[newrank];
  double *rdet = det[newrank], *rpmin = pmin[newrank], *rpplus = pplus[newrank];
    
  int mflag;
  int stride = 1;
  for(int r = 0; r < newrank; r++) stride *= dims[r];
  // Find the interval of the last coordinate, convert to zero-based
  const int imin = findInterval(gr, N, xx,TRUE,TRUE,1,&mflag)-1;

  // Now, find the function values on the two grid points on each side.
  // i.e. imin, imin-1, and imin+1 and imin+2
  // use these to create the stalker splines

  // values of the lower dimensional spline
  double v1=NA_REAL,v2,v3,v4=NA_REAL;
  if(newrank > 0) {
    if(imin > 0)
      v1 = evalstalker(x,newrank,dims,grid,val+(imin-1)*stride,degree,det,pmin,pplus);
    v2 = evalstalker(x,newrank,dims,grid,val+imin*stride,degree,det,pmin,pplus);
    v3 = evalstalker(x,newrank,dims,grid,val+(imin+1)*stride,degree,det,pmin,pplus);
    if(imin < N-2)
      v4 = evalstalker(x,newrank,dims,grid,val+(imin+2)*stride,degree,det,pmin,pplus);
  } else {    
    // save a recursion step
    if(imin > 0) v1 = val[(imin-1)*stride];
    v2 = val[imin*stride];
    v3 = val[(imin+1)*stride];
    if(imin < N-2) v4 = val[(imin+2)*stride];
  }
  // map into normalized coordinates
  //  const double nx = (xx-gr[imin])/(gr[imin+1]-gr[imin]);
  const double nx = (xx-gr[imin])/(gr[imin+1]-gr[imin]);
  double low=NA_REAL;
  double dmin = gr[imin]-gr[imin-1];
  double dplus = gr[imin+1]-gr[imin];
  if(imin > 0) low = v2+stalk1(xx-gr[imin],v1-v2,v3-v2,dmin,dplus,degree[newrank],
			       rdet[imin],rpmin[imin],rpplus[imin]);
  double high=NA_REAL;
  dmin = dplus;
  dplus = gr[imin+2]-gr[imin+1];
  if(imin < N-2) high = v3+stalk1(xx-gr[imin+1],v2-v3,v4-v3,dmin,dplus,degree[newrank],
				  rdet[imin+1],rpmin[imin+1],rpplus[imin+1]);
  /*
  const double low = stalk1(nx,v1,v2,v3,mindeg[newrank],maxdeg[newrank]);
  const double high = stalk1(nx-1,v2,v3,v4,mindeg[newrank],maxdeg[newrank]);
  */
  // combine the basis functions
  double w = 1-nx;
  if(imin == 0) w = 0;  // no basis in the end points
  if(imin == N-2) w = 1;
  // This ensures linear extrapolation:
  if(w < 0) w = 0.0;
  if(w > 1) w = 1.0;

  if(w == 0) return high;
  if(w == 1) return low;
  w = w<0.5 ? 0.5*exp(2-1/w) : (1-0.5*exp(2-1/(1-w)));  // smooth sigmoid blending
  return w*low + (1-w)*high;
}

SEXP R_evalstalker(SEXP Sx, SEXP stalker, SEXP Sdegree, SEXP Sthreads) {
  const double *x = REAL(Sx);
  const int K = nrows(Sx);
  const int N = ncols(Sx);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  const double *val = REAL(VECTOR_ELT(stalker,0));
  SEXP Sgrid = VECTOR_ELT(stalker,1);
  SEXP Sdet = VECTOR_ELT(stalker,2);
  SEXP Spmin = VECTOR_ELT(stalker,3);
  SEXP Spplus = VECTOR_ELT(stalker,4);
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
  if(LENGTH(Sdegree) != nrank) {
    error("degree length must match dimension (%d)",nrank);
  }
  double *degree = REAL(AS_NUMERIC(Sdegree));
  double *grid[nrank];
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));

  double *det[nrank];
  for(int i = 0; i < nrank; i++) det[i] = REAL(VECTOR_ELT(Sdet,i));
  double *pmin[nrank];
  for(int i = 0; i < nrank; i++) pmin[i] = REAL(VECTOR_ELT(Spmin,i));
  double *pplus[nrank];
  for(int i = 0; i < nrank; i++) pplus[i] = REAL(VECTOR_ELT(Spplus,i));

  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalstalker(x+i*nrank, nrank, dims, grid, val, degree, det, pmin, pplus);
  }
  UNPROTECT(1);
  return ret;
}
