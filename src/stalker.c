#include "chebpol.h"

static const double EC = M_E/(M_E*M_E + 1 - 2.0*M_E);
// Stalker spline with uniform grid, normalized coordinates
static R_INLINE double stalk1(double x, double Vmin, double V0, double Vplus,double dmin,double dplus) {
  // Is it endpoint, i.e. linear? Are these NaNs?
  double r = 1;
  r = dmin/dplus;
  double vmin = Vmin-V0;
  double vplus = Vplus-V0;
  if(isnan(vmin) || r==0.0) return V0 + (Vplus-V0)*x/dplus;
  if(isnan(vplus)) return V0 + (V0-Vmin)*x/dmin;


  double a = (vplus*dmin - vmin*dplus)/(dmin*(1-exp(r*dplus)) - dplus*(1-exp(r*dmin)));
  double b = (vplus-a*(1-exp(r*dplus)))/dplus;
  double c = -a;
  return V0 + a + b*x + c*exp(r*x);
  /*
  double c = (vplus+vmin - 2*v0)/(exp(r)+exp(-r)-2);
  double a = v0 - c;
  double b = vplus - a - c*exp(r);
  return a + b*x + c*exp(r*x);
  */
  /*
  double b = 0.5*(vplus-vmin);
  double c = 0.5*(vplus+vmin) - v0;
  
  if(c == 0.0) return v0 + b*x;
  
  // Find the degree. For uniform grids in normalized coordinates, this is quite easy
  double ac = fabs(c), ab = fabs(b), d;
  if(ac <= ab && ab < 2.0*ac) {
    d = ab/ac;
  } else if(ab < ac && ac < 2.0*ab) {
    d = ac/ab;
  } else {
    d = 2.0;
  }
  // map degree [1,2] into [mindeg,maxdeg]
  //  d = mindeg + (maxdeg - mindeg)*(d - 1.0);
  if(d < mindeg) d = mindeg;
  if(d > maxdeg) d = maxdeg;
  if(d == 2.0) {
    return v0 + b*x + c*x*x;
  } else if(d != 1.0) {
    return v0 + b*x + c*pow(fabs(x),d);
  } else {
    return v0 + b*x + c*fabs(x);
  }
  */
}

double evalstalker(const double *x, const int nrank, const int *dims, 
		 double **grid, const double *val,
		 const double *mindeg, const double *maxdeg) {
  if(nrank == 0) return val[0];
  const int newrank = nrank-1;
  double *gr = grid[newrank];
  const int N = dims[newrank];
  const double xx = x[newrank];
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
      v1 = evalstalker(x,newrank,dims,grid,val+(imin-1)*stride,mindeg,maxdeg);
    v2 = evalstalker(x,newrank,dims,grid,val+imin*stride,mindeg,maxdeg);
    v3 = evalstalker(x,newrank,dims,grid,val+(imin+1)*stride,mindeg,maxdeg);
    if(imin < N-2)
      v4 = evalstalker(x,newrank,dims,grid,val+(imin+2)*stride,mindeg,maxdeg);
  } else {    
    // save a recursion step
    if(imin > 0) v1 = val[(imin-1)*stride];
    v2 = val[imin*stride];
    v3 = val[(imin+1)*stride];
    if(imin < N-2) v4 = val[(imin+2)*stride];
  }
  // map into normalized coordinates
  //  const double nx = (xx-gr[imin])/(gr[imin+1]-gr[imin]);
  const double nx = (xx-gr[imin]);
  const double low = stalk1(xx-gr[imin],v1,v2,v3,gr[imin-1]-gr[imin],gr[imin+1]-gr[imin]);
  const double high = stalk1(xx-gr[imin+1],v2,v3,v4,gr[imin]-gr[imin+1],gr[imin+2]-gr[imin+1]);
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
  w = exp(2-w-1/w) - exp(1+w-1/(1-w));
  return w*low + (1-w)*high;
}

SEXP R_evalstalker(SEXP Sx, SEXP stalker, SEXP Smindeg, SEXP Smaxdeg, SEXP Sthreads) {
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
  if(LENGTH(Smindeg) != nrank || LENGTH(Smaxdeg) != nrank) {
    error("degree length must match dimension (%d)",nrank);
  }
  double *mindeg = REAL(AS_NUMERIC(Smindeg));
  double *maxdeg = REAL(AS_NUMERIC(Smaxdeg));
  double *grid[nrank];
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));

  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalstalker(x+i*nrank, nrank, dims, grid, val, mindeg,maxdeg);
  }
  UNPROTECT(1);
  return ret;
}
