#include "chebpol.h"
#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#endif

typedef struct {double dmin, dplus, vmin, vplus; int pos;} rblock;
typedef struct {double *det, *pmin, *pplus;int uniform;} precomp;
#ifdef HAVE_GSL
static double rfun(double r, void *P) {
  rblock *p = (rblock*) P;
  const double dmin=p->dmin, dplus=p->dplus, vmin=p->vmin, vplus=p->vplus;
  int pos=p->pos;
  if(pos) 
    return vmin*pow(dplus,r) - vplus*pow(dmin,r) - r*pow(dplus,r-1)*(vmin*dplus+vplus*dmin);
  else
    return vmin*pow(dplus,r) - vplus*pow(dmin,r) + r*pow(dmin,r-1)*(vmin*dplus+vplus*dmin);
}
#endif

#if 0
//Could be faster, but it runs astray
static double rfun_deriv(double r, void *P) {
  rblock *p = (rblock*) P;
  const double dmin=p->dmin, dplus=p->dplus, vmin=p->vmin, vplus=p->vplus;
  int pos=p->pos;
  if(pos) 
    return vmin*r*pow(dplus,r-1) - vplus*r*pow(dmin,r-1) - r*(r-1)*pow(dplus,r-2)*(vmin*dplus+vplus*dmin);
  else
    return vmin*r*pow(dplus,r-1) - vplus*r*pow(dmin,r-1) + r*(r-1)*pow(dmin,r-2)*(vmin*dplus+vplus*dmin);
}

static void rfun_fdf(double r, void *P, double *y, double *dy) {
  rblock *p = (rblock*) P;
  const double dmin=p->dmin, dplus=p->dplus, vmin=p->vmin, vplus=p->vplus;
  int pos=p->pos;
  if(pos) {
    *y = vmin*pow(dplus,r) - vplus*pow(dmin,r) - r*pow(dplus,r-1)*(vmin*dplus+vplus*dmin);
    *dy = vmin*r*pow(dplus,r-1) - vplus*r*pow(dmin,r-1) - r*(r-1)*pow(dplus,r-2)*(vmin*dplus+vplus*dmin);
  } else {
    *y = vmin*pow(dplus,r) - vplus*pow(dmin,r) + r*pow(dmin,r-1)*(vmin*dplus+vplus*dmin);
    *dy = vmin*r*pow(dplus,r-1) - vplus*r*pow(dmin,r-1) + r*(r-1)*pow(dmin,r-2)*(vmin*dplus+vplus*dmin);
  }
}
#endif

#ifdef HAVE_GSL
static R_INLINE double findmonor(double dmin,double dplus, double vmin, double vplus) {
// We should really cache results here, it could save time when computing r for the same interval
// We can use std::unordered_map for that. We also need some parallel sync, std::mutex. Later.  
  double D2 = (vmin*dplus*dplus - vplus*dmin*dmin)/(vplus*dmin+dplus*vmin);
  double r;
  if(D2 > 2*dplus || D2 < -2*dmin) {
    r = 2.0;
  } else {
    // Find the largest r which keeps it monotonic, i.e. derivative 0 in one
    // of the end points
    rblock rb = {dmin,dplus,vmin,vplus,0};
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    int status;
    F.function = &rfun;
    F.params = &rb;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, 1, 2);
    rb.pos = 1;
    double t1 = rfun(1,&rb), t2 = rfun(2,&rb);
    if(sign(t1*t2) > 0) rb.pos = 0;
    int iter = 0;
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      double lo = gsl_root_fsolver_x_lower(s);
      double hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(lo,hi,0,1e-8);
    } while(status == GSL_CONTINUE && iter < 50);
    if(status != GSL_SUCCESS) {
      // I need to think about what to do here. We're in openmp, so can't do R-calls.
      //      printf("No convergence\n");
    }
    gsl_root_fsolver_free(s);
  }
  return r;
}
#endif

// hyperbolic stalker
static R_INLINE double stalkhyp(double x, double vmin, double vplus,double dmin, double dplus) {
  //  double iD = 1.0/(dmin*pow(dplus,r) + pow(dmin,r)*dplus);
  //  double b = (vplus*pow(dmin,r) - vmin*pow(dplus,r))*iD;
  //  double c = (vplus*dmin + vmin*dplus)*iD;
  if(fabs(x - dplus) < 1e3*DOUBLE_EPS*dplus) return vplus;
  if(fabs(x + dmin) < 1e3*DOUBLE_EPS*dmin) return vmin;
  // hyperbolic function
  // a+bx-a/(cx+1)
  if(fabs(vplus*dmin + vmin*dplus) < 1e3*DOUBLE_EPS*dmin*dplus) return x*(vplus-vmin)/(dplus+dmin); //linear
  if(sign(vplus*vmin) < 0) {
    //monotonic, set b=0
    double a,c;
    if(fabs(vmin-vplus) < 1e3*DOUBLE_EPS) return 0;
#if 1
    c = (vmin*dplus + vplus*dmin)/(dmin*dplus*(vmin-vplus));
    a = vplus + vplus/(c*dplus);
    return a - a/(1+c*x);
#else
    // Alternatively, make one of the end points have zero derivative
    double b;
    double D = vplus*dmin + vmin*dplus;
    double D2hi = dplus*dplus*vmin + dmin*dmin*vplus + 2*dplus*dmin*vplus;
    double D2lo = dplus*dplus*vmin + dmin*dmin*vplus + 2*dplus*dmin*vmin;
    if(sign(D*vplus)<0) {
      // choose zero derivative in dplus:
      a = vplus*(dplus+dmin)*(dplus+dmin) * D*D/(D2hi*D2hi);
      b = vplus*dmin*(vplus-vmin)/D2hi;
      c = -D2hi/(dmin*dplus*dplus*(vplus-vmin));
    } else {
      a = vmin*(dplus+dmin)*(dplus+dmin) * D*D/(D2lo*D2lo);
      b = vmin*dplus*(vplus-vmin)/D2lo;
      c = -D2lo/(dmin*dmin*dplus*(vplus-vmin));
    }
    return a + b*x - a/(1+c*x);
#endif    
  } else {
    // non-monotonic
    double D = dplus*dplus*vmin - dmin*dmin*vplus;
    
    if(fabs(D) < 1e3*DOUBLE_EPS) {
      double a;
      // if c==0, it's a parabola y = a*x^2
      // we have (0,0), (-dmin,vmin), and (dplus,vplus)
      // vmin = a*dmin^2
      // vplus = a*dplus^2
      a = vmin/(dmin*dmin);
      return a*x*x;
    }
    double c = D/(dmin*dplus*(dmin*vplus+dplus*vmin));
    double a,b;
    b = vplus*(1+c*dplus)/(c*dplus*dplus);
    a = -b/c;
    return a + b*x - a/(1+c*x);
  }
}

// Variable degree stalker spline
static R_INLINE double stalk1(double x, double vmin, double v0, double vplus,double dmin,
			      double dplus,double r, double iD,double powmin,double powplus, int uniform) {
  double b,c;
  //  double iD = 1.0/(dmin*pow(dplus,r) + pow(dmin,r)*dplus);
  //  double b = (vplus*pow(dmin,r) - vmin*pow(dplus,r))*iD;
  //  double c = (vplus*dmin + vmin*dplus)*iD;
  // Do linear end points

  if(isnan(vmin)) return v0 + x/dplus * (vplus-v0);
  if(isnan(vplus)) return v0 + x/dmin * (v0-vmin);
  vmin -= v0;
  vplus -= v0;
  if(r == 0) return v0+stalkhyp(x,vmin,vplus,dmin,dplus);
  if(!isnan(r)) {
    // use precomputed degree
    b = (vplus*powmin - vmin*powplus)*iD;
    c = (vplus*dmin + vmin*dplus)*iD;
    return v0 + b*x + c*pow(fabs(x),r);
  } else {
    // Compute r for this basis function
#ifdef HAVE_GSL
    if(uniform) {
#endif
      // Uniform grid
      // This works when dplus==dmin, otherwise we must solve an equation for r.
      x /= dplus;
      b = 0.5*(vplus-vmin);
      c = 0.5*(vplus+vmin);
    
      if(c == 0.0) return v0 + b*x;
      
      // Find the degree. For uniform grids in normalized coordinates, this is quite easy
      double ac = fabs(c), ab = fabs(b);
      if(sign(vplus*vmin) <= 0) {
	// monotonic
	r = (ab < 2.0*ac) ? ab/ac : 2.0;
      } else {
	// non-monotonic
	// pretend sign change of vplus for computation of r, that is, interchange |b| and |c|
	r = (ac < 2.0*ab) ? ac/ab : 2.0;
      }
      
      if(r == 2.0) {
	return v0 + b*x + c*x*x;
      } else if(r != 1.0) {
	return v0 + b*x + c*pow(fabs(x),r);
      } else {
	return v0 + b*x + c*fabs(x);
      }
#ifdef HAVE_GSL
    } else {

      // We must solve for r:
      double b,c;
      if(sign(vmin*vplus) < 0) {
	// monotonic
	r = findmonor(dmin,dplus,vmin,vplus);
      } else {
	// non-monotonic, should we keep f'(0) = 0? No, won't work.
	// We can use r=2, but if we just precisely is non-monotonic
	// we should use a smaller r to avoid overshoot. But which?
	// Keep it in some way "symmetric" with the monotonic case?  How?
	// Change sign on one of the values, compute the r which
	// fits these monotonic points.
	r = findmonor(dmin,dplus,-vmin,vplus);
      }
      c = (dplus*vmin + dmin*vplus)/(pow(dmin,r)*dplus + pow(dplus,r)*dmin);
      b = (vplus - c*pow(dplus,r))/dplus;
      return v0 + b*x + c*pow(fabs(x),r);
    }
#endif
  }
}
static R_INLINE double blendfun(double w,int blend) {
  switch(blend) {
  case 1:
    w = w<0.5 ? 0.5*exp(2-1/w) : 1-0.5*exp(2-1/(1-w));  // smooth sigmoid blending
    break;
  case 2:
    w = w<0.5 ? 0.5*exp(4-1/(w*w)) : 1-0.5*exp(4-1/((1-w)*(1-w)));  // parodic sigmoid blending
    break;
  case 3:
    w = (-2*w + 3)*w*w; // cubic blending
    break;
  case 4:
    w = (w < 0.5) ? 0 : 1;  // discontinuous blending
  }    
  return w;
}

#if 1
double evalstalker(const double *x, const int nrank, const int *dims, 
		   double **grid, const double *val, int blend,
		   const double *degree, precomp *pcomp) {
		   
  if(nrank == 0) return val[0];
  const int newrank = nrank-1;
  double *gr = grid[newrank];
  const int N = dims[newrank];
  const double xx = x[newrank];
  //  double *rdet = det[newrank], *rpmin = pmin[newrank], *rpplus = pplus[newrank];
  double *rdet = pcomp[newrank].det, *rpmin = pcomp[newrank].pmin, *rpplus = pcomp[newrank].pplus;
  int uniform = pcomp[newrank].uniform;
  int mflag;
  int stride = 1;
  for(int r = 0; r < newrank; r++) stride *= dims[r];
  // Find the interval of the last coordinate, convert to zero-based
  const int imin = findInterval(gr, N, xx,TRUE,TRUE,1,&mflag)-1;

  // Now, find the function values on the two grid points on each side.
  // i.e. imin, imin-1, and imin+1 and imin+2
  // use these to create the stalker splines across this dimension

  // values of the lower dimensional spline
  double v1=NA_REAL,v2,v3,v4=NA_REAL;
  if(newrank > 0) {
    if(imin > 0)
      v1 = evalstalker(x,newrank,dims,grid,val+(imin-1)*stride,blend,degree,pcomp);
    v2 = evalstalker(x,newrank,dims,grid,val+imin*stride,blend,degree,pcomp);
    v3 = evalstalker(x,newrank,dims,grid,val+(imin+1)*stride,blend,degree,pcomp);
    if(imin < N-2)
      v4 = evalstalker(x,newrank,dims,grid,val+(imin+2)*stride,blend,degree,pcomp);
  } else {    
    // save a recursion step
    if(imin > 0) v1 = val[(imin-1)*stride];
    v2 = val[imin*stride];
    v3 = val[(imin+1)*stride];
    if(imin < N-2) v4 = val[(imin+2)*stride];
  }
  double low=NA_REAL;
  double dmin=NA_REAL;
  if(imin > 0) dmin = gr[imin]-gr[imin-1];
  double dplus = gr[imin+1]-gr[imin];
  low = stalk1(xx-gr[imin],v1,v2,v3,dmin,dplus,degree[newrank],
    	       rdet[imin],rpmin[imin],rpplus[imin],uniform);
  double high=NA_REAL;
  dmin = dplus;
  if(imin < N-2) dplus = gr[imin+2]-gr[imin+1]; else dplus = NA_REAL;
  high = stalk1(xx-gr[imin+1],v2,v3,v4,dmin,dplus,degree[newrank],
		rdet[imin+1],rpmin[imin+1],rpplus[imin+1],uniform);

  // combine the basis functions
  // map into normalized coordinates
  const double nx = (xx-gr[imin])/(gr[imin+1]-gr[imin]);
  double w = 1-nx;
  //  if(imin == 0) w = 0;  // no basis in the end points
  //  if(imin == N-2) w = 1;

  if(w == 0) return high;
  if(w == 1) return low;
  w = blendfun(w,blend);
  return w*low + (1-w)*high;
}
#else

// more like the multilinear, use only stalkers on the grid lines,
// do linear interpolation between.
double evalstalker(const double *x, const int nrank, const int *dims, 
		   double **grid, const double *val, int blend,
		   const double *degree, precomp *pcomp) { 
  // loop over the corners
  // for each corner, stalk the value in each dimension.
  // take weighted averages
  // First, find the weights
  int stride = 1;
  double weight[nrank];
  double xx[nrank];
  int gpos[nrank];
  int llcorner = 0;
  for(int g = 0; g < nrank; g++) {
    int mflag;
    double *gr = grid[g];
    int gp = findInterval(gr, dims[g], x[g], TRUE,TRUE,1,&mflag)-1;
    // gp is at or below x
    gpos[g] = gp;
    xx[g] = x[g]-gr[gp];
    weight[g] = xx[g]/(gr[gp+1]-gr[gp]);
    llcorner += stride*gp;
    stride *= dims[g];
  }
  // Fine. That's the weights.
  // Now, loop over corners
  double V = 0.0;
  double sumw = 0.0;
  for(int i = 0; i < (1<<nrank); i++) {
    // i is a corner represented as nrank bits
    // bit j is 1 if above in dimension j, 0 if below
    // Find the position and weight of the corner
    int stride = 1;
    int corner = llcorner;
    double cw = 1.0;
    for(int g = 0; g < nrank; g++) {
      if(i & (1<<g)) {
	cw *= weight[g];
	corner += stride;
      } else {
	cw *= (1.0-weight[g]);
      }
      stride *= dims[g];
    }
    // Now, loop over the dimensions, stalk it.
    stride = 1;
    double cV = 0;
    double ccw = cw;
    for(int g = 0; g < nrank; g++) {
      double *gr = grid[g];
      double *rdet = pcomp[g].det, *rpmin = pcomp[g].pmin, *rpplus = pcomp[g].pplus;
      int uniform = pcomp[g].uniform;
      //      double *rdet = det[g], *rpmin=pmin[g], *rpplus=pplus[g];
      int pos = gpos[g];
      const double *vptr = &val[corner];
      double vmin=NA_REAL,v0=vptr[0],vplus=NA_REAL,dmin=NA_REAL,dplus=NA_REAL;
      double w;
      int above = i & (1<<g);
      if(above) {
	// corner is above point, stalk down
	pos++;
	vmin = vptr[-stride];
	dmin = gr[pos] - gr[pos-1];
	if(pos < dims[g]-1) {
	  vplus = vptr[stride];
	  dplus = gr[pos+1]-gr[pos];
	}
	w = weight[g];
      } else {
	// corner is below point, stalk up
	vplus = vptr[stride];
	dplus = gr[pos+1]-gr[pos];
	if(pos > 0) {
	  vmin = vptr[-stride];
	  dmin = gr[pos]-gr[pos-1];
	}
	w = 1-weight[g];
      }
      double bw = blendfun(w,blend);
      ccw = ((w>0)?bw/w:1);
      cV += ccw*stalk1(x[g]-gr[pos],vmin,v0,vplus,dmin,dplus,degree[g],rdet[pos],rpmin[pos],rpplus[pos],uniform);
      stride *= dims[g];
    }
    V += cw*cV;
    sumw += cw;
  }
  return V/nrank;
}
#endif
SEXP R_evalstalker(SEXP Sx, SEXP stalker, SEXP Sdegree, SEXP Sblend, SEXP Sthreads) {
  const double *x = REAL(Sx);
  const int K = nrows(Sx);
  const int N = ncols(Sx);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  const double *val = REAL(VECTOR_ELT(stalker,0));
  SEXP Sgrid = VECTOR_ELT(stalker,1);
  SEXP Sdet = VECTOR_ELT(stalker,2);
  SEXP Spmin = VECTOR_ELT(stalker,3);
  SEXP Spplus = VECTOR_ELT(stalker,4);
  SEXP Suniform = VECTOR_ELT(stalker,5);
  int blend = INTEGER(AS_INTEGER(Sblend))[0];
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
  precomp pcomp[nrank];
  for(int i = 0; i < nrank; i++) {
    pcomp[i].det = REAL(VECTOR_ELT(Sdet,i));
    pcomp[i].pmin = REAL(VECTOR_ELT(Spmin,i));
    pcomp[i].pplus = REAL(VECTOR_ELT(Spplus,i));
    pcomp[i].uniform = LOGICAL(Suniform)[i];
#ifndef HAVE_GSL
    if(!pcomp[i].uniform && isnan(degree[i]))
      error('Non-uniform grid in dimension %d not supported without GSL. Recompile with GSL.',i+1);
#endif
  }
#ifdef HAVE_GSL
  gsl_error_handler_t *gslh = gsl_set_error_handler_off();
#endif
  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);
#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalstalker(x+i*nrank, nrank, dims, grid, val, blend, degree, pcomp);
  }
#ifdef HAVE_GSL
  gsl_set_error_handler(gslh);
#endif
  UNPROTECT(1);
  return ret;
}

SEXP havegsl() {
#ifdef HAVE_GSL
  return ScalarLogical(TRUE);
#else
  return ScalarLogical(FALSE);
#endif
}
