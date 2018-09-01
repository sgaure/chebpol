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

// Find the parameters for the hyperbolic stalker in each grid point
// use normalized coordinates as in the formulas, f(0) = 0
#define MYEPS (1e3*sqrt(DOUBLE_EPS))
static void makehyp(const int nrank, const int *dims, double **grid, const double *val,
		double *a, double *b, double *c, double *d, int depth, int *point) {
  // Loop over every grid point
  // In each grid point, find the function living there, i.e. an a,
  // and b,c,d for each dimension.
  // how many grid points are there? prod(dims)
  // We need nested loops to nrank depth. I.e. recursion.
  // record point
  if(depth < nrank) {
    for(int i = 0; i < dims[depth]; i++) {
      point[depth] = i;
      makehyp(nrank, dims, grid, val, a, b, c, d, depth+1,point);
    }
    return;
  }

    // What is the index of the point?  (We could do this in the recursion as well).
  R_xlen_t index = 0,stride = 1;  

  for(int i = 0; i < nrank; i++) {
    index += point[i]*stride;
    stride *=  dims[i];
  }

  // Inner loop, we have a grid point in 'point'
  // loop over the dimensions, collecting b, c, and d

  // Where to store them.
  double *bb = &b[nrank*index], *cc=&c[nrank*index], *dd=&d[nrank*index];
  a[index] = 0.0;
  stride = 1; // Redo this for each dimension
  for(int i = 0; i < nrank; i++) {
    double *gr = grid[i];
    int idx = point[i];

    // Check border, we do it linearly there
    if(idx == 0) {
      // Nothing below us, we are linear upwards
      bb[i] = (val[index+stride]-val[index])/(gr[idx+1]-gr[idx]);
      cc[i] = dd[i] = 0.0;

    } else if(idx == dims[i]-1) {
      // Nothing above us, linear downwards
      bb[i] = (val[index-stride] - val[index])/(gr[idx-1]-gr[idx]);
      cc[i] = dd[i] = 0.0;

    } else {
      // inside somewhere
      double v0 = val[index];
      double vplus = val[index+stride]-v0;
      double vmin = val[index-stride]-v0;
      double kplus = gr[idx+1] - gr[idx];
      double kmin = gr[idx-1] - gr[idx];
      // Many common expressions below. Fix later. Compiler does anyway.
      double D = kmin*vplus - kplus*vmin;
      double D2 = kmin*kmin*vplus - kplus*kplus*vmin;
      if(sign(vplus*vmin) <= 0 || fabs(vplus) <= MYEPS || fabs(vmin) <= MYEPS) {
	// Monotonic
	if(fabs(vplus-vmin) <= MYEPS*kplus) {
	  // constant zero
	  bb[i] = cc[i] = dd[i] = 0.0;
	} else if(fabs(D) <= MYEPS*kplus) {
	  // Linear
	  bb[i] = vplus/kplus;
	  cc[i] = dd[i] = 0.0;
	} else {
	  // general
	  bb[i] = 0.0;
	  cc[i] = -D/(kmin*kplus*(vplus-vmin));
	  dd[i] = -vplus*vmin*(kmin-kplus)/D;
	}
      } else {
	// Non-monotonic
	if(fabs(D2) <= MYEPS*kplus) {
	  // Parabola
	  bb[i] = vplus/(kplus*kplus);
	  dd[i] = 0.0;
	  cc[i] = NA_REAL; // signals parabola
	} else {
	  // general
	  bb[i] = vplus*vmin*(kmin-kplus)/D2;
	  cc[i] = -D2/(kmin*kplus*D);
	  dd[i] = -kmin*vplus*kplus*vmin*(kmin-kplus)*D/(D2*D2);
	}
      }
      /*
      if(isnan(dd[i]) || fabs(dd[i]) > 1e12) {
	printf("dd is na, i=%d, k-=%.2f, k+=%.2f, v-=%.2f, v+=%.2f D=%.2e D2=%.2e\n",
	       i,kmin,kplus,vmin,vplus,D,D2);
      }
      */
    }
    a[index] -= dd[i];
    stride *= dims[i];
  }
}

double evalhyp(const double *x, const int nrank, const int *dims, 
		  double **grid, const double *val, int blend,
		  double *a, double *b, double *c, double *d) {
  // loop over the corners of the hypercube surrounding x.
  // For each corner, evaluate the function there.
  // Make weighted sum.

  int stride = 1;
  double weight[nrank];
  int gbase[nrank];
  int llcorner = 0;
  // Find the lower left corner for each dimension.
  // Store in gpos
  for(int g = 0; g < nrank; g++) {
    int mflag;
    double *gr = grid[g];
    int gp = findInterval(gr, dims[g], x[g], TRUE,TRUE,1,&mflag)-1;
    // gp is below x
    gbase[g] = gp;
    weight[g] = 1-(x[g]-gr[gbase[g]])/(gr[gbase[g]+1]-gr[gbase[g]]);
    // This ensures extrapolation is linear
    if(weight[g] < 0) weight[g] = 0;
    if(weight[g] > 1) weight[g] = 1;
    llcorner += stride*gp;
    stride *= dims[g];
  }

  // Loop over the corners
  double V=0.0;
  double zx[nrank];
  double sumw = 0.0;
  for(int i = 0; i < (1<<nrank); i++) {
    // i is a corner represented as nrank bits
    // bit j is 1 if above in dimension j, 0 if below
    // Find the position and weight of the corner
    int stride = 1;
    int corner = llcorner;
    double cw = 1.0;
    for(int g = 0; g < nrank; g++) {
      double *gr = grid[g];
      if(i & (1<<g)) {
	cw *= blendfun(1-weight[g],blend);
	corner += stride;
	// zero based in this corner
	zx[g] = x[g] - gr[gbase[g]+1];
      } else {
	zx[g] = x[g] - gr[gbase[g]];
	cw *= blendfun(weight[g],blend);
      }
      stride *= dims[g];
    }
    if(cw == 0.0) continue;
    //    cw = blendfun(cw,blend);
    // Nice. Now 'corner' is the index into our arrays
    // Find a, b, c, and d for the corner. These have stride 'nrank'
    double *bb = &b[nrank*corner], *cc=&c[nrank*corner], *dd=&d[nrank*corner];
    // These are for a zero-based function, i.e. based on the corner being 0,
    // and f(0)=0
    // Compute the hyperbolic stalker in x
    // The zero-based x is in zx
    // 
    double fval = val[corner] + a[corner];
    for(int j = 0; j < nrank; j++) {
      // It might be a parabola, but typically a hyperbola
      if(!isnan(cc[j])) {
	fval += bb[j]*zx[j] + dd[j]/(1.0+cc[j]*zx[j]);
      } else {
	fval += bb[j]*zx[j]*zx[j];
      }
    }
    // That's the function value. The weight is in cw
    //    printf("corner %d fval %.2f weight %.2f\n",corner,fval,cw);
    sumw += cw;
    V += cw*fval;
  }
  return V/sumw;
}



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

SEXP R_makehyp(SEXP Sval, SEXP Sgrid) {
  const int nrank = LENGTH(Sgrid);
  int dims[nrank];
  int len = 1;
  for(int r = 0; r < nrank; r++) {
    dims[r] = LENGTH(VECTOR_ELT(Sgrid,r));
    len *= dims[r];
  }
  if(len != LENGTH(Sval)) {
    error("number of values (%d) does not match grid size (%d)\n",
	  LENGTH(Sval), len);
  }
  double *val = REAL(Sval);
  double *grid[nrank];
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));
  int point[nrank];
  SEXP Sa = PROTECT(NEW_NUMERIC(len));
  SEXP Sb = PROTECT(allocMatrix(REALSXP,nrank,len));
  SEXP Sc = PROTECT(allocMatrix(REALSXP,nrank,len));
  SEXP Sd = PROTECT(allocMatrix(REALSXP,nrank,len));
  SEXP ret = PROTECT(NEW_LIST(6));
  SEXP names = PROTECT(allocVector(STRSXP,6));
  SET_NAMES(ret,names);
  SET_VECTOR_ELT(ret,0,Sa); SET_STRING_ELT(names, 0, mkChar("a"));
  SET_VECTOR_ELT(ret,1,Sb); SET_STRING_ELT(names, 1, mkChar("b"));
  SET_VECTOR_ELT(ret,2,Sc);  SET_STRING_ELT(names, 2, mkChar("c"));
  SET_VECTOR_ELT(ret,3,Sd); SET_STRING_ELT(names, 3, mkChar("d"));
  SET_VECTOR_ELT(ret,4,Sval);  SET_STRING_ELT(names, 4, mkChar("val"));
  SET_VECTOR_ELT(ret,5,Sgrid); SET_STRING_ELT(names, 5, mkChar("grid"));
  
  makehyp(nrank, dims, grid, val, REAL(Sa), REAL(Sb), REAL(Sc), REAL(Sd), 0, point);
  UNPROTECT(6);
  return ret;
}

SEXP R_evalhyp(SEXP Sx, SEXP stalker, SEXP Sblend, SEXP Sthreads) {
  const double *x = REAL(Sx);
  //const int K = nrows(Sx);
  const int N = ncols(Sx);
  int blend = INTEGER(Sblend)[0];
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  double *a = REAL(VECTOR_ELT(stalker,0));
  double *b = REAL(VECTOR_ELT(stalker,1));
  double *c = REAL(VECTOR_ELT(stalker,2));
  double *d = REAL(VECTOR_ELT(stalker,3));
  double *val = REAL(VECTOR_ELT(stalker,4));
  SEXP Sgrid = VECTOR_ELT(stalker,5);
  const int nrank = LENGTH(Sgrid);
  double *grid[nrank];
  for(int i = 0; i < nrank; i++) grid[i] = REAL(VECTOR_ELT(Sgrid,i));
  int dims[nrank];
  int len = 1;
  for(int r = 0; r < nrank; r++) {
    dims[r] = LENGTH(VECTOR_ELT(Sgrid,r));
    len *= dims[r];
  }

  SEXP ret = PROTECT(NEW_NUMERIC(N));
  double *out = REAL(ret);

#pragma omp parallel for num_threads(threads) schedule(guided) if(N > 1 && threads > 1)
  for(int i = 0; i < N; i++)  {
    out[i] = evalhyp(x+i*nrank, nrank, dims, grid, val, blend, a,b,c,d);
  }
  
  UNPROTECT(1);
  return ret;
}


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
      error("Non-uniform grid in dimension %d not supported without GSL. Recompile with GSL.",i+1);
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
