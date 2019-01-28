//---------------------------------------------------------------------------
//
// $Id: SWSHSpheroid.cc,v 1.28 2019/01/26 17:57:06 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// The functions in this file are used to compute the spin weighted
// spheroidal harmonics.  See my notes for mathematical details.
// These functions are used by the class SWSH.
//
// Scott Hughes, 24 July 1998
//

#include <cmath>
#include <fftw3.h>
#include "Globals.h"
#include "Tensors.h"
#include "SWSH.h"

SWSH::SWSH(const int ss, const int ll, const int mm, const Real spintimesfreq) :
  Gl(ll), s(ss), Gm(mm), aw(spintimesfreq)
{
  lmin = Max(abs(s), abs(Gm));

  expand(&E, b, &N);
  lambda = E - 2.*Gm*aw + aw*aw - ((Real)(s*(s + 1)));

  cheb_spheroid();
  //
  // For the s == -2 code, we also need some operators on the spheroid
  // in order to build the Teukolsky source.
  //
  if (s == -2) {
    cheb_l2dagspheroid();
    cheb_l1dagl2dagspheroid();
  }
}

//
// Wrapper function for load_M and findeig.  This function starts out
// keeping a small number of terms in b, and then allows b to grow until
// findeig says that the b vector is OK.
//
void SWSH::expand(Real *E, Real *b, int *n)
{
  Real *M, norm = 0.0, maxcof;
  int okflag = 0, i, K, maxind = 1;

  if (Gl < abs(s) || Gl < abs(Gm)) {
    cerr << "Illegal values of m or s in expand" << endl;
    exit(0);
  }
  K = Gl - lmin + 5;
  while(!okflag) {
    K += 2;
    M = Tensor<Real>::vector(0, K*K - 1);
    load_M(M, K);
    okflag = findeig(M, K, E, b);
    Tensor<Real>::free_vector(M, 0, K*K - 1);
  }
  maxcof = 0.0;
  for (i = 1; i <= K; i++) {
    norm += b[i]*b[i];
    if (fabs(b[i]) > maxcof) {
      maxcof = fabs(b[i]);
      maxind = i;
    }
  }
  norm = sqrt(norm);
  if (b[maxind] < 0.0) norm *= -1.0;
  // Only normalize if the norm is nonzero...
  if (fabs(norm) > 1e-15)
    for (i = 1; i <= K; i++)
      b[i] /= norm;
  *n = K;
}

void SWSH::load_M(Real *M, const int K)
{
  Clebsch cg;
  int i, j, p, q;

  for (i = 0; i < K; i++) {
    p = i + lmin;
    for (j = 0; j < K; j++) {
      q = j + lmin;
      if (i == j) {
	M[j + i*K] = aw*aw*cg.xsqrbrac(s, p, q, Gm) - 
	  2.0*aw*((Real)s)*cg.xbrac(s, p, q, Gm)
	   - ((Real)(p*(p + 1)));
      }
      else if (i == j - 1 || i == j + 1)
	M[j + i*K] = aw*aw*cg.xsqrbrac(s, p, q, Gm) - 
	  2.0*aw*((Real)s)*cg.xbrac(s, p, q, Gm);
      else if (i == j - 2 || i == j + 2)
	M[j + i*K] = aw*aw*cg.xsqrbrac(s, p, q, Gm);
      else M[j + i*K] = 0.0;
    }
  }
}

//
// This function finds the eigenvalues and eigenvectors of M.
//
// Housekeeping note: it's unclear a priori how many terms need
// to be kept in the vector b.  Note that this function is integer:
// it returns 1 if the last few entries are less than 10^(-15), and
// zero otherwise.  This allows the calling function to determine
// if it is keeping enough terms.
//
int SWSH::findeig(Real *M, const int K, Real *E, Real b[])
{
  int i,j;
  
  gsl_matrix_view emm = gsl_matrix_view_array(M, K, K);
  gsl_vector *eigenvals = gsl_vector_alloc(K);
  gsl_matrix *eigenvecs = gsl_matrix_alloc(K, K);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(K);
  gsl_eigen_symmv(&emm.matrix, eigenvals, eigenvecs, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eigenvals, eigenvecs, GSL_EIGEN_SORT_ABS_ASC);
  //
  *E = -gsl_vector_get(eigenvals, Gl - lmin);
  for(i = 1; i <= K; i++) {
    b[i] = gsl_matrix_get(eigenvecs, i-1, Gl - lmin);
  }
  //
  gsl_vector_free(eigenvals);
  gsl_matrix_free(eigenvecs);
  if (fabs(b[K]) < 1.e-15 && fabs(b[K - 1]) < 1.e-15) {
    return(1);
  } else {
    return(0);
  }
}

Real SWSH::error(const Real E, const Real b[], const int n)
{
  Real *M, *berr;
  Real err = 0.0;
  int i, j;

  M = Tensor<Real>::vector(0, n*n - 1);
  load_M(M, n);
  berr = Tensor<Real>::vector(0, n - 1);
  for (i = 0; i < n; i++) {
    berr[i] = 0.0;
    for (j = 0; j < n; j++)
      berr[i] += M[j + i*n]*b[j];
    berr[i] += E*b[i];
  }
  for (i = 0; i < n; i++)
    err += berr[i]*berr[i];
  Tensor<Real>::free_vector(berr, 0, n - 1);
  Tensor<Real>::free_vector(M, 0, n*n - 1);
  return(sqrt(err)/E);
}

//
// Returns spheroid.
//
Real SWSH::spheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  y2 = 2.0*(y = x);
  for (j = N_spheroid; j >= 1; j--) {
    sv = d;
    d = y2*d - dd + spheroid_cofs[j];
    dd = sv;
  }
  if (Gm%2)
    ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*spheroid_cofs[0]);
  else
    ans = y*d - dd + 0.5*spheroid_cofs[0];

  if (fabs(ans) < 5.e-16) ans = 0.;
  return ans;
}

//
// Returns edthbar^2 on spheroid.  Only does
// spheroid of spin weight +2.
//
Complex SWSH::edthbaredthbarspheroid(const Real a, const Real x)
{
  Complex ans = Complex(0., 0.);

  Real rp = Kerr::rplus(a);
  Real OmegaH = a/(2.*rp);
  Real sinthsq = (1. - x)*(1. + x);
  Real sinth = sqrt(sinthsq);
  Complex rpmiact = rp - II*a*x;
  Complex A1, A2;

  int k;

  A1 = -2.*a*sinth*(Gm*OmegaH + 2.*II/rpmiact);
  A2 = a*a*sinthsq*(Gm*Gm*OmegaH*OmegaH + 4.*Gm*OmegaH*II/rpmiact
		    -6./(rpmiact*rpmiact));

  if (s == 2) {
    for (k = lmin; k <= lmin + N - 1; k++) {
      ans += b[k - lmin + 1]*sqrt((k + 2)*(k + 1)*k*(k - 1))*zeroY(k, Gm, x);
      ans += b[k - lmin + 1]*A1*sqrt((k + 2)*(k - 1))*pos1Y(k, Gm, x);
      ans += b[k - lmin + 1]*A2*pos2Y(k, Gm, x);
    }
    ans /= (2.*rpmiact*rpmiact);
  } else {
    // This function only works for the s == 2 case!
    // Do nothing if called for another spin weight.
  }
  return(ans);
}

//
// Returns operator L^\dag_2 on spheroid.
//
Real SWSH::l2dagspheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  if (s == -2) {
    y2 = 2.0*(y = x);
    for (j = N_l2dagspheroid; j >= 1; j--) {
      sv = d;
      d = y2*d - dd + l2dagspheroid_cofs[j];
      dd = sv;
    }
    if (Gm%2)
      ans = y*d - dd + 0.5*l2dagspheroid_cofs[0];
    else
      ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*l2dagspheroid_cofs[0]);
    if (fabs(ans) < 5.e-16) ans = 0.;
    return ans;
  } else return(0.);
}

//
// Returns operator L^\dag_1 L^\dag_2 on spheroid.
//
Real SWSH::l1dagl2dagspheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  if (s == -2) {
    y2 = 2.0*(y = x);
    for (j = N_l1dagl2dagspheroid; j >= 1; j--) {
      sv = d;
      d = y2*d - dd + l1dagl2dagspheroid_cofs[j];
      dd = sv;
    }
    if (Gm%2)
      ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*l1dagl2dagspheroid_cofs[0]);
    else
      ans = y*d - dd + 0.5*l1dagl2dagspheroid_cofs[0];
    if (fabs(ans) < 5.e-16) ans = 0.;
    return ans;
  } else return(0.);
}

void SWSH::cheb_spheroid()
{
  int k;
  Real y, fac;
  Real input[CHEBCOFS];

  fac = 1.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_spheroid(y)/sqrt((1. - y)*(1. + y));
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_spheroid(y);
    }
  }
  fftw_plan p;
  p = fftw_plan_r2r_1d(CHEBCOFS, input, spheroid_cofs, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(spheroid_cofs[k]) > 5.e-15) N_spheroid = k;
}

void SWSH::cheb_l2dagspheroid()
{
  int k;
  Real y, fac;
  Real input[CHEBCOFS];

  fac = 1.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_l2dagspheroid(y);
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_l2dagspheroid(y)/sqrt((1. - y)*(1. + y));
    }
  }
  fftw_plan p;
  p = fftw_plan_r2r_1d(CHEBCOFS, input, l2dagspheroid_cofs, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(l2dagspheroid_cofs[k]) > 5.e-15) N_l2dagspheroid = k;
}

void SWSH::cheb_l1dagl2dagspheroid()
{
  int k;
  Real y, fac;
  Real input[CHEBCOFS];

  fac = 1.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_l1dagl2dagspheroid(y)/sqrt((1. - y)*(1. + y));
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      input[k] = fac*raw_l1dagl2dagspheroid(y);
    }
  }
  fftw_plan p;
  p = fftw_plan_r2r_1d(CHEBCOFS, input, l1dagl2dagspheroid_cofs, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(l1dagl2dagspheroid_cofs[k]) > 5.e-15) N_l1dagl2dagspheroid = k;
}

Real SWSH::raw_spheroid(const Real x)
{
  int k;
  Real ans = 0.0;

  if (s != 2 && s!= 1 && s != 0 && s != -1 && s != -2) {
    cerr << "This code can only handle spin weights 2, 1, 0, -1, -2." << endl;
    exit(0);
  }

  for (k = lmin; k <= lmin + N - 1; k++) {
    if (s == 2)
      ans += b[k - lmin + 1]*pos2Y(k, Gm, x);
    else if (s == 1)
      ans += b[k - lmin + 1]*pos1Y(k, Gm, x);
    else if (s == 0)
      ans += b[k - lmin + 1]*zeroY(k, Gm, x);
    else if (s == -1)
      ans += b[k - lmin + 1]*neg1Y(k, Gm, x);
    else if (s == -2)
      ans += b[k - lmin + 1]*neg2Y(k, Gm, x);
  }
  return ans;
}

Real SWSH::raw_l2dagspheroid(const Real x)
{
  int k, s = -2;
  Real ans = 0.0;

  for (k = lmin; k <= lmin + N - 1; k++) {
    ans += aw*sqrt((1. - x)*(1. + x))*b[k - lmin + 1]*neg2Y(k, Gm, x);
    ans -= b[k - lmin + 1]*sqrt((k - 1)*(k + 2))*neg1Y(k, Gm, x);
  }

  return ans;
}

Real SWSH::raw_l1dagl2dagspheroid(const Real x)
{
  int k, s = -2;
  Real ans = 0.0;

  for (k = lmin; k <= lmin + N - 1; k++) {
    ans += sqrt((k - 1)*k*(k + 1)*(k + 2))*b[k - lmin + 1]*zeroY(k, Gm, x);
    ans -= aw*aw*(1. - x)*(1. + x)*b[k - lmin + 1]*neg2Y(k, Gm, x);
    //
    // These next two lines are the $2 aw sin(theta) L_2^\dag S$ term
    //
    ans += 2.0*aw*aw*(1. - x)*(1. + x)*b[k - lmin + 1]*neg2Y(k, Gm, x);
    ans -= 2.0*aw*sqrt((1. - x)*(1. + x))*b[k - lmin + 1]*
      sqrt((k - 1)*(k + 2))*neg1Y(k, Gm, x);
  }
  return ans;
}


