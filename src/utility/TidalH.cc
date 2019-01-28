//---------------------------------------------------------------------------
//
// $Id: TidalH.cc,v 1.13 2019/01/27 18:24:36 sahughes Exp $
//
//---------------------------------------------------------------------------
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include "Globals.h"
#include "Tensors.h"
#include "SWSH.h"
#include "TidalH.h"

#define JMAX 15
#define JMAXP (JMAX+1)
#define K 5
#define TINY 1.0e-25
#define EPS 2.e-15

TidalH::TidalH(const Real spin, const int ellmax)
{
  a = spin;
  rp = Kerr::rplus(a);
  const Real oma2 = (1. - a)*(1. + a);
  eps = sqrt(oma2)/(4.*rp);
  if (fabs(a) < 1.e-10)
    Kph = -a/2.;
  else if (fabs(a) < 1.e-8)
    Kph = -a/2*(1. + (2.*log(a/2.) - 1)*a*a/4.);
  else
    Kph = 0.5*a/(rp - a*a)*
      (a*a - rp + 2.*atanh(sqrt(oma2)) + sqrt(oma2)*log(0.25*a*a/oma2));
  lmax = ellmax;
  ZH = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);
  Clm = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);
  swshp = Tensor<SWSH>::matrix(2, lmax, -lmax, lmax);
  w = Tensor<Real>::vector(-lmax, lmax);
  p = Tensor<Real>::vector(-lmax, lmax);
  lamb = Tensor<Real>::matrix(2, lmax, -lmax, lmax);
  Epslm = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);
  int lc, mc;
  for (lc = 2; lc <= lmax; lc++) {
    for (mc = -lmax; mc <= lmax; mc++) {
      Epslm[lc][mc] = 0.;
    }
  }
}

TidalH::~TidalH()
{
  Tensor<Complex>::free_matrix(ZH, 2, lmax, -lmax, lmax);
  Tensor<Complex>::free_matrix(Clm, 2, lmax, -lmax, lmax);
  Tensor<SWSH>::free_matrix(swshp, 2, lmax, -lmax, lmax);
  Tensor<Real>::free_vector(w, -lmax, lmax);
  Tensor<Real>::free_vector(p, -lmax, lmax);
  Tensor<Real>::free_matrix(lamb, 2, lmax, -lmax, lmax);
  Tensor<Complex>::free_matrix(Epslm, 2, lmax, -lmax, lmax);
}

//
// The following function assumes that w[m] and lamb[l][m] have been
// loaded.  Segfaults will result if used at the wrong moment.
//
void TidalH::loadClm()
{
  int l, m;

  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      const Real lp2 = lamb[l][m] + 2.;

      const Real tmp1 = lp2*lp2 + 4.*a*w[m]*(((Real)m) - a*w[m]);
      const Real tmp2 = lamb[l][m]*lamb[l][m] + 36.*a*w[m]*(((Real)m) - a*w[m]);
      const Real tmp3 = 48.*a*w[m]*(2.*lamb[l][m] + 3.)*(2.*a*w[m] - ((Real)m));
      Real abs_sqr_clm = tmp1*tmp2 + tmp3 + 144.*w[m]*w[m]*(1. - a*a);
      const Real sign = ((l + m)%2 ? -1. : 1.);
      const Real Im_clm = sign*12.*w[m];
      const Real Re_clm = sqrt(abs_sqr_clm - 144.*w[m]*w[m]);
      const Complex clm = Re_clm + II*Im_clm;

      Clm[l][m] = 256.*rp*rp*(p[m] + 4.*II*eps)*(II*p[m] - 2.*eps)/clm;
    }
  }
}

void TidalH::loadEpslm()
{
  Real *mDqell;
  Real *mRq_re, *mRq_im;
  int q, ell, m, ellmin;
  //
  // Notice some of the shenanigans with offset.  This is done in
  // order for arrays that should be on the range (lmin, lmax) map
  // to (0, lmax - lmin) in order to work nicely with the linear
  // algebra package.
  //
  int n, offset;
  Real d;
  int *indx;
  //
  for (m = -lmax; m <= lmax; m++) {
    ellmin = Max(2, abs(m));
    offset = ellmin - 1;
    n = lmax - offset;
    mRq_re = Tensor<Real>::vector(0, n - 1);
    mRq_im = Tensor<Real>::vector(0, n - 1);
    mDqell = Tensor<Real>::vector(0, n*n - 1);
    for (q = ellmin; q <= lmax; q++) {
      const Complex tmp = mR_vector(q, m);
      mRq_re[q-ellmin] = tmp.real();
      mRq_im[q-ellmin] = tmp.imag();
      for (ell = ellmin; ell <= lmax; ell++)
 	mDqell[(ell - ellmin) + (q - ellmin)*n] = mD_matrix(q, ell, m);
    }
    gsl_matrix_view emm = gsl_matrix_view_array(mDqell, n, n);
    gsl_vector_view mRq_re_vec = gsl_vector_view_array(mRq_re, n);
    gsl_vector_view mRq_im_vec = gsl_vector_view_array(mRq_im, n);
    gsl_vector *solre = gsl_vector_alloc(n);
    gsl_vector *solim = gsl_vector_alloc(n);
    int s;
    gsl_permutation *p = gsl_permutation_alloc(n);

    gsl_linalg_LU_decomp(&emm.matrix, p, &s);
    gsl_linalg_LU_solve (&emm.matrix, p, &mRq_re_vec.vector, solre);
    gsl_linalg_LU_solve (&emm.matrix, p, &mRq_im_vec.vector, solim);

    for (ell = ellmin; ell <= lmax; ell++)
      Epslm[ell][m] = gsl_vector_get(solre, ell - ellmin) +
	II*gsl_vector_get(solim, ell - ellmin);

    gsl_permutation_free(p);
    gsl_vector_free(solre);
    gsl_vector_free(solim);
    Tensor<Real>::free_vector(mDqell, 0, n*n - 1);
    Tensor<Real>::free_vector(mRq_im, 0, n - 1);
    Tensor<Real>::free_vector(mRq_re, 0, n - 1);
  }
}

Real TidalH::R1lm(const int l, const int m, const Real x, const Real psi)
{
  Real Rm = (Real)m;
  Complex phase = exp(II*Rm*(psi - Kph));

  Complex ans = Clm[l][m]*ZH[l][m]*phase*
    swshp[l][m].edthbaredthbarspheroid(a, x);

  return(ans.imag());
}

Real TidalH::epsr(const int l, const int m, const Real x, const Real psi)
{
  Real Rm = (Real)m;
  Complex phase = exp(II*Rm*(psi - Kph));
  const Complex tmp = phase*Epslm[l][m];

  return(rp*tmp.imag()*swshp[l][m].zeroY(l, m, x));
}

Complex TidalH::mR_vector(const int q, const int m)
{
  Complex answer = 0.;
  if (q >= m) {
    Gq = q; Gm = m;
    const int lmin = Max(2, abs(Gm));
    INTEGRAND = mRq_IGRND;
    for (Gl = lmin; Gl <= lmax; Gl++) {
      Real igralre, igralim;
      gsl_integration_romberg_workspace *wre = gsl_integration_romberg_alloc(12);
      gsl_integration_romberg_workspace *wim = gsl_integration_romberg_alloc(12);
      size_t nevals;
      gsl_function Igrndre, Igrndim;
      Igrndre.function = &TidalH::IGRNDre_wrapper;
      Igrndre.params = this;
      Igrndim.function = &TidalH::IGRNDim_wrapper;
      Igrndim.params = this;
      gsl_integration_romberg(&Igrndre, 0, 1. - EPS, 0, 1e-14, &igralre, &nevals, wre);
      if (fabs(igralre) <= 1.e-14) igralre = 0.;

      gsl_integration_romberg(&Igrndim, 0, 1. - EPS, 0, 1e-14, &igralim, &nevals, wim);
      if (fabs(igralim) <= 1.e-14) igralim = 0.;

      answer += 2.*M_PI*Clm[Gl][Gm]*ZH[Gl][Gm]*(igralre + II*igralim);
      gsl_integration_romberg_free(wre);
      gsl_integration_romberg_free(wim);
    }
  }
  return(answer);
}

Real TidalH::mD_matrix(const int q, const int ell, const int m)
{
  Real answer = 0.;
  //
  // Only compute answer if k and ell are both greater or equal to m;
  // the harmonics are zero otherwise.
  //
  if (q >= abs(m) && ell >= abs(m)) {
    //
    // Only compute answer if q and ell are either both even or both odd;
    // the parity of C and D are such that the answer is zero otherwise.
    //
    if (!((q-ell)%2)) {
      Gq = q; Gl = ell; Gm = m;
    
      Real igralC, igralD;
      gsl_integration_romberg_workspace *wC = gsl_integration_romberg_alloc(12);
      gsl_integration_romberg_workspace *wD = gsl_integration_romberg_alloc(12);
      size_t nevals;
      gsl_function IgrndC, IgrndD;
      //
      INTEGRAND = mD_C_IGRND;
      IgrndC.function = &TidalH::IGRNDre_wrapper;
      IgrndC.params = this;
      gsl_integration_romberg(&IgrndC, 0, 1., 0, 1e-14, &igralC, &nevals, wC);
      if (fabs(igralC) <= 1.e-14)
	igralC = 0.;
      answer = 2.*M_PI*igralC;
      gsl_integration_romberg_free(wC);
      //
      INTEGRAND = mD_D_IGRND;
      IgrndD.function = &TidalH::IGRNDre_wrapper;
      IgrndD.params = this;
      gsl_integration_romberg(&IgrndD, 0, 1. - EPS, 0, 1e-14, &igralD, &nevals, wD);
      if (fabs(igralD) <= 1.e-14)
	igralD = 0.;
      answer += 2.*M_PI*igralD;
      gsl_integration_romberg_free(wD);
    }
  }
  return(answer);
}

Real TidalH::C0(const Real x)
{
  const Real x2 = x*x;
  const Real x3 = x*x2;
  const Real x4 = x*x3;
  const Real x5 = x*x4;
  const Real x6 = x*x5;
  const Real x7 = x*x6;
  const Real x8 = x*x7;

  const Real rp2 = rp*rp;
  const Real rp3 = rp*rp2;
  const Real rp4 = rp*rp3;
  const Real rp5 = rp*rp4;
  const Real rp6 = rp*rp5;

  const Real eig = ((Real)(Gl*(Gl+1)));
  if (fabs(a) > sqrt(3.)/2.) {
    cerr << "This code is only code for a <= sqrt(3)/2" << endl;
    exit(0);
  }
  const Real H = sqrt((rp4*rp4) - 6.*a*a*a*a*rp4*x2 -
		      4.*a*a*a*a*a*a*rp2*x2*(1. + x2) -
		      a*a*a*a*a*a*a*a*x2*(1. + x2 + x4));

  const Real c16 = x4*(-10. - x2 - x4 + 6.*x6);
  const Real c14 = rp2*x2*(3. - 71.*x2 - 8.*x4 + 18.*x6 + 18.*x8);
  const Real c12 = rp4*x2*(21. - 217.*x2 + 6.*x4 + 60.*x6 + 18.*x8);
  const Real c10 = rp*x2*(2.*H*x6*(3. + eig*x2) +
			     rp5*(63. - 355.*x2 + 56.*x4 + 62.*x6 + 6.*x8));
  const Real c8 = rp2*x2*(36.*H*x2 +
			     2.*H*rp*x4*(8. + (4. + 5.*eig)*x2) +
			     rp6*(104. - 332.*x2 + 67.*x4 + 21.*x6));
  const Real c6 = rp4*(4.*H*x2*(-4. + 27.*x2) +
			  4.*H*rp*x4*(3. + (6. + 5.*eig)*x2) +
			  rp6*(-2. + 103.*x2 - 181.*x4 + 24.*x6));
  const Real c4 = rp6*(4.*(6. + 5.*eig)*H*rp*x4 +
			  12.*H*x2*(-4. + 9.*x2) +
			  rp6*(-6. + 63.*x2 - 57.*x4));
  const Real c2 = 8.*eig*H*H*rp3*x4 + (rp6*rp6*rp2)*(-6. + 23.*x2 - 9.*x4)
    + 2.*H*rp5*(-16.*eig*x4 + 6.*rp3*x2*(-4. + 3.*x2)
		+ rp4*(-1. + (4. + 5.*eig)*x2));
  const Real c0 = 2.*rp5*(-(rp5*rp6) + 2.*rp5*(-4.*H + rp6)*x2
			     + eig*H*(rp6 + 4.*(H - 4.*rp2)*x2));

  const Real numer = c0 + a*a*(c2 +
				  a*a*(c4 +
				       a*a*(c6 +
					    a*a*(c8 +
						 a*a*(c10 +
						      a*a*(c12 +
							   a*a*(c14 +
								a*a*c16)))))));

  const Real denom = 2.*H*rp*pow(rp2 + a*a*x*x, 5.5);

  return(numer/denom);
}

Real TidalH::C1(const Real x)
{
  const Real omx2 = (1. - x)*(1. + x);
  const Real x2 = x*x;
  const Real x3 = x*x2;
  const Real x4 = x*x3;
  const Real x5 = x*x4;
  const Real x6 = x*x5;
  const Real x7 = x*x6;
  const Real x8 = x*x7;

  const int msq = Gm*Gm;

  const Real rp2 = rp*rp;
  const Real rp3 = rp*rp2;
  const Real rp4 = rp*rp3;
  const Real rp5 = rp*rp4;
  const Real rp6 = rp*rp5;
  const Real rp7 = rp*rp6;
  const Real rp8 = rp*rp7;
  const Real rp9 = rp*rp8;

  if (fabs(a) > sqrt(3.)/2.) {
    cerr << "This code is only code for a <= sqrt(3)/2" << endl;
    exit(0);
  }
  const Real H = sqrt((rp4*rp4) - 6.*a*a*a*a*rp4*x2 -
		      4.*a*a*a*a*a*a*rp2*x2*(1. + x2) -
		      a*a*a*a*a*a*a*a*x2*(1. + x2 + x4));

  Real numer;
  if (omx2 < 1.e-9) {
    numer = 128.*a*a*msq*(64.*rp -
			  a*a*(16.*(3.*rp + 2.) +
			       a*a*(8.*(rp - 2.) -
				    a*a*(5.*rp + 6. -
					 a*a))));
  } else {
    const Real c10 = -2.*H*msq*rp*x8*x2;
    const Real c8 = -2.*msq*rp2*x8*(H*(6. + 5.*rp) -
				       6.*(H - 4.*rp2)*x2);
    const Real c6 = -4.*msq*rp4*x6*(H*(8. + 5.*rp) -
				       8.*(H - 4.*rp2)*x2);
    const Real c4 = -4.*msq*rp6*x4*(H*(6. + 5.*rp) -
				       6.*(H - 4.*rp2)*x2);
    const Real c2 = -2.*H*msq*rp3*x2*(5*rp6 + 4.*(H - 4.*rp2)*x2);
    const Real c0 = -2.*msq*rp3*(-2.*H*rp7 + H*rp8 + 4.*H*H*rp2*x2 -
				    16.*H*rp4*x2 + 2.*H*rp7*x2 - 8.*rp9*x2);

    numer = c0 + a*a*(c2 +
		      a*a*(c4 + 
			   a*a*(c6 +
				a*a*(c8 +
				     a*a*c10))));
    numer /= omx2;
  }
  const Real denom = 2.*H*rp*pow(rp2 + a*a*x*x, 5.5);

  return(numer/denom);
}

Real TidalH::D(const Real x)
{
  const Real x2 = x*x;
  const Real x3 = x*x2;
  const Real x4 = x*x3;
  const Real x5 = x*x4;
  const Real x6 = x*x5;
  const Real x7 = x*x6;
  const Real x8 = x*x7;
  const Real x9 = x*x8;

  const Real rp2 = rp*rp;
  const Real rp3 = rp*rp2;
  const Real rp4 = rp*rp3;
  const Real rp5 = rp*rp4;
  const Real rp6 = rp*rp5;
  const Real rp7 = rp*rp6;
  const Real rp8 = rp*rp7;

  if (fabs(a) > sqrt(3.)/2.) {
    cerr << "This code is only code for a <= sqrt(3)/2" << endl;
    exit(0);
  }
  const Real H = sqrt((rp4*rp4) - 6.*a*a*a*a*rp4*x2 -
		      4.*a*a*a*a*a*a*rp2*x2*(1. + x2) -
		      a*a*a*a*a*a*a*a*x2*(1. + x2 + x4));

  const Real c14 = x5*(-8. + x2 + x4 + 6.*x6);
  const Real c12 = rp2*x3*(5. - 47.*x2 + 7.*x4 + 23.*x6 + 12.*x8);
  const Real c10 = rp4*x3*(30. - 114.*x2 + 35.*x4 + 43.*x6 + 6.*x8);
  const Real c8 = 2.*H*rp*x9 + rp6*x3*(75. - 149.*x2 + 53.*x4 + 21.*x6);
  const Real c6 = rp2*x3*(rp6*(89. - 113.*x2 + 24.*x4) +
			     4.*H*(2.*rp*x4 -10. + 9.*x2));
  const Real c4 = -4.*rp4*x*(rp6*(1. - 13.*x2 + 12.*x4) -
				3.*H*(1. - 8.*x2 + (6. + rp)*x4));
  const Real c2 = rp6*x*(rp6*(-8. + 17.*x2 - 9.*x4) +
			    4.*H*(6. + 2.*(-9. + rp)*x2 + 9.*x4));
  const Real c0 = 2.*rp8*x*(2.*rp6*(-1. + x2) + H*(rp + (6. - 8.*x2)));

  const Real numer = c0 + a*a*(c2 +
				  a*a*(c4 +
				       a*a*(c6 +
					    a*a*(c8 +
						 a*a*(c10 +
						      a*a*(c12 +
							   a*a*c14))))));
  const Real denom = H*pow(rp2 + a*a*x*x, 5.5);

  return(numer/denom);
}

Complex TidalH::IGRND(const Real x)
{
  Complex igrnd;

  if (INTEGRAND == mRq_IGRND) {
    igrnd = (swshp[Gl][Gm].edthbaredthbarspheroid(a, x)*
	     swshp[Gl][Gm].zeroY(Gq, Gm, x)) +
      (swshp[Gl][Gm].edthbaredthbarspheroid(a, -x)*
       swshp[Gl][Gm].zeroY(Gq, Gm, -x));
  } else if (INTEGRAND == mD_C_IGRND) {
    igrnd = swshp[Gl][Gm].zeroY(Gq, Gm, x)*(C0(x) + C1(x))*
      swshp[Gl][Gm].zeroY(Gl, Gm, x) +
      swshp[Gl][Gm].zeroY(Gq, Gm, -x)*(C0(-x) + C1(-x))*
      swshp[Gl][Gm].zeroY(Gl, Gm, -x);
  } else if (INTEGRAND == mD_D_IGRND) {
    igrnd = swshp[Gl][Gm].zeroY(Gq, Gm, x)*D(x)*
      swshp[Gl][Gm].dzeroYdx(Gl, Gm, x) +
      swshp[Gl][Gm].zeroY(Gq, Gm, -x)*D(-x)*
      swshp[Gl][Gm].dzeroYdx(Gl, Gm, -x);
  } else {
    igrnd = 0.;
  }
  return igrnd;
}
