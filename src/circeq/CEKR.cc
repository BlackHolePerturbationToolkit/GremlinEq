//---------------------------------------------------------------------------
//
// $Id: CEKR.cc,v 1.21 2019/01/26 14:39:18 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods of class CEKR --- the class containing objects relevant
// to radiation for circular, equatorial orbits.
//
// Scott Hughes, 11 January 1999.
//
#include <cmath>
#include "Globals.h"
#include "CEKR.h"
#include "SWSH.h"
#include "FT.h"

#define SW -2

// Constructor.
CEKR::CEKR(SWSH *swsh_in, FT *ft_in, CEKG *cekg_in) :
  swsh(swsh_in), ft(ft_in), cekg(cekg_in)
{
  //
  // Local storage of orbital characteristics.
  //
  r = cekg->r;
  a = cekg->a;
  E = cekg->E;
  Lz = cekg->Lz;
  //
  // Local storage of mode characteristics.
  //
  l = ft->l;
  m = ft->m;
  wm = ft->omega;
  pm = wm - m*a/(2.*Kerr::rplus(a));
  //
  // Local storage of useful quantities that appear all over the
  // place.
  //
  Delta = Kerr::Delta(r, a);
  dr_Delta = Kerr::dr_Delta(r);
  ddr_Delta = Kerr::ddr_Delta();
  r_plus = Kerr::rplus(a);
  SN_Kay = ft->K(r);
  dr_SN_Kay = ft->dr_K(r);
  //
  // Compute amplitudes ZI, ZH.
  //
  ZI = ft->c_trans*Zed(ft->rteuk_in, ft->drteuk_in, ft->ddrteuk_in);
  ZH = ft->b_trans*Zed(ft->rteuk_up, ft->drteuk_up, ft->ddrteuk_up);
}

// The complex amplitudes ZH, ZI that are used to construct Edot, Lzdot,
// and the waveform.  To get ZH, pass in Rup and its derivatives; for ZI,
// pass in Rin and its derivatives.
Complex CEKR::Zed(const Complex R,
		  const Complex dr_R,
		  const Complex ddr_R)
{
  Complex term1 = R*(Ann0() + Anmbar0() + Ambarmbar0());
  Complex term2 = -dr_R*(Anmbar1() + Ambarmbar1());
  Complex term3 = ddr_R*Ambarmbar2();
  //
  Complex Wronskian;

  Real rp = Kerr::rplus(a);
  Real rm = Kerr::rminus(a);
  if (m == 0) {
    Wronskian = -2.*exp(gsl_sf_lngamma(2.*l + 2.) - gsl_sf_lngamma(l + 1.)
			- gsl_sf_lngamma(l + 3.));
    Wronskian *= pow(rp - rm, 2 - l);
  } else {
    Wronskian = 2.*II*wm*ft->b_inc*ft->c_trans;
  }
  Complex tmp;
  Real dfactor = 1.;
  //
  // We've been running into some underflow problems associated with
  // the denominator getting huge.  If it exceeds 1e200 or 1e100,
  // factor that out and reinsert later.
  if (abs(Wronskian) > 1.e200) {
    dfactor = 1.e-200;
    Wronskian *= dfactor;
  } else if (abs(Wronskian) > 1.e100) {
    dfactor = 1.e-100;
    Wronskian *= dfactor;
  }
  tmp = (term1 + term2 + term3)/Wronskian;
  tmp *= dfactor;

  if (isnan(Wronskian.real())) tmp = 0.;
  //
  return(-2.*M_PI*tmp);
}

// The Newman-Penrose quantity rho (with sign as in
// Teukolsky's papers).  Note this is the *equatorial*
// code, so the imaginary part is zero.  (Note added
// after nearly suffering a heart attack.)
Complex CEKR::rho_func()
{
  const Complex tmp = Complex(r, 0);

  return(-1./tmp);
}

// The C_{ab} functions are projections of the orbiting
// particle's stress-energy tensor onto the Newman-Penrose
// tetrad (modulo some factors).  Here they are, plus
// their derivatives.
Complex CEKR::Cnn()
{
  const Real tmp = E*(r*r + a*a) - a*Lz;
  const Real sig = Kerr::Sigma(r, a, 0);

  return(Complex(0.25*tmp*tmp/(sig*sig*gkg.TFunc(r, a, 0., E, Lz, 0.)), 0.));
}

Complex CEKR::Cnmbar()
{
  const Real st = 1.;
  const Complex rho = rho_func();
  const Real sig = Kerr::Sigma(r, a, 0.);

  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  const Complex tmp2 = Complex(sqrt(gkg.ThetaFunc(r, a, 0., E, Lz, 0.)),
			       st*a*E - Lz/st);
  const Complex tmp3 = tmp1*tmp2;

  return(rho*tmp3/(2.*sqrt(2.)*sig*gkg.TFunc(r, a, 0., E, Lz, 0.)));
}

Complex CEKR::Cmbarmbar()
{
  const Real st = 1.;
  const Complex rho = rho_func();

  const Complex tmp = Complex(sqrt(gkg.ThetaFunc(r, a, 0., E, Lz, 0.)),
			      a*E*st - Lz/st);

  return(rho*rho*tmp*tmp/(2.*gkg.TFunc(r, a, 0., E, Lz, 0.)));
}

Complex CEKR::Ann0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;

  const Complex tmp1 = 2.*II*a*rho*st*swsh->l2dagspheroid(0.);
  const Real tmp2 = swsh->l1dagl2dagspheroid(0.);

  const Complex pref = -2.*Cnn()/(Delta*Delta*rhob*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CEKR::Anmbar0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;
  const Real S = swsh->spheroid(0.);
  const Real L2dagS = swsh->l2dagspheroid(0.);

  const Complex tmp1 = (II*SN_Kay/Delta - rho - rhob)*L2dagS;
  const Complex tmp2 = (SN_Kay/Delta)*a*st*S*(rhob - rho);

  const Complex pref = -2.*sqrt(2.)*Cnmbar()/(Delta*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CEKR::Anmbar1()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;
  const Real S = swsh->spheroid(0.);
  const Real L2dagS = swsh->l2dagspheroid(0.);

  const Complex tmp = L2dagS + II*a*st*(rho - rhob)*S;

  const Complex pref = -2.*sqrt(2.)*Cnmbar()/(Delta*rho3);

  return(pref*tmp);
}
    
Complex CEKR::Ambarmbar0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  const Complex tmp1 = (SN_Kay/Delta)*(SN_Kay/Delta);
  const Complex tmp2 = 2.*rho*II*SN_Kay/Delta;
  const Complex tmp3 = II*(dr_SN_Kay - SN_Kay*dr_Delta/Delta)/Delta;

  const Complex pref = S*Cmbarmbar()*rhob/rho3;

  return(pref*(tmp1 + tmp2 + tmp3));
}
    
Complex CEKR::Ambarmbar1()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  const Complex tmp = -II*SN_Kay/Delta;

  const Complex pref = 2.*S*rhob*Cmbarmbar()/rho3;

  return(pref*(rho + tmp));
}

Complex CEKR::Ambarmbar2()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  return(-Cmbarmbar()*rhob*S/rho3);
}
