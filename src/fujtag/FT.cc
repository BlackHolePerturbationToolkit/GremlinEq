//---------------------------------------------------------------------------
//
// $Id: FT.cc,v 1.6 2014/04/17 00:07:00 sahughes Exp $
//
//---------------------------------------------------------------------------

/* FT.cc
 * A class for interacting with the Teukolsky equation solver described
 * by Fujita and Tagoshi
 */

#include <cstdio>
#include "FT.h"

/* try to get the best accuracy to within this factor */
#define ACCURACY_BOUND 10

#define NUM_UP_SOLVERS 2
#define NUM_IN_SOLVERS 3

#define MAX_NUM_SOLVERS \
  (NUM_UP_SOLVERS > NUM_IN_SOLVERS ? NUM_UP_SOLVERS : NUM_IN_SOLVERS)

#define xofr(r) ((1+(1-(r))/Sqrt(1-a*a))/2)
#define zofr(r) (omega*((r) - 1 + Sqrt(1-a*a)))

FT::FT(const int l, const int m, const Real r, const Real a, const Real omega,
       const Real lambda, const Real tolerance)
  : l(l), m(m), a(a), p(r), e(0), omega(omega), lambda(lambda),
    tolerance(tolerance) {
  accuracy_in = accuracy_up = tolerance;
  rin_request_precision = rup_request_precision = REAL_EPSILON;
  hypergeom_terms = 0;
#ifdef USE_GMP
  gmp_on = 1;
#else /* USE_GMP */
  gmp_on = 0;
#endif
  //
  // Code introduced by SAH 23 August 2013
  //
  if (m == 0) {
    Real rp = 1. + sqrt((1. - a)*(1. + a));
    Real rm = 1. - sqrt((1. - a)*(1. + a));
    Real rpmrm = rp - rm;
    Real x = (r - rp)/rpmrm;
    Real oox = 1./x;
    //
    // First three args of hypergeometric function for Rin.
    Real ain = 2. - l, bin = l + 3., cin = 3.;
    //
    // Rin depends on 2F1(ain, bin, cin, -x) ... use identities
    // to map this to a form with a smaller argument.
    Complex hyperin;
    hyperin = pow(1. + x, -ain)*hypergeom2F1(ain, cin - bin, cin, -x/(-x - 1.));
    //
    // Also need derivative:
    // d[2F1(ain, bin, cin, -x)]/dx = -(ain bin/cin)*
    //   2F1(ain + 1, bin + 1, cin + 1, -x).
    // Use identities to get this too.
    Complex dhyperin;
    dhyperin = -(ain*bin/cin)*pow(1. + x, -(ain + 1.))*
      hypergeom2F1(ain + 1, cin - bin, cin + 1., -x/(-x - 1.));

    rteuk_in = x*x*(1. + x)*(1. + x)*hyperin;
    rteuk_in *= rpmrm*rpmrm*rpmrm*rpmrm;
    drteuk_in = x*(1. + x)*((2.*x + 2*(1. + x))*hyperin +
			    x*(1. + x)*dhyperin);
    drteuk_in *= rpmrm*rpmrm*rpmrm;
    ddrteuk_in = get_ddrteuk(r, rteuk_in, drteuk_in);
    //
    // First three args of hypergeometric function for Rup.
    Real aup = l + 3., bup = l + 1., cup = 2.*l + 2.;
    //
    // Rup depends on 2F1(aup, bup, cup, -1/x) ... remap this
    // one too since 1/x gets large for orbits near the horizon.
    //
    Complex hyperup;
    hyperup = pow(1. + oox, -aup)*
      hypergeom2F1(aup, cup - bup, cup, oox/(oox + 1.));
    //
    // Also need derivatives; identity similar to the one above applies
    // here.
    //
    Complex dhyperup;
    dhyperup = (aup*bup/cup)*oox*oox*pow(1. + oox, -aup - 1.)*
      hypergeom2F1(aup + 1., cup - bup, cup + 1., oox/(oox + 1.));
    rteuk_up = pow(oox, l - 1.)*(1. + oox)*(1. + oox)*hyperup;
    rteuk_up *= pow(rpmrm, 1. - l);
    drteuk_up = (1. + oox)*pow(oox, l - 1.)*
      (((1. + oox)*oox*(1. - l) - 2.*oox*oox)*hyperup +
       (1. + oox)*dhyperup);
    drteuk_up *= pow(rpmrm, -l);
    ddrteuk_up = get_ddrteuk(r, rteuk_up, drteuk_up);
    //
    // These amplitudes are irrelevant for the m = 0 solutions,
    // set them to 1.
    b_inc = 1.;
    b_trans = 1.;
    c_trans = 1.;
    return;
  } else {
    for(int i=0;i<TEST_LENGTH;i++)
      test_renangmoms[i] = CNAN;
    
    if(!renangmom(2*omega, a, l, m, lambda, &nu)) {
      fprintf(stderr, "1: Couldn't find nu at l=%d m=%d omega=%.17f\n",
	      l, m, omega);
      // throw "Couldn't find nu!";
      rteuk_in = 0.;
      drteuk_in = 0.;
      ddrteuk_in = 0.;
      rteuk_up = 0.;
      drteuk_up = 0.;
      ddrteuk_up = 0.;
      b_inc = 1.;
      b_trans = 1.;
      c_trans = 1.;
    } else {
      asympt_amps(nu, 2*omega, a, m, lambda, &b_trans, &b_inc, NULL, &c_trans);
      
      int numsolvers = MAX_NUM_SOLVERS;
      Real *errs = new Real[numsolvers];
      Real minerr;
      
      minerr = INFINITY;
      geterrs(p, errs, 0);
      //
      // Modified by SAH 2 Aug 2012: Make sure up_choices[0] and
      // in_choices[0] get set, even if the error condition is not met.
      //
      up_choices[0] = 0;
      for(int n=0;n<NUM_UP_SOLVERS;n++) {
	if(errs[n] < minerr) {
	  up_choices[0] = n;
	  minerr = errs[n];
	}
      }
      up_choice_locs[0] = INFINITY;
      accuracy_up = minerr;
      
      minerr = INFINITY;
      geterrs(p, errs, 1);
      in_choices[0] = 0;
      for(int n=0;n<NUM_IN_SOLVERS;n++) {
	if(errs[n] < minerr) {
	  in_choices[0] = n;
	  minerr = errs[n];
	}
      }
      in_choice_locs[0] = INFINITY;
      accuracy_in = minerr;
      
      delete[] errs;
      
      CalcRFields(r, 0);
      CalcRFields(r, 1);
    }
  }
}

FT::FT(const int l, const int m, const Real p, const Real e,
       const Real a, const Real omega, const Real lambda, const Real tolerance)
  : l(l), m(m), a(a), p(p), e(e), omega(omega), lambda(lambda),
    tolerance(tolerance) {
  accuracy_in = accuracy_up = tolerance;
  rin_request_precision = rup_request_precision = REAL_EPSILON;
  hypergeom_terms = 0;
#ifdef USE_GMP
  gmp_on = 1;
#else /* USE_GMP */
  gmp_on = 0;
#endif
  for(int i=0;i<TEST_LENGTH;i++)
    test_renangmoms[i] = CNAN;

  if(!renangmom(2*omega, a, l, m, lambda, &nu)) {
    fprintf(stderr, "2: Couldn't find nu at l=%d m=%d omega=%.17f\n",
	    l, m, omega);
    //    throw "Couldn't find nu!";
  }
  asympt_amps(nu, 2*omega, a, m, lambda, &b_trans, &b_inc, NULL, &c_trans);

  rteuk_in = drteuk_in = ddrteuk_in = 0;
  rteuk_up = drteuk_up = ddrteuk_up = 0;

  up_choice_locs[0] = -1;
  in_choice_locs[0] = -1;
}

void FT::CalcRFields(const Real rad, const int DoH) {
  if(DoH) {
    rteuk_in = callin(rad, &drteuk_in);
    ddrteuk_in = get_ddrteuk(rad, rteuk_in, drteuk_in);
  } else {
    rteuk_up = callup(rad, &drteuk_up);
    ddrteuk_up = get_ddrteuk(rad, rteuk_up, drteuk_up);
  }
}

inline Complex FT::get_ddrteuk(const Real rad,
			       const Complex rteuk, const Complex drteuk) {
  return (teuk_potential(rad) * rteuk + 2*(rad-1)*drteuk)
    / (rad*rad - 2*rad + a*a);
}

/* choose which solvers to use for different parts of the orbit */
int FT::choose_solvers(const Real low, const Real high,
		   const Real *lowerrs, const Real *higherrs,
		   const int isin, const int oldchoice,
		   Real **locptr, int **choiceptr, const int *maxchoiceptr) {
  int choice = -1;
  Real choiceerr = isin ? accuracy_in : accuracy_up;

  int numsolvers = isin ? NUM_IN_SOLVERS : NUM_UP_SOLVERS;
  for(int n=0;n<numsolvers;n++) {
    if(lowerrs[n] < choiceerr && higherrs[n] < choiceerr) {
      choice = n;
      choiceerr = Fmax(lowerrs[n],higherrs[n]);
    }
  }

  if(choice != -1) {
    /* we've found an acceptable solver */
    if(choice != oldchoice && oldchoice != -1) {
      if(*choiceptr >= maxchoiceptr) {
	throw("Too many solver divisions in FT");
      }
      *((*choiceptr)++) = oldchoice;
      *((*locptr)++) = low;
    }
    return choice;
  }

  /* no good choices, need to subdivide */
  Real med = (low + high)/2;
  Real *mederrs = new Real[numsolvers];
  geterrs(med, mederrs, isin);

  choice = choose_solvers(low, med, lowerrs, mederrs, isin, oldchoice,
			  locptr, choiceptr, maxchoiceptr);
  choice = choose_solvers(med, high, mederrs, higherrs, isin, choice,
			  locptr, choiceptr, maxchoiceptr);
  delete[] mederrs;
  return choice;
}

/* estimate the numerical error in the different solvers */
void FT::geterrs(const Real r, Real *errs, int isin) {
  int numsolvers = isin ? NUM_IN_SOLVERS : NUM_UP_SOLVERS;
  const Real testpattern_r[]   = { 1+REAL_EPSILON, 1-REAL_EPSILON, 1, 1 };
  const Real testpattern_eps[] = { 1, 1, 1+REAL_EPSILON, 1-REAL_EPSILON };

  for(int n=0;n<numsolvers;n++) {
    Complex deriv;
    Complex val = callsolver(isin, n, r, 2*omega, nu, REAL_EPSILON, &deriv);
    if(Cisnan(val) || Cisnan(deriv)) {
      errs[n] = NAN;
      continue;
    }
    errs[n] = 0;
    for(int test=0;test<TEST_LENGTH;test++) {
      Complex tmpnu;
      if(Cisnan(test_renangmoms[test])) {
	if(testpattern_eps[test] == 1 ||
	   !renangmom(2*omega*testpattern_eps[test], a, l, m, lambda,
		      &tmpnu)) {
	  tmpnu = nu;
	}
	test_renangmoms[test] = tmpnu;
      } else {
	tmpnu = test_renangmoms[test];
      }
      Complex epsderiv;
      Complex epsval = callsolver(isin, n, r*testpattern_r[test],
				  2*omega*testpattern_eps[test], tmpnu,
				  REAL_EPSILON, &epsderiv);
      if(Cisnan(epsval) || Cisnan(epsderiv)) {
	errs[n] = NAN;
	break;
      }
      Real err = Fmax(abs(epsval/val-1),abs(epsderiv/deriv-1));
      if(err > errs[n])
	errs[n] = err;
    }
  }
  Real *accuracy_var = isin ? &accuracy_in : &accuracy_up;
  Real minerr = INFINITY;
  for(int n=0;n<numsolvers;n++) {
    if(errs[n] < *accuracy_var)
      return;
    if(errs[n] < minerr)
      minerr = errs[n];
  }
  *accuracy_var = ACCURACY_BOUND * minerr;
}

Complex FT::callsolver(const int isin, const int which, const Real r,
		       const Real epsilon, const Complex nu,
		       const Real precision, Complex *deriv) {
  if(isin) {
    switch(which) {
    case 0:
      return rin_small(nu, epsilon, a, m, lambda, xofr(r),
		       precision, deriv);
    case 1:
      return rin_coulomb(nu, epsilon, a, m, lambda, zofr(r),
			 precision, deriv);
    case 2:
      return rin_hyper(nu, epsilon, a, m, lambda, xofr(r),
		       precision, deriv);
    }
  } else {
    switch(which) {
    case 0:
      return rup_tricomi(nu, epsilon, a, m, lambda, zofr(r),
			 precision, deriv);
    case 1:
      return rup_hyper(nu, epsilon, a, m, lambda, xofr(r),
		       precision, deriv);
    }
  }
  throw("Bad solver in FT");
}

Complex FT::callin(const Real r, Complex *deriv) {
  int choicenum = 0;
  if(in_choice_locs[0] == -1) getlocs(1);
  while(in_choice_locs[choicenum] < r)
    choicenum++;
  return callsolver(1, in_choices[choicenum], r, 2*omega, nu,
		    rin_request_precision, deriv);
}

Complex FT::callup(const Real r, Complex *deriv) {
  int choicenum = 0;
  if(up_choice_locs[0] == -1) getlocs(0);
  while(up_choice_locs[choicenum] < r)
    choicenum++;
  return callsolver(0, up_choices[choicenum], r, 2*omega, nu,
		    rup_request_precision, deriv);
}

void FT::getlocs(const int isin) {
  Real r_peri = p/(1+e);
  Real r_ap = p/(1-e);

  int choice;
  Real *errs_peri = new Real[MAX_NUM_SOLVERS];
  Real *errs_ap = new Real[MAX_NUM_SOLVERS];
  Real *locptr;
  int *choiceptr;

  if(isin) {
    locptr = in_choice_locs;
    choiceptr = in_choices;
  } else {
    locptr = up_choice_locs;
    choiceptr = up_choices;
  }
  geterrs(r_peri, errs_peri, isin);
  geterrs(r_ap, errs_ap, isin);
  choice = choose_solvers(r_peri, r_ap, errs_peri, errs_ap, isin, -1,
			  &locptr, &choiceptr, choiceptr + (FT_MAXDIVS-1));
  *choiceptr = choice;
  *locptr = INFINITY;

  delete[] errs_peri;
  delete[] errs_ap;
}

void FT::set_gmp(int flag) {
#if USE_GMP
  if(gmp_on == flag)
    return;
  gmp_on = flag;

  accuracy_in = accuracy_up = tolerance;
  rin_request_precision = rup_request_precision = REAL_EPSILON;
  in_choice_locs[0] = up_choice_locs[0] = -1;
#else /* USE_GMP */
  throw "gmp switching not enabled!";
#endif /* USE_GMP */
}
