//---------------------------------------------------------------------------
//
// $Id: fsum.cc,v 1.1 2013/05/17 15:39:53 sahughes Exp $
//
//---------------------------------------------------------------------------

/* fsum.cc
 * Perform sums of the expansion coefficients of the radial Teukolsky equation
 */

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "fsum.h"
#  include "fractions.h"
#  include "specialfuncs.h"
#endif
#include "fraction_macros.h"

/* the number of coefficients to generate in a block
 * setting this to 1 triggers some bug in the gcc optimization code as of
 * gcc 4.1.2 - not sure exactly which flag sets it off, but it works with no
 * optimization */
#define COEFFICIENT_BLOCK_SIZE 8

/* the absolute value of the sum of certain terms in the loop code must
 * exceed this fraction of the terms to not trigger a full recalculation of
 * the continued fraction */
#define LOOP_CANCEL_TOL 1e-1

/* the number of consecutive terms that must leave the loop sum stationary
 * before the loop exits */
#define LOOP_STAT_COUNT 10

/* sanity test for runaway sums */
#define LOOP_MAX 1000

/* the number of values of kfactor() to cache
 * kfactor() tends to be called in pairs, so there is probably no benefit
 * from cacheing if this is less than 2 */
#define K_CACHE_SIZE 2

/* this should be a list of K_CACHE_SIZE zeros */
#define K_CACHE_INIT { 0, 0 }

void FT_PRFX asympt_amps(Complex nu, Real epsilon, Real q, int m, Real lambda,
		     Complex *b_trans, Complex *b_inc, Complex *b_ref,
		     Complex *c_trans) {
  const Real kappa = Sqrt(1-q*q);
  const Real tau = (epsilon - m*q) / kappa;

  const Complex ieloge = I*epsilon*log((Complex)epsilon);

  Complex straightsum;
  Complex k_nu;
  Complex k__nu__1;
  Complex aminussum;

  Complex dummy;
  if(b_ref && !c_trans)
    c_trans = &dummy;
  if(b_ref || b_inc) {
    k_nu = kfactor(nu, epsilon, q, m, lambda);
    k__nu__1 = kfactor(-nu-1.0, epsilon, q, m, lambda);
  }
  if(c_trans)
    aminussum = aminus(nu, epsilon, q, m, lambda);
  if(b_trans || b_inc)
    straightsum = fsum(nu, epsilon, q, m, lambda);

  /* The expressions in Fujuta and Tagoshi, and also in Mino et al. are
   * incorrect. */
  //*b_trans = pow(omega / (epsilon * kappa), 4)
  //  * exp(I*0.5*(epsilon+tau)*log(kappa)) * straightsum;
  if(b_trans)
    *b_trans = 1 / pow(2 * kappa,4)
      * exp(I*(epsilon+tau)*kappa*(0.5 + log(kappa)/(1+kappa))) * straightsum;

  //*b_inc = (k_nu - I*exp(-I*M_PI*nu) * sin(M_PI*(nu + 2.0 + I*epsilon))
  //    / sin(M_PI*(nu - 2.0 - I*epsilon)) * k__nu__1) / omega 
  //  * pow(2.0,-3.0 - I*epsilon)
  //  * exp(M_PI*0.5*(-epsilon+I*(nu+3.0)) - ieloge
  //  + gammln(nu + 3.0 + I*epsilon) - gammln(nu - 1.0 - I*epsilon))
  //  * straightsum;
  if(b_inc)
    *b_inc = 2.0 * (k_nu - I*exp(-I*M_PI*nu) * sin(M_PI*(nu + 2.0 + I*epsilon))
		  / sin(M_PI*(nu - 2.0 - I*epsilon)) * k__nu__1) / epsilon 
      * pow(2.0,-3.0 - I*epsilon)
      * exp(M_PI*0.5*(-epsilon+I*(nu+3.0)) - ieloge + 0.5*I*(1-kappa)*epsilon
	    + gammln(nu + 3.0 + I*epsilon) - gammln(nu - 1.0 - I*epsilon))
      * straightsum;

  //*c_trans = omega * omega * omega * aminussum * exp(ieloge);
  if(c_trans)
    *c_trans = epsilon * epsilon * epsilon / 8 * aminussum
      * exp(ieloge - 0.5*I*(1-kappa)*epsilon);

  /* b_ref is correct *relative* to c_trans */
  if(b_ref)
    *b_ref = *c_trans * (k_nu + I*exp(I*M_PI*nu)*k__nu__1);
}

#define LOOP_POSITIVE(commands)						\
  do {									\
    lastn = 0;								\
    int stattimes = 0;							\
    for(;;) {								\
      fracs[COEFFICIENT_BLOCK_SIZE-1] =					\
	cplusfrac(lastn+COEFFICIENT_BLOCK_SIZE-1+nu,			\
		  epsilon, q, m, lambda, NULL);				\
      for(n=COEFFICIENT_BLOCK_SIZE-2;n>=0;n--) {			\
	Complex LOOP_denom = beta(lastn+n+1) + fracs[n+1];		\
	if(norm(LOOP_denom) >						\
	   norm(fracs[n+1])*(LOOP_CANCEL_TOL*LOOP_CANCEL_TOL)) {	\
	  fracs[n] = - alphagamma(lastn+n) / LOOP_denom;		\
	} else {							\
	  fracs[n] = cplusfrac(lastn+n+nu, epsilon, q, m, lambda, NULL); \
	}								\
      }									\
      									\
      for(n=0;n<COEFFICIENT_BLOCK_SIZE;n++,lastn++) {			\
	lastcoeff *= fracs[n]/alpha(lastn);				\
	commands;							\
      }									\
      if(lastn > LOOP_MAX)						\
	return NAN;							\
    }									\
  } while(0)

#define LOOP_NEGATIVE(commands)						\
  do {									\
    lastn = 0;								\
    int stattimes = 0;							\
    for(;;) {								\
      fracs[COEFFICIENT_BLOCK_SIZE-1] =					\
	cminusfrac(lastn-COEFFICIENT_BLOCK_SIZE+1+nu,			\
		   epsilon, q, m, lambda, NULL);			\
      for(n=COEFFICIENT_BLOCK_SIZE-2;n>=0;n--) {			\
	Complex LOOP_denom = beta(lastn-n-1) + fracs[n+1];		\
	if(norm(LOOP_denom) >						\
	   norm(fracs[n+1])*(LOOP_CANCEL_TOL*LOOP_CANCEL_TOL)) {	\
	  fracs[n] = - alphagamma(lastn-n-1) / LOOP_denom;		\
	} else {							\
	  fracs[n] = cminusfrac(lastn-n+nu, epsilon, q, m, lambda, NULL); \
	}								\
      }									\
      									\
      for(n=0;n<COEFFICIENT_BLOCK_SIZE;n++,lastn--) {			\
	lastcoeff *= fracs[n]/gamma(lastn);				\
	commands;							\
      }									\
      if(lastn < -LOOP_MAX)						\
	return NAN;							\
    }									\
  } while(0)

Complex FT_PRFX fsum(Complex nu, Real epsilon, Real q, int m, Real lambda) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q)/(1 - q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  int n;
  Complex sum;
  Complex lastcoeff;
  int lastn;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  /* I'm fairly sure generating a higher fraction value and then applying
   * the recursive definitions can't lose us any precision, and it'll be
   * faster than evaluating the entire thing for each f
   * We want to use the + fraction for positive and - for negative because
   * they don't cancel as much with the betas */
  sum = 1;
  lastcoeff = 1;

  LOOP_POSITIVE(if(sum == sum + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += lastcoeff;
		if(Cisnan(sum)) return CNAN;
		);
 plusfrac_done:

  lastcoeff = 1;

  LOOP_NEGATIVE(if(sum == sum + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += lastcoeff;
		if(Cisnan(sum)) return CNAN;
		);
 minusfrac_done:

  return sum;
}

Complex FT_PRFX kfactor(Complex nu, Real epsilon, Real q, int m, Real lambda) {
  /* in some applications, this will be called repeatedly on the same
   * values, so implement caching */
  static int cache_next = 0;
  static Complex cache[K_CACHE_SIZE] = K_CACHE_INIT;
  static Complex cache_nu[K_CACHE_SIZE] = K_CACHE_INIT;
  static Real cache_epsilon[K_CACHE_SIZE] = K_CACHE_INIT;
  static Real cache_q[K_CACHE_SIZE] = K_CACHE_INIT;
  static int cache_m[K_CACHE_SIZE] = K_CACHE_INIT;
  static Real cache_lambda[K_CACHE_SIZE] = K_CACHE_INIT;

  int i;
  for(i=0;i<K_CACHE_SIZE;i++) {
    if(cache_nu[i] == nu &&
       cache_epsilon[i] == epsilon &&
       cache_q[i] == q &&
       cache_m[i] == m &&
       cache_lambda[i] == lambda)
      return cache[i];
  }

  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex prefact, num, denom;
  Complex gammahold;

  Complex lastcoeff;
  int lastn, n;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  const Complex nu2_1 = 2.0*nu + 1.0;
  const Complex nu_3_ie = nu + 3.0 + I*epsilon;
  const Complex nu_3__ie = nu + 3.0 - I*epsilon;
  const Complex nu__1_ie = nu - 1.0 + I*epsilon;
  const Complex nu_1_it = nu + 1.0 + I*tau;
  const Complex nu_1__it = nu + 1.0 - I*tau;
  const Complex nu_2_ie = nu + 2.0 + I*epsilon;
  const Complex nu__2__ie = nu - 2.0 - I*epsilon;

  /* taking N=0 */

  prefact = gammln(3.0 - I*(epsilon+tau)) + gammln(2.0*nu + 2.0)
    - gammln(nu_3_ie) + I*epsilon*kappa + gammln(nu2_1)
    - gammln(nu_3__ie) - gammln(nu_1__it);
  prefact = 4.0 * pow(2*epsilon*kappa, -2.0-nu) * exp(prefact);

  lastcoeff = 1;
  num = 1;
  LOOP_POSITIVE(lastcoeff *=
		- ((lastn + nu2_1) * (lastn + nu__1_ie) * (lastn + nu_1_it))
		/ ((lastn + 1) * (lastn + nu_3__ie) * (lastn + nu_1__it));
		if(num == num + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		num += lastcoeff;
		if(Cisnan(num)) return CNAN;
		);
 plusfrac_done:

  lastcoeff = 1;
  denom = 1;
  LOOP_NEGATIVE(lastcoeff *=
		((lastn + nu2_1) * (lastn + nu_2_ie))
		/ ((lastn - 1) * (lastn + nu__2__ie));
		if(denom == denom + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		denom += lastcoeff;
		if(Cisnan(denom)) return CNAN;
		);
 minusfrac_done:

  Complex ans;
  cache_nu[cache_next] = nu;
  cache_epsilon[cache_next] = epsilon;
  cache_q[cache_next] = q;
  cache_m[cache_next] = m;
  cache_lambda[cache_next] = lambda;
  cache[cache_next] = ans = prefact * num / denom;
  if(++cache_next >= K_CACHE_SIZE) cache_next = 0;

  return ans;
}

Complex FT_PRFX aminus(Complex nu, Real epsilon, Real q, int m, Real lambda) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex sum;

  int n;
  Complex lastcoeff;
  int lastn;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  const Complex nu__1__ie = nu - 1.0 - I*epsilon;
  const Complex nu_3_ie = nu + 3.0 + I*epsilon;
  const Complex nu_2_ie = nu + 2.0 + I*epsilon;
  const Complex nu__2__ie = nu - 2.0 - I*epsilon;

  sum = 1;
  lastcoeff = 1;

  LOOP_POSITIVE(lastcoeff *= - (lastn + nu__1__ie) / (lastn + nu_3_ie);
		if(sum == sum + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += lastcoeff;
		if(Cisnan(sum)) return CNAN;
		);
 plusfrac_done:

  lastcoeff = 1;

  LOOP_NEGATIVE(lastcoeff *= - (lastn + nu_2_ie) / (lastn + nu__2__ie);
		if(sum == sum + lastcoeff) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += lastcoeff;
		if(Cisnan(sum)) return CNAN;
		);
 minusfrac_done:

  return pow(2.0, 1.0 + I*epsilon)
    * exp(-0.5 * M_PI * (epsilon + I*(nu - 1.0))) * sum;
}

Complex FT_PRFX rzero(Complex nu, Real epsilon, Real q, int m, Real lambda, Real x,
		      Complex *deriv) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex prefact,sum,derivsum;

  LongComplex lastcoeff; // to avoid overflow
  int lastn, n;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  const Complex nu2_1 = nu*2.0 + 1.0;
  const Complex nu2__1 = nu*2.0 - 1.0;
  const Complex nu2__3 = nu*2.0 - 3.0;
  const Complex nu_it = nu + I*tau;
  const Complex nu__it = nu - I*tau;
  const Complex nu_1 = nu + 1.0;
  const Complex nu_1_it = nu_1 + I*tau;
  const Complex nu_1__it = nu_1 - I*tau;
  const Complex nu__1 = nu - 1.0;
  const Complex nu__1_ie = nu__1 + I*epsilon;
  const Complex nu__1_it = nu__1 + I*tau;
  const Complex nu_2__ie = nu + 2.0 - I*epsilon;
  const Complex nu__2_ie = nu - 2.0 + I*epsilon;
  const Complex nu__2_it = nu - 2.0 + I*tau;
  const Complex nu_3__ie = nu + 3.0 - I*epsilon;
  const Complex nu__3_ie = nu - 3.0 + I*epsilon;
  const Complex nu__4_ie = nu - 4.0 + I*epsilon;

  const Real xflip = 1 - x;
  const Real xflipinv = 1 / xflip;

  prefact = exp(I*epsilon*kappa*x
		+ gammln(3.0 - I*(epsilon+tau)) + gammln(nu2_1)
		- gammln(nu_1__it) - gammln(nu_3__ie))
    * pow(-x, 2.0 - 0.5*I*(epsilon+tau))
    * pow(xflip, nu + 0.5*I*(epsilon+tau));

  lastcoeff = 1;
  sum = hypergeom2F1(-nu_it, -nu__2_ie, -2.0*nu, xflipinv);
  derivsum = 0;
  if(deriv)
    derivsum = xflipinv * (nu_it)*(nu__2_ie)/(2.0*nu)
      * hypergeom2F1(-nu__1_it, -nu__3_ie, -nu2__1, xflipinv);

  LOOP_POSITIVE(lastcoeff *= (2.0*lastn + nu2_1) * (2.0*(lastn + nu_1)) * xflip
		/ ((lastn + nu_1__it) * (lastn + nu_3__ie));

		Complex term;
		term = lastcoeff
		* hypergeom2F1(-lastn - nu_1_it, -lastn - nu__1_ie,
			       -2.0*(lastn + nu_1), xflipinv);
		Complex derivterm = 0;
		if(deriv)
		  derivterm =
		    (lastn + 1.0) * term
		     + lastcoeff
		     * xflipinv * (lastn + nu_1_it) * (lastn + nu__1_ie)
		     / (2.0*(lastn + nu_1))
		     * hypergeom2F1(-lastn - nu_it,-lastn - nu__2_ie,
				    -2.0*lastn - nu2_1, xflipinv);
		if(sum == sum + term &&
		   derivsum == derivsum + derivterm) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 plusfrac_done:

  lastcoeff = 1;
  LOOP_NEGATIVE(lastcoeff *= (lastn + nu__it) * (lastn + nu_2__ie)
		/ ((2.0*(lastn + nu)) * (2.0*lastn + nu2__1) * xflip);

		Complex term;
		term = lastcoeff
		* hypergeom2F1(-lastn - nu__1_it, -lastn - nu__3_ie,
			       -2.0*(lastn + nu__1), xflipinv);
		Complex derivterm = 0;
		if(deriv)
		  derivterm = (lastn - 1.0) * term
		    + lastcoeff * xflipinv * (lastn + nu__1_it)
		    * (lastn + nu__3_ie) / (2.0*(lastn + nu__1))
		    * hypergeom2F1(-lastn - nu__2_it, -lastn - nu__4_ie,
				   -2.0*lastn - nu2__3, xflipinv);
		if(sum == sum + term &&
		   derivsum == derivsum + derivterm) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 minusfrac_done:

  if(deriv)
    *deriv = 0.5 * prefact / kappa
      * (xflipinv * derivsum
	 - (I*epsilon*kappa - (I*0.5*(epsilon+tau) - 2.0)/x
	    - (I*0.5*(epsilon+tau) + nu)*xflipinv) * sum);
  return prefact * sum;
}

Complex FT_PRFX rin_hyper(Complex nu, Real epsilon, Real q, int m, Real lambda,
			  Real x, Real /*precision*/, Complex *deriv) {
  if(!deriv)
    return
      rzero(nu, epsilon, q, m, lambda, x, NULL) +
      rzero(-nu - 1.0, epsilon, q, m, lambda, x, NULL);

  Complex val,deriv1,deriv2;
  val = rzero(nu, epsilon, q, m, lambda, x, &deriv1) +
    rzero(-nu - 1.0, epsilon, q, m, lambda, x, &deriv2);
  *deriv = deriv1 + deriv2;
  return val;
}

Complex FT_PRFX rcoulomb(Complex nu, Real epsilon, Real q, int m, Real lambda,
			 Real z, Complex *deriv) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex prefact,sum,derivsum;

  LongComplex lastcoeff; // to avoid overflow
  int lastn, n;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  Complex iz2 = I*z*2.0;

  const Complex nu2_1 = nu * 2.0 + 1.0;
  const Complex nu2_3 = nu * 2.0 + 3.0;
  const Complex nu2_5 = nu * 2.0 + 5.0;
  const Complex nu_1 = nu + 1.0;
  const Complex nu__1__ie = nu - 1.0 - I*epsilon;
  const Complex nu_2 = nu + 2.0;
  const Complex nu_2_ie = nu_2 + I*epsilon;
  const Complex nu__2__ie = nu - 2.0 - I*epsilon;
  const Complex nu_3_ie = nu + 3.0 + I*epsilon;
  const Complex nu_4_ie = nu + 4.0 + I*epsilon;
  const Complex nu_5_ie = nu + 5.0 + I*epsilon;

  prefact = pow(2.0,nu) * pow(z,nu_2)
    * pow(1-epsilon*kappa/z,2.0-I*0.5*(epsilon+tau))
    * exp(gammln(nu_3_ie) - gammln(2.0*nu_1) - I*z);

  lastcoeff = 1;
  sum = hypergeom1F1(nu_3_ie, 2.0*nu_1, iz2);
  derivsum = 0;
  if(deriv)
    derivsum = I * z * nu_3_ie / nu_1 * hypergeom1F1(nu_4_ie,nu2_3,iz2);

  LOOP_POSITIVE(lastcoeff *= -iz2 * (lastn + nu__1__ie)
		/ ((2.0*lastn + nu2_3) * (2.0*(lastn + nu_1)));

		Complex term;
		term = lastcoeff
		* hypergeom1F1(lastn + nu_4_ie, 2.0*(lastn + nu_2), iz2);
		Complex derivterm = 0;
		if(deriv)
		  derivterm = (lastn+1) * term
		    + lastcoeff * I * z * (lastn + nu_4_ie) / (lastn+nu_2)
		    * hypergeom1F1(lastn + nu_5_ie, 2*lastn + nu2_5, iz2);
		if(sum == sum + term &&
		   derivsum == derivsum + derivterm) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 plusfrac_done:

  lastcoeff = 1;
  LOOP_NEGATIVE(lastcoeff *= (2.0*lastn+nu2_1) * (2.0*(lastn + nu))
		/ (-iz2 * (lastn + nu__2__ie));

		Complex term;
		term = lastcoeff
		* hypergeom1F1(lastn + nu_2_ie, 2.0*(lastn + nu), iz2);
		Complex derivterm = 0;
		if(deriv)
		  derivterm = (lastn-1) * term
		    + lastcoeff * I * z * (lastn + nu_2_ie) / (lastn + nu)
		    * hypergeom1F1(lastn + nu_3_ie, 2*lastn + nu2_1, iz2);
		if(sum == sum + term &&
		   derivsum == derivsum + derivterm) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 minusfrac_done:

  if(deriv)
    *deriv = 0.5 * epsilon * prefact
      * (derivsum / z +
	 (-I + (nu+2) / z +
	  (2 - 0.5*I*(epsilon+tau)) / (z * (z / (epsilon*kappa) - 1))) * sum);
  return prefact * sum;
}

Complex FT_PRFX rin_coulomb(Complex nu, Real epsilon, Real q, int m,
			    Real lambda, Real z, Real /*precision*/,
			    Complex *deriv) {
  Complex k_nu = kfactor(nu, epsilon, q, m, lambda);
  Complex k__nu__1 = kfactor(-nu - 1.0, epsilon, q, m, lambda);
  if(!deriv)
    return
      k_nu * rcoulomb(nu, epsilon, q, m, lambda, z, NULL) +
      k__nu__1 * rcoulomb(-nu - 1.0, epsilon, q, m, lambda, z, NULL);

  Complex val, deriv1, deriv2;
  val = k_nu * rcoulomb(nu, epsilon, q, m, lambda, z, &deriv1) +
    k__nu__1 * rcoulomb(-nu - 1.0, epsilon, q, m, lambda, z, &deriv2);
  *deriv = k_nu * deriv1 + k__nu__1 * deriv2;
  return val;
}

/*
Complex rin(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
	    Complex *deriv) {
  Real x = (1+(1-r)/Sqrt(1-q*q))/2;

  / * This seems to be a reasonable cutoff * /
  return x > -10 ?
    rin_small(nu, epsilon, q, m, lambda, x, REAL_EPSILON, deriv) :
    rin_coulomb(nu, epsilon, q, m, lambda,
		0.5*epsilon*(r - 1 + Sqrt(1-q*q)), deriv);
}
*/

Complex FT_PRFX rup_tricomi(Complex nu, Real epsilon, Real q, int m,
			    Real lambda, Real z, Real precision,
			    Complex *deriv) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex prefact,sum,derivsum;

  LongComplex lastcoeff; // to avoid overflow
  int lastn, n;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  Complex iz2 = I*z*2.0;

  const Complex nu2_1 = nu*2.0 + 1.0;
  const Complex nu2_2 = nu*2.0 + 2.0;
  const Complex nu2_5 = nu*2.0 + 5.0;
  const Complex nu__ie = nu - I*epsilon;
  const Complex nu_1 = nu + 1.0;
  const Complex nu_1__ie = nu__ie + 1.0;
  const Complex nu__1__ie = nu__ie - 1.0;
  const Complex nu_2 = nu + 2.0;
  const Complex nu_2_ie = nu_2 + I*epsilon;
  const Complex nu__2__ie = nu__ie - 2.0;
  const Complex nu_3_ie = nu + 3.0 + I*epsilon;

  prefact = pow(2.0,nu) * pow(z,nu_2)
    * pow(1-epsilon*kappa/z, 2.0 - I*0.5*(epsilon+tau))
    * exp(I*z - M_PI*epsilon - I*M_PI*(nu-1.0));

  lastcoeff = pow(-2*I*z,-2*nu-1) * exp(gammln(nu2_1) - gammln(nu__1__ie));
  sum = hypergeomU(nu-1.0-I*epsilon, 2.0*(nu+1.0), -iz2);
  derivsum = 0;
  if(deriv)
    derivsum = iz2 * (nu - 1.0 - I*epsilon)
      * hypergeomU(nu - I*epsilon,2.0*nu+3.0,-iz2);

  LOOP_POSITIVE(lastcoeff *= (2*lastn + nu2_1) * 2*(lastn + nu_1)
		/ ((lastn + nu_3_ie) * iz2);

		Complex term;
		term = lastcoeff
		* hypergeomU_reduced(-lastn - nu_3_ie, -2.0*(lastn + nu_1),
				     -iz2);
		Complex derivterm = 0;
		if(deriv)
		  derivterm = -(lastn + nu2_2) * term
		    - lastcoeff * iz2 * (lastn + nu_3_ie)
		    / (2.0*(lastn + nu_1))
		    * hypergeomU_reduced(-lastn - nu_2_ie, -2.0*lastn - nu2_1,
					 -iz2);
		if(abs(sum) * precision >= abs(term) &&
		   abs(derivsum) * precision >= abs(derivterm)) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 plusfrac_done:

  lastcoeff = exp(gammln(-nu2_1) - gammln(-nu_2_ie));
  LOOP_NEGATIVE(lastcoeff *= -(2*lastn + nu2_1) * 2*(lastn + nu)
		/ (iz2 * (lastn + nu__2__ie));

		Complex term;
		term = lastcoeff
		* hypergeomU_reduced(lastn + nu__2__ie, 2.0*(lastn + nu),
				     -iz2);
		Complex derivterm = 0;
		if(deriv)
		  derivterm = (lastn - 1.0) * term
		    - lastcoeff * iz2 * (lastn + nu__2__ie)
		    / (2.0*(lastn + nu))
		    * hypergeomU_reduced(lastn + nu__1__ie, 2.0*lastn + nu2_1,
					 -iz2);
		if(abs(sum) * precision >= abs(term) &&
		   abs(derivsum) * precision >= abs(derivterm)) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 minusfrac_done:

  if(deriv)
    *deriv = 0.5 * epsilon * prefact
      * (derivsum / z +
	 (I + (nu + I*0.5*(epsilon+tau)
	       - (I*0.5*(epsilon+tau) - 2.0)/(1 - epsilon*kappa/z))/z) * sum);
  return prefact * sum;
}

/*
Complex rup(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
	    Complex *deriv) {
  if(r > 3) {
    Real z = epsilon * (r - 1 + Sqrt(1-q*q))/2;
    return rup_tricomi(nu,epsilon,q,m,lambda,z,deriv);
  } else {
    Real x = ((1-r)/Sqrt(1-q*q) + 1) / 2;
    return rup_hyper(nu,epsilon,q,m,lambda,x,deriv);
  }
}
*/

Complex FT_PRFX rin_small(Complex nu, Real epsilon, Real q, int m, Real lambda,
			  Real x, Real precision, Complex *deriv) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q) / (1-q*q);

  PREP_ALPHAGAMMA(Complex);
  PREP_BETA(Complex);
  PREP_ALPHA_GAMMA(); /* always complex */

  Complex prefact,sum,derivsum;

  LongComplex lastcoeff; // to avoid overflow
  int lastn, n;
  Complex fracs[COEFFICIENT_BLOCK_SIZE];

  const Complex nu__it = nu - I*tau;
  const Complex nu_1_it = nu + 1.0 + I*tau;
  const Complex nu_1__it = nu + 1.0 - I*tau;
  const Complex nu__1_it = nu - 1.0 + I*tau;
  const Complex nu_2__ie = nu + 2.0 - I*epsilon;
  const Complex nu_2__it = nu + 2.0 - I*tau;
  const Complex nu_3__ie = nu + 3.0 - I*epsilon;
  const Complex nu_3__it = nu + 3.0 - I*tau;
  const Complex nu_4__ie = nu + 4.0 - I*epsilon;
  const Complex nu_5__ie = nu + 5.0 - I*epsilon;
  const Complex nu0_3__ie__it = 3.0 - I*(epsilon+tau);
  const Complex nu0_4__ie__it = 4.0 - I*(epsilon+tau);

  const Real xflip = 1 - x;
  const Real xreduced = - x / xflip;

  prefact = exp(I*epsilon*kappa*x)
    * pow(-x, 2.0 - 0.5*I*(epsilon+tau))
    * pow(xflip, -nu - 1 + 0.5*I*(epsilon+tau));

  lastcoeff = 1;
  sum = hypergeom2F1(nu_1__it, nu_3__ie, nu0_3__ie__it, xreduced);
  derivsum = 0;
  if(deriv)
    derivsum = nu_1__it * (nu + I*tau)
      * hypergeom2F1(nu_2__it, nu_3__ie, nu0_4__ie__it, xreduced);

  LOOP_POSITIVE(lastcoeff /= xflip;

		Complex term;
		term = lastcoeff
		* hypergeom2F1(lastn + nu_2__it, lastn + nu_4__ie,
			       nu0_3__ie__it, xreduced);
		Complex derivterm = 0;
		if(deriv)
		  derivterm =
		    lastcoeff * (lastn + nu_2__it) * (lastn + nu_1_it)
		    * hypergeom2F1(lastn + nu_3__it, lastn + nu_4__ie,
				   nu0_4__ie__it, xreduced);
		if(abs(sum) * precision >= abs(term) &&
		   abs(derivsum) * precision >= abs(derivterm)) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto plusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 plusfrac_done:

  lastcoeff = 1;
  LOOP_NEGATIVE(lastcoeff *= xflip;

		Complex term;
		term = lastcoeff
		* hypergeom2F1(lastn + nu__it, lastn + nu_2__ie,
			       nu0_3__ie__it, xreduced);
		Complex derivterm = 0;
		if(deriv)
		  derivterm =
		    lastcoeff * (lastn + nu__it) * (lastn + nu__1_it)
		    * hypergeom2F1(lastn + nu_1__it, lastn + nu_2__ie,
				   nu0_4__ie__it, xreduced);
		if(abs(sum) * precision >= abs(term) &&
		   abs(derivsum) * precision >= abs(derivterm)) {
		  if(++stattimes >= LOOP_STAT_COUNT) goto minusfrac_done;
		} else { stattimes = 0; }
		sum += term;
		derivsum += derivterm;
		if(Cisnan(sum) || Cisnan(derivsum)) {
		  if(deriv) *deriv = CNAN;
		  return CNAN;
		}
		);
 minusfrac_done:

  if(deriv)
    *deriv = 0.5 * prefact / kappa
      * (derivsum / (xflip * nu0_3__ie__it)
	 - (I*epsilon*kappa - (I*0.5*(epsilon+tau) - 2.0)/x
	    - I*0.5*(epsilon-tau)/xflip) * sum);
  return prefact * sum;
}

Complex FT_PRFX rup_hyper(Complex nu, Real epsilon, Real q, int m, Real lambda,
			  Real x, Real /*precision*/, Complex *deriv) {
  Complex k_nu = kfactor(nu, epsilon, q, m, lambda);
  Complex k__nu__1 = kfactor(-nu-1.0, epsilon, q, m, lambda);

  Complex prefact = exp(-M_PI * epsilon) / sin(2*M_PI*nu);
  Complex coeff_nu = exp(-I*M_PI*nu) * sin(M_PI*(nu+2+I*epsilon)) / k_nu;
  Complex coeff__nu__1 =
    -I * sin(M_PI*(nu-2-I*epsilon)) / k__nu__1;

  Complex deriv1, deriv2;

  Complex val1 = coeff_nu * rzero(nu, epsilon, q, m, lambda, x, &deriv1);
  Complex val2 = coeff__nu__1
    * rzero(-nu - 1.0, epsilon, q, m, lambda, x, &deriv2);
  *deriv = prefact * (coeff_nu * deriv1 + coeff__nu__1 * deriv2);

  return prefact * (val1 + val2);
}
