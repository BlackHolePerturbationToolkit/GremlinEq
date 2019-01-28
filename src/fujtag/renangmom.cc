//---------------------------------------------------------------------------
//
// $Id: renangmom.cc,v 1.3 2014/04/08 19:36:44 sahughes Exp $
//
//---------------------------------------------------------------------------

/* renangmom.cc
 * calculates the renormalized angular momentum nu
 */

#include <iostream>

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "renangmom.h"
#  include "radialfrac.h"
#  include "specialradial.h"
#endif
#include <cstdio>
#include <cmath>

#define RENANGMOM_PRECISION 1e-15

/* number of points on either side of the main value for get_deciders()
 * noise calculations.  Set to a false value (or undef) to use only term
 * estimation. */
#define DECIDERS_NOISE_POINTS 2

/* similar for renangmom_half(), but no fallback */
#define HALF_NOISE_POINTS 2

/* distance a root must be from a half integer to be loacated */
#define RENANGMOM_DETECT_LIMIT 1e-12

/* first separation from a half integer to try for the real searcher */
#define RENANGMOM_REAL_DETECT_ATTEMPT 1e-3

/* the first point to try to try to account for numerical problems for high
 * epsilon (where Im(nu) is probably large anyway)
 * RENENGMOM_DETECT_LIMIT will be tried if this fails */
#define RENANGMOM_IINT_DETECT_ATTEMPT .1

/* The upper bound (as a multiple of epsilon) for the _iint solver */
#define RENANGMOM_MAX_DIST 2

/* The maximum number of divisions the secant code will divide the interval
 * into performing its search.  Should be a power of 2 */
#define RENANGMOM_IINT_MAX_DIVISIONS 1024

/* maximum number of steps the secant code will execute before giving up */
#define RENANGMOM_SECANT_MAX_STEPS 100

/* relative step size that will trigger success */
#define RENANGMOM_SECANT_MIN_REL_STEP 1e-15

/* largest allowed ratio of imag(step)/real(step) */
#define SECANT_IMAG_LIMIT 1e4

/* difference between two secant solutions that will allow passing through
 * the skeptic check */
#define RENANGMOM_SKEPTIC_REPEAT_SPACING 1e-9

/* initial step size for the last ditch step search in the real solver */
#define DIVIDE_SEARCH_STEP 0.01

/* scaling down ratio for the divide searcher */
#define DIVIDE_SEARCH_STEP_SCALE 0.25

/* minimum search scale for the divide search solver */
#define MIN_DIVIDE_SEARCH_STEP 1e-4

/* the fraction of the predictor function below which we become skeptical
 * of the results of the secant search */
#define SECANT_SKEPTIC_MULT 0.1

/* the assumed width of the positive region for the half case.
 * Only used if the first attmepts give garbage. */
#define HALF_POS_SEARCH_LIMIT 2

/* the minimum distance from noise for which a point will be accepted
 * this check is disabled for any proposed solution less than twice this
 * distance form the real axis */
#define IINT_NOISE_SEP_MIN 1

/* the spacing of the points used in the singularity fitting routine */
#define SING_FIT_SPACINGS { 1e-3, 5e-3, 1e-2 }

/* the maximum difference between the two singularity estimates which
 * will allow the value to be acccepted */
#define SING_FIT_THRESHOLD 1e-5

/* calculates nu
 * returns 1 if it suceeds, 0 else */
int FT_PRFX renangmom(Real epsilon, Real q, int l, int m, Real lambda,
		      Complex *nu) {
  Real nupart;

  if(epsilon == 0.0) {
    *nu = 1;
    return 1;
  }

  /* Should be able to locate roots by:
   * l-1+eps l-0.5-eps l-0.5+eps l-eps
   * -       +         -         +     real
   * -       -         +         +     half also frac_half(eps) > 0
   * +       +         -         -     iint also Im(frac_iint(eps)) > 0
   *                                             or d/d(nu) frac(1) > 0
   * Last line is incorrect.  I think the correct test is that 
   * d/d(nu) frac(1) has the same sign as epsilon - m*q
   */

  Real halfval, oneslope;
  if(!get_deciders(epsilon,q,m,lambda,&halfval,&oneslope))
    return renangmom_far(epsilon,q,l,m,lambda,nu);

  /* this checks if an odd number of the inequalities hold */
  int iinttest = ((oneslope > 0) == (epsilon >= m*q)) == (epsilon > 0);

  if(halfval < 0) {
    if(iinttest) {
      if(renangmom_iint(epsilon, q, l, m, lambda, &nupart)) {
	*nu = 1.0 + I*nupart;
	return 1;
      }
    } else {
      if(renangmom_real_search_procedure(epsilon, q, l, m, lambda, nu)) {
	return 1;
      }
    }
  } else {
    if(iinttest) {
      /* Root placement tests inconsistent */
    } else {
      if(renangmom_half(epsilon, q, l, m, lambda, &nupart)) {
	*nu = - 0.5 + I*nupart;
	return 1;
      }
    }
  }

  /* Search failed - might as well throw everything at it */
  return renangmom_far(epsilon,q,l,m,lambda,nu);
}

int FT_PRFX get_deciders(Real epsilon, Real q, int m, Real lambda,
			 Real *halfval, Real *oneslope) {
  Real term;
  *halfval = radial_half_value(epsilon,q,m,lambda,&term);
  if(term * SPECIAL_TERM_NOISE > Fabs(*halfval))
    return 0;
#if DECIDERS_NOISE_POINTS
  int i;
  Real testval, maxval, minval;
  maxval = minval = *halfval;
  for(i=0;i<DECIDERS_NOISE_POINTS;i++) {
    testval = radial_half_value(epsilon*(1+i*REAL_EPSILON),q,m,lambda,NULL);
    if(testval < minval) {
      minval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    } else if(testval > maxval) {
      maxval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    }
    testval = radial_half_value(epsilon*(1-i*REAL_EPSILON),q,m,lambda,NULL);
    if(testval < minval) {
      minval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    } else if(testval > maxval) {
      maxval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    }
  }
#endif /* DECIDERS_NOISE_POINTS */
  *oneslope = radial_int_slope(epsilon,q,m,lambda,&term);
  if(term * SPECIAL_TERM_NOISE > Fabs(*oneslope))
    return 0;
#if DECIDERS_NOISE_POINTS
  minval = maxval = *oneslope;
  for(i=0;i<DECIDERS_NOISE_POINTS;i++) {
    testval = radial_int_slope(epsilon*(1+i*REAL_EPSILON),q,m,lambda,NULL);
    if(testval < minval) {
      minval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    } else if(testval > maxval) {
      maxval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    }
    testval = radial_int_slope(epsilon*(1-i*REAL_EPSILON),q,m,lambda,NULL);
    if(testval < minval) {
      minval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    } else if(testval > maxval) {
      maxval = testval;
      if(maxval-minval > Fmin(Fabs(minval),Fabs(maxval)))
	return 0;
    }
  }
#endif /* DECIDERS_NOISE_POINTS */
  return 1;
}

int FT_PRFX renangmom_real_search_procedure(Real epsilon, Real q, int l,
					    int m, Real lambda, Complex *nu) {
  Real nupart;
  Real guess = fractions_poly_approx(epsilon, q, m, lambda);

  if(!isnan(guess)) {
    if(renangmom_real_guess(epsilon, q, m, lambda, guess, &nupart)) {
      *nu = nupart;
      return 1;
    }
  } else {
    guess = l - 0.75;
  }

  int guessint = (int) ceil(guess);
  int guessint2 = guessint + (guess>ceil(guess)-.5 ? 1 : -1);

  if(renangmom_real(epsilon, q, guessint, m, lambda,
		    RENANGMOM_REAL_DETECT_ATTEMPT, &nupart) ||
     renangmom_real(epsilon, q, guessint2, m, lambda,
		    RENANGMOM_REAL_DETECT_ATTEMPT, &nupart) ||
     renangmom_real_divide_search(epsilon, q, guessint, m, lambda,
				  &nupart) ||
     renangmom_real(epsilon, q, guessint, m, lambda,
		    RENANGMOM_DETECT_LIMIT, &nupart) ||
     renangmom_real(epsilon, q, guessint2, m, lambda,
		    RENANGMOM_DETECT_LIMIT, &nupart)) {
    *nu = nupart;
    return 1;
  }

  /* This often fails near a singularity.  Try fitting to the sides */
  const static Real spacings[] = SING_FIT_SPACINGS;
  const static int numspacings = sizeof(spacings)/sizeof(Real);
  for(int i=0;i<numspacings;i++) {
    if(singularity_fit(epsilon, q, l, m, lambda, guessint, guessint2,
		       spacings[i], nu))
      return 1;
  }
  return 0;
}

int FT_PRFX singularity_fit(Real epsilon, Real q, int l, int m, Real lambda,
			    int guessint, int guessint2, Real spacing,
			    Complex *nu) {
  Real leftval[3];
  Real rightval[3];
  int step;

  for(step=0;step<=2;step++) {
    if(!(renangmom_real(epsilon-(step+1)*spacing, q, guessint, m, lambda,
			RENANGMOM_REAL_DETECT_ATTEMPT, leftval+step) ||
	 renangmom_real(epsilon-(step+1)*spacing, q, guessint2, m, lambda,
			RENANGMOM_REAL_DETECT_ATTEMPT, leftval+step) ||
	 renangmom_real_divide_search(epsilon-(step+1)*spacing, q, guessint,
				      m, lambda, leftval+step) ||
	 renangmom_real(epsilon-(step+1)*spacing, q, guessint, m, lambda,
			RENANGMOM_DETECT_LIMIT, leftval+step) ||
	 renangmom_real(epsilon-(step+1)*spacing, q, guessint2, m, lambda,
			RENANGMOM_DETECT_LIMIT, leftval+step)) ||
       !(renangmom_real(epsilon+(step+1)*spacing, q, guessint, m, lambda,
			RENANGMOM_REAL_DETECT_ATTEMPT, rightval+step) ||
	 renangmom_real(epsilon+(step+1)*spacing, q, guessint2, m, lambda,
			RENANGMOM_REAL_DETECT_ATTEMPT, rightval+step) ||
	 renangmom_real_divide_search(epsilon+(step+1)*spacing, q, guessint,
				      m, lambda, rightval+step) ||
	 renangmom_real(epsilon+(step+1)*spacing, q, guessint, m, lambda,
			RENANGMOM_DETECT_LIMIT, rightval+step) ||
	 renangmom_real(epsilon+(step+1)*spacing, q, guessint2, m, lambda,
			RENANGMOM_DETECT_LIMIT, rightval+step)))
      return 0;
    leftval[step]  += nearbyint(leftval[0] - leftval[step]);
    rightval[step] += nearbyint(leftval[0] - rightval[step]);
  }

  Real leftest = leftval[2] + 3*(leftval[0]-leftval[1]);
  Real rightest = rightval[2] + 3*(rightval[0]-rightval[1]);
  if(Fabs(leftest-rightest) < SING_FIT_THRESHOLD) {
    *nu = 0.5*(leftest+rightest);
    fprintf(stderr,
	    "WARNING: Used singularity fiting at l=%d m=%d epsilon=%.17f to find nu=%f (delta %e)\n",
	    l, m, epsilon, real(*nu), Fabs(leftest-rightest));
    return 1;
  }
  return 0;
}

/* Brent's method
 * argument should be f(c)
 * before calling, make sure (a,c) brackets the root, and aval, cval are set
 *
 * note that we know we won't converge to zero, so can stick to
 * relative precision */
#define BRENT(func)							\
  ({									\
    Real b, bval;							\
    Real step, oldstep;							\
    Real num, denom, bcdiff, acvdiff, bcvdiff, tol;			\
    b = a;								\
    bval = aval;							\
    oldstep = step = c - a;						\
									\
    for(;;) {								\
      if((bval > 0) == (cval > 0)) {					\
	b = a;								\
	bval = aval;							\
	oldstep = step = Fabs(c - a);					\
      }									\
      if(Fabs(bval) < Fabs(cval)) {					\
	a = c;								\
	aval = cval;							\
	c = b;								\
	cval = bval;							\
	b = a;								\
	bval = aval;							\
      }									\
									\
      bcdiff = b - c;							\
      tol = Fabs(c)*RENANGMOM_PRECISION;				\
      if(cval == 0 || Fabs(bcdiff) <= tol) break;			\
									\
      if(oldstep > tol && Fabs(cval) < Fabs(aval)) {			\
	acvdiff = aval - cval;						\
	if(a == b) {							\
	  /* linear interpolation */					\
	  num = bcdiff * cval;						\
	  denom = acvdiff;						\
	} else {							\
	  /* inverse quadratic interpolation */				\
	  bcvdiff = bval - cval;					\
	  num = (bcdiff * acvdiff * aval + (c - a) * bcvdiff * bval) * cval; \
	  denom = acvdiff * bcvdiff * (aval - bval);			\
	}								\
	if(num > 0) { denom *= -1; }					\
	else { num *= -1; }						\
									\
	if(2*num < 1.5*bcdiff*denom - tol*Fabs(denom) &&		\
	   2*num < oldstep * Fabs(denom)) {				\
	  /* interpolation worked */					\
	  oldstep = Fabs(step);						\
	  step = num / denom;						\
	} else {							\
	  step = bcdiff / 2;						\
	  oldstep = Fabs(step);						\
	}								\
      } else {								\
	step = bcdiff / 2;						\
	oldstep = Fabs(step);						\
      }									\
      a = c;								\
      aval = cval;							\
									\
      if(Fabs(step) > tol) {						\
	c += step;							\
      } else if(c < b) {						\
	c += tol/2;							\
      } else {								\
	c -= tol/2;							\
      }									\
      cval = (func);							\
    }									\
  })

/* calulates nu assuming it is real.  returns 1 if it suceeds, 0 else
 * Uses Brent's algorithm */
int FT_PRFX renangmom_real(Real epsilon, Real q, int l, int m, Real lambda,
		       Real detect_limit, Real *nu) {
  Real a, c;
  Real aval, cval;

  Real sizeupper,sizelower;
  int upper, interval_forced;
  Real search_scale;

  /* First we need a bracketing interval
   * we assume here that the continued fraction function is generally
   * increasing (i.e. the slope is positive at the roots) on (l-1, l), and
   * that we can evaluate it without difficulty near the endpoints (which
   * can be wrong, particularly for large epsilon).  An alternate aproach is
   * attempted if that one fails
   * To choose which of (l-1,l-.5) or (l-.5,l) to search on, evaluate at the
   * midpoint of each and choose the one with the smaller value (this
   * hopefully picks out the better behaved half) */

  interval_forced = 0;
  sizeupper = Fabs(radialfrac(l-0.25,epsilon,q,m,lambda,NULL));
  sizelower = Fabs(radialfrac(l-0.75,epsilon,q,m,lambda,NULL));
  upper = (sizeupper < sizelower);

  /* Particularly for small values of epsilon, searching close to the roots
   * is important */
  search_scale = Fmin(epsilon*epsilon,1);

 choose_endpoints:
  if(upper) {
    a = (l - 0.5) + search_scale * detect_limit;
    c =  l        - search_scale * detect_limit;
  } else {
    a = (l - 1)   + search_scale * detect_limit;
    c = (l - 0.5) - search_scale * detect_limit;
  }

  aval = radialfrac(a, epsilon, q, m, lambda, NULL);
  if(aval > 0) {
    /* whoops.  Maybe we made a bad choice */
    if(interval_forced++) return 0;
    upper = !upper;
    goto choose_endpoints;
  }

  cval = radialfrac(c, epsilon, q, m, lambda, NULL);
  if(cval < 0) {
    /* whoops.  Maybe we made a bad choice */
    if(interval_forced++) return 0;
    upper = !upper;
    goto choose_endpoints;
  }

  BRENT(radialfrac(c, epsilon, q, m, lambda, NULL));

  /* sanity check: does this look like a root? */
  if(Fabs(cval) > Fabs(upper ? sizeupper : sizelower)) {
    /* no.  Maybe the other interval will work better */
    if(interval_forced++) return 0;
    upper = !upper;
    goto choose_endpoints;
  }

  *nu = c;
  return 1;
}

int FT_PRFX renangmom_real_divide_search(Real epsilon, Real q, int l, int m,
				     Real lambda, Real *nu) {
  Real lb, ub, lb_val, ub_val;
  Real current_search_step;
  for(current_search_step=DIVIDE_SEARCH_STEP;
      current_search_step>=MIN_DIVIDE_SEARCH_STEP;
      current_search_step*=DIVIDE_SEARCH_STEP_SCALE) {
    lb = l - 2;
    lb_val = 1;
    ub = l - 1 + current_search_step/2;
    ub_val = radialfrac(ub, epsilon, q, m, lambda, NULL);
    while(1) {
      lb = ub;
      ub += current_search_step;
      if(ub > l) break;
      lb_val = ub_val;
      ub_val = radialfrac(ub, epsilon, q, m, lambda, NULL);

      /* don't look at half integers */
      if((int) (2*lb) != (int) (2*ub)) continue;

      /* check if we're bracketing anything */
      if((lb_val > 0) || (ub_val < 0)) continue;

      Real a, c, aval, cval;
      a = lb;
      c = ub;
      aval = lb_val;
      cval = ub_val;
      BRENT(radialfrac(c, epsilon, q, m, lambda, NULL));

      /* sanity check: does this look like a root? */
      if(Fabs(cval) < ub_val) {
	*nu = c;
	return 1;
      }
    }
  }
  return 0;
}

/* calulates nu assuming it is on the special line.
 * Returns 1 if it suceeds, 0 else
 * Uses Brent's algorithm */
int FT_PRFX renangmom_half(Real epsilon, Real q, int l, int m, Real lambda,
			   Real *imnu) {
  Real a, c;
  Real aval, cval;
  Real badbound, term;

  /* First we need a bracketing interval */
  /* When the function starts out negative there is no root, and when it
   * starts out positive there is one
   * Of course, this function should only be called if there is one */
  a = RENANGMOM_IINT_DETECT_ATTEMPT;
  aval = radialfrac_half(a, epsilon, q, m, lambda, &term);
  if(aval < 0 && term * TERM_NOISE < -aval) {
    /* maybe we started too far out trying to be safe */
    a = RENANGMOM_DETECT_LIMIT;
    aval = radialfrac_half(a, epsilon, q, m, lambda, &term);
    if(aval < 0)
      return 0;
  }
  badbound = 0;
  /* search outward until we clear the noise */
  Real noiseest = INFINITY;
  while(noiseest > Fabs(aval)) {
    while(term * TERM_NOISE > Fabs(aval)) {
      badbound = a;
      a *= 2;
      if(a > 2*Fabs(epsilon))
	return 0;
      aval = radialfrac_half(a, epsilon, q, m, lambda, &term);
    }
    Real testval, minval, maxval;
    minval = maxval = aval;
    for(int i=0;i<HALF_NOISE_POINTS;i++) {
      testval = radialfrac_half(a*(1+i*REAL_EPSILON),epsilon,q,m,lambda,NULL);
      if(testval < minval) {
	minval = testval;
      } else if(testval > maxval) {
	maxval = testval;
      }
      testval = radialfrac_half(a*(1-i*REAL_EPSILON),epsilon,q,m,lambda,NULL);
      if(testval < minval) {
	minval = testval;
      } else if(testval > maxval) {
	maxval = testval;
      }
    }
    noiseest = maxval - minval;
    term = INFINITY; // make sure the inner loop starts
  }
  /* search inward to find the positive region */
  while(aval < 0) {
    if(a - badbound < HALF_POS_SEARCH_LIMIT)
      return 0;
    Real attempt = 0.5*(a + badbound);
    Real attemptval = radialfrac_half(attempt, epsilon, q, m, lambda, &term);
    if(term * TERM_NOISE > Fabs(attemptval)) {
      badbound = attempt;
    } else {
      a = attempt;
      aval = attemptval;
    }
  }

  c = Fabs(epsilon);
  for(;;) {
    cval = radialfrac_half(c, epsilon, q, m, lambda, NULL);

    if(cval < 0)
      break;
    c *= 2;
  }

  BRENT(radialfrac_half(c, epsilon, q, m, lambda, NULL));
  *imnu = c;
  return 1;
}


#define SECANT_SEARCH()							\
  do {									\
    for(n=0;n<RENANGMOM_SECANT_MAX_STEPS;n++) {				\
      Complex valrat;							\
      if(cval == aval) {						\
	valrat = 0;							\
	next = (c + a) / 2;						\
      } else {								\
	valrat = aval / (cval - aval);					\
	next = a - real(valrat) * (c - a);				\
      }									\
									\
      if(next < badbound || next > maxdist)				\
	break;								\
									\
      c = a;								\
      cval = aval;							\
      a = next;								\
      aval = cradialfrac(1.0 + a*I, epsilon, q, m, lambda, NULL);	\
									\
      /* better test? */						\
      if(Fabs(c-a) < Fabs(a)*RENANGMOM_SECANT_MIN_REL_STEP) {		\
	/* sanity check: are we trying to shoot off into complex space? */ \
	if(Fabs(imag(valrat)) > SECANT_IMAG_LIMIT*Fabs(real(valrat)))	\
	  break;							\
	found = 1;							\
	break;								\
      }									\
    }									\
  } while(0)

/* calculates nu assuming it is on the line Re(nu) = 1 (so any int is
 * then a solution)
 * Returns 1 if it suceeds, 0 else
 * Uses a projective secant method */
int FT_PRFX renangmom_iint(Real epsilon, Real q, int l, int m, Real lambda,
		       Real *imnu) {
  const Real maxdist = RENANGMOM_MAX_DIST * Fabs(epsilon);
  const Real skeptic = SECANT_SKEPTIC_MULT * nu_predictor(epsilon,q,m,lambda);

  Real badbound = RENANGMOM_DETECT_LIMIT;
  /* SECANT_SEARCH() variables */
  int n;
  Real a, c, next;
  Complex aval, cval;

  Real holdval = -1;

  int divisions, location;

  for(divisions=2;divisions<=RENANGMOM_IINT_MAX_DIVISIONS;divisions*=2) {
    for(location=1;location<divisions;location+=2) {
      c = (maxdist-badbound)*(location+1)/divisions + badbound;
      cval = cradialfrac(1.0+I*c, epsilon, q, m, lambda, NULL);
      a = (maxdist-badbound)*location/divisions + badbound;
      aval = cradialfrac(1.0+I*a, epsilon, q, m, lambda, NULL);

      int found = 0;
      SECANT_SEARCH();
      if(found &&
	 (a < 2*IINT_NOISE_SEP_MIN || a - badbound > IINT_NOISE_SEP_MIN)) {
	/* check that this looks like an expected answer */
	if(a < skeptic) {
	  if(Fabs(a - holdval) < Fabs(a)*RENANGMOM_SKEPTIC_REPEAT_SPACING) {
	    /* found the same result twice, pass it to the noise checker */
	  } else {
	    found = 0;
	    /* we don't tend to mistakenly return answers too large */
	    if(a > holdval)
	      holdval = a;
	  }
	}
	if(found) {
	  /* require there to ba a good point half way between the proposed
	   * root and the last known bad point.  This can fail on the right
	   * answer, but if the routine returns it again the check will
	   * probably pass the second time. */ 
	  Real testpoint = 0.5*(a + badbound);
	  Real term;
	  Complex testval =
	    cradialfrac(1.0 + testpoint*I, epsilon, q, m, lambda, &term);
	  if(term * TERM_NOISE < abs(testval)) {
	    *imnu = a;
	    return 1;
	  } else {
	    badbound = testpoint;
	  }
	}
      }
    }
  }

  return 0;
}

Real FT_PRFX nu_predictor(Real epsilon,Real q,int m,Real lambda) {
  Real epsilonsq = epsilon*epsilon;
  Real lamb_term = lambda + 4 - epsilonsq - epsilon*(epsilon-m*q);
  Real eps_term = 4*epsilon*(epsilon-m*q)*(4+epsilonsq);

  return imag(sqrt(0.25+0.5*(lamb_term +
			     sqrt((Complex)(lamb_term*lamb_term-eps_term)))));
}

int FT_PRFX renangmom_far(Real epsilon, Real q, int l, int m, Real lambda,
			  Complex *nu) {
  Real imnu;
  if(renangmom_half(epsilon, q, l, m, lambda, &imnu)) {
    *nu = - 0.5 + I*imnu;
    return 1;
  } else if(renangmom_iint(epsilon, q, l, m, lambda, &imnu)) {
    *nu = 1.0 + I*imnu;
    return 1;
  } else if(renangmom_real_search_procedure(epsilon, q, l, m, lambda, nu)) {
    return 1;
  }
  return 0;
}

/* search for the root assuming the function is well behaved near guess and
 * that the root is closer to guess than to a half integer */
int FT_PRFX renangmom_real_guess(Real epsilon, Real q, int m, Real lambda,
			     Real guess, Real *nu) {
  Real a, c;
  Real aval, cval;

  Real guessval = radialfrac(guess, epsilon, q, m, lambda, NULL);

  if(guessval < 0) {
    a = guess;
    aval = guessval;

    c = (2*guess + ceil(2*guess))/4;
    cval = radialfrac(c, epsilon, q, m, lambda, NULL);
    if(cval < 0) return 0;
  } else {
    c = guess;
    cval = guessval;

    a = (2*guess + floor(2*guess))/4;
    aval = radialfrac(a, epsilon, q, m, lambda, NULL);
    if(aval > 0) return 0;
  }

  BRENT(radialfrac(c, epsilon, q, m, lambda, NULL));

  /* sanity check: does this look like a root? */
  if(Fabs(cval) > Fabs(guessval)) return 0;

  *nu = c;
  return 1;
}

Real FT_PRFX fractions_poly_approx(Real epsilon, Real q, int m, Real lambda) {
  Real epsilonsq = epsilon * epsilon;
  Real kappasq = 1 - q*q;

  Real lambdastuff = lambda + 2 - 2*epsilonsq + epsilon*m*q;
  Real forth = -1.5 - 0.5*epsilonsq*kappasq - lambdastuff;
  Real second =
    0.5625
    + (epsilonsq - epsilon*m*q) * (4 + 0.5 * (epsilonsq + epsilon*m*q))
    - epsilonsq * kappasq * (epsilonsq - 3.75 + 0.125 * epsilonsq * kappasq)
    + (1.25 - 0.5 * epsilonsq*kappasq) * lambdastuff;

  Real discriminant = forth*forth - 4*second;
  if(discriminant < 0) return NAN;

  Real rootsq = (sqrt(discriminant) - forth) / 2;
  if(rootsq < 0) return NAN;

  return sqrt(rootsq) - 0.5;
}
