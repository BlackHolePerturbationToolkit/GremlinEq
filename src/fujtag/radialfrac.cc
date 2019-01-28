//---------------------------------------------------------------------------
//
// $Id: radialfrac.cc,v 1.2 2014/04/08 19:36:44 sahughes Exp $
//
//---------------------------------------------------------------------------

/* radialfrac.cc
 * function to evaluate the continued fraction
 */

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "radialfrac.h"
#endif

#include <cmath>
#include "fraction_macros.h"

/* set to 1 if the fraction evaluators should use more precision
 * internally to avoid some loss of precision */
#define RADIALFRAC_MORE_PRECISION 1

#define RADIALFRAC_SMALLNUMBER 1e-30

#define MAKE_FRAC_FUNCS(type,internal_prefix,prefix,smalltest,badtest,absfunc,badans) \
  type FT_PRFX prefix ## plusfrac(type nu, Real epsilon, Real q, int m,	\
				  Real lambda, Real *term) {		\
    internal_prefix ## type frac;					\
    int n;								\
    internal_prefix ## type a,b,d,delta;				\
    									\
    const internal_prefix ## Real epsilonsq = epsilon*epsilon;		\
    const internal_prefix ## Real tausq =				\
      (epsilon - m*q)*(epsilon - m*q)/(1 - q*q);			\
    									\
    PREP_ALPHAGAMMA(internal_prefix ## type);				\
    PREP_BETA(internal_prefix ## type);					\
    									\
    a = alphagamma(0);							\
    b = beta(1);							\
    d = 1.0 / b;							\
    delta = -a*d;							\
    									\
    frac = delta;							\
    if(term) *term = absfunc(delta);					\
    									\
    for(n=2;;n++) {							\
      if(smalltest) return (type) frac;					\
      									\
      a = alphagamma(n-1);						\
      b = beta(n);							\
      									\
      delta *= a*d;							\
      d = 1.0 / (b - a * d);						\
      delta *= d;							\
      									\
      if(frac == frac + delta) return (type) frac;			\
      frac += delta;							\
      if(term) {							\
	Real absdelta = absfunc(delta);					\
	Real absfrac = absfunc(frac);					\
	if(absdelta > *term) *term = absdelta;				\
	if(absfrac > *term) *term = absfrac;				\
      }									\
      if(badtest) /* something went horribly wrong */			\
	return badans;							\
    }									\
  }									\
									\
  type FT_PRFX prefix ## minusfrac(type nu, Real epsilon, Real q, int m, \
				   Real lambda, Real *term) {		\
    internal_prefix ## type frac;					\
    int n;								\
    internal_prefix ## type a,b,d,delta;				\
									\
    const internal_prefix ## Real epsilonsq = epsilon*epsilon;		\
    const internal_prefix ## Real tausq =				\
      (epsilon - m*q)*(epsilon - m*q)/(1 - q*q);			\
									\
    PREP_ALPHAGAMMA(internal_prefix ## type);				\
    PREP_BETA(internal_prefix ## type);					\
									\
    a = alphagamma(-1);							\
    b = beta(-1);							\
    d = 1.0 / b;							\
    delta = -a*d;							\
									\
    frac = delta;							\
    if(term) *term = absfunc(delta);					\
									\
    for(n=-2;;n--) {							\
      if(smalltest) return (type) frac;					\
									\
      a = alphagamma(n);						\
      b = beta(n);							\
									\
      delta *= a*d;							\
      d = 1.0 / (b - a * d);						\
      delta *= d;							\
									\
      if(frac == frac + delta) return (type) frac;			\
      frac += delta;							\
      if(term) {							\
	Real absdelta = absfunc(delta);					\
	Real absfrac = absfunc(frac);					\
	if(absdelta > *term) *term = absdelta;				\
	if(absfrac > *term) *term = absfrac;				\
      }									\
      if(badtest) /* something went horribly wrong */			\
	return badans;							\
    }									\
  }									\
									\
  type FT_PRFX prefix ## radialfrac(type nu, Real epsilon, Real q, int m, \
				    Real lambda, Real *term) {		\
    /* answer */							\
    type frac;								\
    Real fracterm;							\
    Real *fractermptr = term ? &fracterm : NULL;			\
									\
    const Real epsilonsq = epsilon*epsilon;				\
									\
    PREP_BETA(type);							\
									\
    frac = beta(0);							\
    if(term) *term = absfunc(frac);					\
									\
    frac += prefix ## plusfrac(nu, epsilon, q, m, lambda, fractermptr);	\
    if(term) {								\
      if(fracterm > *term) *term = fracterm;				\
      Real absfrac = absfunc(frac);					\
      if(absfrac > *term) *term = absfrac;				\
    }									\
    if(badtest)								\
      return badans;							\
    									\
    frac += prefix ## minusfrac(nu, epsilon, q, m, lambda, fractermptr); \
    if(term) {								\
      if(fracterm > *term) *term = fracterm;				\
      Real absfrac = absfunc(frac);					\
      if(absfrac > *term) *term = absfrac;				\
    }									\
									\
    return frac;							\
  }

#if RADIALFRAC_MORE_PRECISION
#  if LREAL_MANT_DIG <= REAL_MANT_DIG
#    warning "Increased internal precision was requested for fraction" \
  "evaluators, but external precision is the highest supported by the compiler"
#  endif

MAKE_FRAC_FUNCS(Real,Long,,
		Fabs(a) < RADIALFRAC_SMALLNUMBER,
		isnan(frac),Fabs,NAN)
MAKE_FRAC_FUNCS(Complex,Long,c,
		norm(a) < RADIALFRAC_SMALLNUMBER*RADIALFRAC_SMALLNUMBER,
		Cisnan(frac),abs,CNAN)
#else /* RADIALFRAC_MORE_PRECISION */
MAKE_FRAC_FUNCS(Real,,,
		Fabs(a) < RADIALFRAC_SMALLNUMBER,
		isnan(frac),Fabs,NAN)
MAKE_FRAC_FUNCS(Complex,,c,
		norm(a) < RADIALFRAC_SMALLNUMBER*RADIALFRAC_SMALLNUMBER,
		Cisnan(frac),abs,CNAN)
#endif /* RADIALFRAC_MORE_PRECISION */

/* Evaluate the continued fraction on the special line Re(nu) = -1/2
 * It is guaranteed to be real there */
Real FT_PRFX radialfrac_half(Real imnu, Real epsilon, Real q, int m,
			     Real lambda, Real *term) {
  const Real epsilonsq = epsilon*epsilon;

  /* constants for beta_n */
  /* offsets for the inner n calculation */
  const Complex b_off_in1 = imnu*I - 0.5;
  const Complex b_off_in2 = imnu*I + 0.5;
  /* the outer roots */
  const Complex b_rt1 = imnu*I + 1.0;
  const Complex b_rt2 = imnu*I - 1.0;
  /* the additive things for the inner calulation */
  const Real b_add1 = 2*epsilonsq - epsilon*m*q - lambda - 2;
  const Real b_add2 = epsilon*(epsilon - m*q)*(4 + epsilonsq);


  /* We know the fractions are complex conjugates of one another, so can just
   * do one evaluation */

  Real b0 = real(beta(0));
  if(term) *term = Fabs(b0);
  Real fracterm;
  Real *fractermptr = term ? &fracterm : NULL;
  Real frac =
    b0 + 2.0*real(cplusfrac(-0.5 + imnu*I,epsilon,q,m,lambda,fractermptr));
  if(term) {
    if(2*fracterm > *term) *term = 2*fracterm;
    Real absfrac = Fabs(frac);
    if(absfrac > *term) *term = absfrac;
  }
  return frac;
}
