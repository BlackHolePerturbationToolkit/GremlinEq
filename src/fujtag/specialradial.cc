//---------------------------------------------------------------------------
//
// $Id: specialradial.cc,v 1.3 2014/08/01 17:16:59 sahughes Exp $
//
//---------------------------------------------------------------------------

/* specialradial.cc
 * functions to compute properties of the radial continued fraction
 * at special values
 */

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "specialradial.h"
#endif
#include <cmath>

#define SPECIALRADIAL_SMALLNUMBER 1e-30


/* useful macros */
#define DECLARE_ALPHAGAMMA()			\
  Real macroa1,macroa2,macroa3
#define alphagamma(n)							\
  (macroa1 = (n)-1.5,							\
   macroa2 = (n)+2.5,							\
   macroa3 = (n)+0.5,							\
   (ag_mult*((n)-0.5)*((n)+1.5)*((n)-1)*((n)+2)				\
    *(macroa1*macroa1 + epsilonsq)*(macroa2*macroa2 + epsilonsq)	\
    *(macroa3*macroa3 + tausq)))
#define DECLARE_BETA()			\
  Real macrob1
#define beta(n)							\
  (macrob1 = ((n)*(n)-0.25),					\
   (4*(macrob1*(macrob1 + b_add1) + b_add2)*((n)*(n) - 1)))

/* Edited SAH 1 Aug 2014 --- getting compiler warnings */
// #define DECLARE_DALPHAGAMMA()					\
//   Real macroda1,macroda2,macroda3,macroda4,macroda5,macroda6
#define dalphagamma(n)							\
  (2*ag_mult*(2*((n) + 0.5)*((n)*((n) + 1) - 1.375)*			\
	      (((n) - 1.5)*((n) - 1.5) + epsilonsq)*			\
	      (((n) + 2.5)*((n) + 2.5) + epsilonsq)*			\
	      (((n) + 0.5)*((n) + 0.5) + tausq) +			\
	      ((n) - 0.5)*((n) + 1.5)*((n) - 1)*((n) + 2)		\
	      *(((n) - 1.5)*(((n) + 2.5)*((n) + 2.5) + epsilonsq)	\
		*(((n) + 0.5)*((n) + 0.5) + tausq) +			\
		(((n) - 1.5)*((n) - 1.5) + epsilonsq)*			\
		((n) + 2.5)*(((n) + 0.5)*((n) + 0.5) + tausq) +		\
		(((n) - 1.5)*((n) - 1.5) + epsilonsq)*                  \
		(((n) + 2.5)*((n) + 2.5) + epsilonsq)*((n) + 0.5))))
// #define dalphagamma(n)							\
//   (macroda1 = (n)-1.5,							\
//    macroda2 = (n)+2.5,							\
//    macroda3 = (n)+0.5,							\
//    macroda4 = macroda1*macroda1 + epsilonsq,				\
//    macroda5 = macroda2*macroda2 + epsilonsq,				\
//    macroda6 = macroda3*macroda3 + tausq,				\
//    2*ag_mult*(2*macroda3*((n)*((n)+1)-1.375)*macroda4*macroda5*macroda6 + \
// 	      ((n)-0.5)*((n)+1.5)*((n)-1)*((n)+2)			\
// 	      *(macroda1*macroda5*macroda6 +				\
// 		macroda4*macroda2*macroda6 +				\
// 		macroda4*macroda5*macroda3)))

#define DECLARE_DBETA()
#define dbeta(n)				\
  (8*(n)*(3*((n)*(n)-0.75)*((n)*(n)-0.25) +	\
	  2*((n)*(n)-0.625)*b_add1 + b_add2))

/* The followint two macros are only for n=1 */
/* Edited SAH 1 Aug 2014 --- getting compiler warnings */
// #define DECLARE_DDALPHAGAMMA1()			\
//   Real macrodda4,macrodda5,macrodda6
#define ddalphagamma1					\
  (0.25*ag_mult*(41*(0.25 + epsilonsq)*                 \
(12.25 + epsilonsq)*(2.25 + tausq)			\
 + 15.*(-(12.25 + epsilonsq)*(2.25 + tausq)		\
       + 7.*(0.25 + epsilonsq)*(2.25 + tausq)		\
	+ 3.*(0.25 + epsilonsq)*(12.25 + epsilonsq))))
// #define ddalphagamma1					\
//   (macrodda4 = 0.25 + epsilonsq,			\
//    macrodda5 = 12.25 + epsilonsq,			\
//    macrodda6 = 2.25 + tausq,				\
//    0.25*ag_mult*(41*macrodda4*macrodda5*macrodda6	\
// 		 + 15*(-  macrodda5*macrodda6		\
// 		       +7*macrodda4*macrodda6		\
// 		       +3*macrodda4*macrodda5)))
#define DECLARE_DDBETA1()
#define ddbeta1 (26.25 + 19*b_add1 + 4*b_add2)

#define GET_PARTIAL(n)					\
  do {							\
    Real a,b,d,delta;					\
    int loop;						\
    partial_frac = beta(n);				\
							\
    a = alphagamma(n);					\
    b = beta((n)+1);					\
    d = 1/b;						\
    delta = -a*d;					\
							\
    partial_frac += delta;				\
							\
    for(loop=(n)+2;;loop++) {				\
      if(Fabs(a) < SPECIALRADIAL_SMALLNUMBER) break;	\
      							\
      a = alphagamma(loop-1);				\
      b = beta(loop);					\
      							\
      delta *= a*d;					\
      d = 1/(b - a * d);				\
      delta *= d;					\
      							\
      if(partial_frac == partial_frac + delta) break;	\
      partial_frac += delta;				\
      if(isnan(partial_frac))				\
	return partial_frac;				\
    }							\
  } while(0)


/* calculate the value of the function at nu = -1/2 */
Real FT_PRFX radial_half_value(Real epsilon, Real q, int m, Real lambda,
			       Real *term) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q)/(1 - q*q);

  /* constants for alpha_n gamma_{n+1} */
  const Real ag_mult = 4*epsilonsq*(1 - q*q);
 
  /* constants for beta_n */
  /* the additive things for the inner calulation */
  const Real b_add1 = 2*epsilonsq - epsilon*m*q - lambda - 2;
  const Real b_add2 = epsilon*(epsilon - m*q)*(4 + epsilonsq);

  DECLARE_ALPHAGAMMA();
/* Edited SAH 1 Aug 2014 --- getting compiler warnings */
  // DECLARE_DALPHAGAMMA();
  // DECLARE_DDALPHAGAMMA1();
  DECLARE_BETA();
  DECLARE_DBETA();
  DECLARE_DDBETA1();

  int n;
  Real frac;

  Real termmult;
  Real partial_frac;

  /* add beta(0) at the end to avoid calculating needless precision
   * (the full value vanishes more often than the perturbation term) */

  frac = 0;
  GET_PARTIAL(2);
  partial_frac = 1/partial_frac;

  termmult = -2.0/(dbeta(1)-dalphagamma(1)*partial_frac);
  frac += termmult * dalphagamma(0);
  if(term) *term = frac*frac;

  termmult *= 0.5 * alphagamma(0) * termmult; /* this is not a typo */
  frac += termmult * ddbeta1;
  if(term) *term = termmult*termmult*ddbeta1*ddbeta1;

  termmult *= -partial_frac;
  frac += termmult * ddalphagamma1;
  if(term) *term = termmult*termmult*ddalphagamma1*ddalphagamma1;

  termmult *= -dalphagamma(1)*partial_frac;
  //if(term) *term = Fabs(frac);

  for(n=2;;n++) {
    Real delta;

    delta = termmult * dbeta(n);
    if(term) *term += delta*delta;

    GET_PARTIAL(n+1);
    partial_frac = 1/partial_frac;
    termmult *= -partial_frac;

    delta += termmult * dalphagamma(n);
    if(term) *term += termmult*termmult*dalphagamma(n)*dalphagamma(n);

    if(frac == frac + delta) break;
    frac += delta;
    //if(term && Fabs(delta)>*term) *term = Fabs(delta);
    if(isnan(frac)) /* something went horribly wrong */
      return frac;

    termmult *= -alphagamma(n)*partial_frac;
  }

  Real b0 = beta(0);
  if(term) *term = sqrt(*term + b0*b0);
  return frac + b0;
}


/* redefine the macros for nu = 1 */
#undef DECLARE_ALPHAGAMMA
#undef alphagamma
#undef DECLARE_BETA
#undef beta
#undef DECLARE_DALPHAGAMMA
#undef dalphagamma
#undef DECLARE_DBETA
#undef dbeta
#undef DECLARE_DDALPHAGAMMA1
#undef ddalphagamma1
#undef DECLARE_DDBETA1
#undef ddbeta1

#define DECLARE_ALPHAGAMMA()
#define alphagamma(n)						\
  (ag_mult*((n)+1)*((n)+3)*((n)+0.5)*((n)+3.5)			\
   *((n)*(n) + epsilonsq)*(((n)+4)*((n)+4) + epsilonsq)		\
   *(((n)+2)*((n)+2) + tausq))
#define DECLARE_BETA()				\
  Real macrob1
#define beta(n)								\
  (macrob1 = ((n) + 1)*((n) + 2),					\
   (4*(macrob1*(macrob1 + b_add1) + b_add2)*((n) + 2.5)*((n) + 0.5)))

#define DECLARE_DALPHAGAMMA()				\
  Real macroda2,macroda3,macroda4,macroda5,macroda6
#define dalphagamma(n)							\
  (macroda2 = (n)+4,							\
   macroda3 = (n)+2,							\
   macroda4 = (n)*(n) + epsilonsq,					\
   macroda5 = macroda2*macroda2 + epsilonsq,				\
   macroda6 = macroda3*macroda3 + tausq,				\
   2*ag_mult*(2*macroda3*((n)*macroda2+2.375)*macroda4*macroda5*macroda6 + \
	      ((n)+1)*((n)+3)*((n)+0.5)*((n)+3.5)			\
	      *((n)*macroda5*macroda6 +					\
		macroda4*macroda2*macroda6 +				\
		macroda4*macroda5*macroda3)))
#define DECLARE_DBETA()
#define dbeta(n)					\
  (8*((n)+1.5)*(3*((n)*((n)+3)+1.5)*((n)+1)*((n)+2) +	\
		2*((n)*((n)+3)+1.625)*b_add1 + b_add2))

#define DECLARE_DDALPHAGAMMA_1()
#define ddalphagamma_1						\
  (0.25*ag_mult*(11*(1+epsilonsq)*(1+epsilonsq)*(1+tausq)	\
		 + 48*(1+epsilonsq)*(1+tausq)			\
		 - 20*(1+epsilonsq)*(1+epsilonsq)		\
		 +160*(tausq-epsilonsq)))
#define DECLARE_DDALPHAGAMMA_2()
#define ddalphagamma_2							\
  (0.25*ag_mult*(9*(4+epsilonsq)*(4+epsilonsq)				\
		 - (epsilonsq*(13*epsilonsq + 86) + 280)*tausq))
#define DECLARE_DDBETA_1()
#define ddbeta_1 (-3 + b_add1 + 4*b_add2)


/* calculate the slope of the function at nu = 1 */
Real FT_PRFX radial_int_slope(Real epsilon, Real q, int m, Real lambda,
			      Real *term) {
  const Real epsilonsq = epsilon*epsilon;
  const Real tausq = (epsilon - m*q)*(epsilon - m*q)/(1 - q*q);

  /* constants for alpha_n gamma_{n+1} */
  const Real ag_mult = 4*epsilonsq*(1 - q*q);

  /* constants for beta_n */
  /* the additive things for the inner calulation */
  const Real b_add1 = 2*epsilonsq - epsilon*m*q - lambda - 2;
  const Real b_add2 = epsilon*(epsilon - m*q)*(4 + epsilonsq);

  DECLARE_ALPHAGAMMA();
  DECLARE_BETA();
  DECLARE_DALPHAGAMMA();
  DECLARE_DBETA();
  DECLARE_DDALPHAGAMMA_1();
  DECLARE_DDALPHAGAMMA_2();
  DECLARE_DDBETA_1();

  int n;
  Real frac;

  Real termmult;
  Real partial_frac;

  /* add the closed form term at the end to avoid calculating needless
   * precision
   * (the full value vanishes more often than the perturbation term) */

  frac = 0;
  termmult = 2;

  for(n=0;;n++) {
    Real delta;

    delta = termmult * dbeta(n);

    GET_PARTIAL(n+1);
    partial_frac = 1/partial_frac;
    termmult *= -partial_frac;

    delta += termmult * dalphagamma(n);

    if(frac == frac + delta) break;
    frac += delta;
    if(isnan(frac)) /* something went horribly wrong */
      return frac;

    termmult *= -alphagamma(n)*partial_frac;
  }

  Real dag_1 = dalphagamma(-1);
  Real b_1i = 1/beta(-1);
  Real db_1 = dbeta(-1);

  GET_PARTIAL(0);

  Real nonfrac =
    -b_1i * dag_1
    + partial_frac * (2 * b_1i * db_1 - 2 * ddalphagamma_1 / dag_1
		      + partial_frac / dag_1 *
		      (2*ddbeta_1 - b_1i * (ddalphagamma_2 + db_1*db_1)));

  if(term) *term = Fabs(nonfrac);
  return frac + nonfrac;
}
