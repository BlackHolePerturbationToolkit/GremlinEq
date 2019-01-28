//---------------------------------------------------------------------------
//
// $Id: fraction_macros.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------

/* fractions_macros.h
 * macros related to the continued fractions
 */

#define PREP_ALPHAGAMMA(type)						\
  /* constants for alpha_n gamma_{n+1} */				\
  const Real ag_mult = 4*epsilonsq*(1 - q*q);				\
  /* the rt variables are minus the linear roots of the poly */		\
  const type ag_rt1 = nu;						\
  const type ag_rt2 = nu + 2.0;						\
  const type ag_rt3 = nu - 0.5;						\
  const type ag_rt4 = nu + 2.5;						\
  /* the qoff variables are the offsets for the quadratic terms */	\
  const type ag_qoff_e1 = nu - 1.0;					\
  const type ag_qoff_e2 = nu + 3.0;					\
  const type ag_qoff_t = nu + 1.0
#define alphagamma(n)							\
  ({ typeof(ag_qoff_e1) macroa1 = (n)+ag_qoff_e1;			\
    typeof(ag_qoff_e2) macroa2 = (n)+ag_qoff_e2;			\
    typeof(ag_qoff_t) macroa3 = (n)+ag_qoff_t;				\
    ag_mult * ((n)+ag_rt1) * ((n)+ag_rt2) * ((n)+ag_rt3) * ((n)+ag_rt4)	\
      *(macroa1*macroa1 + epsilonsq) * (macroa2*macroa2 + epsilonsq)	\
      *(macroa3*macroa3 + tausq); })

#define PREP_BETA(type)							\
  /* constants for beta_n */						\
  /* offsets for the inner n calculation */				\
  const type b_off_in1 = nu;						\
  const type b_off_in2 = nu + 1.0;					\
  /* the outer roots */							\
  const type b_rt1 = nu + 1.5;						\
  const type b_rt2 = nu - 0.5;						\
  /* the additive things for the inner calulation */			\
  const Real b_add1 = 2*epsilonsq - epsilon*m*q - lambda - 2;		\
  const Real b_add2 = epsilon*(epsilon - m*q)*(4 + epsilonsq)
#define beta(n)								\
  ({ typeof(b_off_in1) macrob1 = ((n) + b_off_in1) * ((n) + b_off_in2);	\
    4.0 * (macrob1*(macrob1 + b_add1) + b_add2)				\
      * ((n) + b_rt1) * ((n) + b_rt2); })

#define PREP_ALPHA()					\
  const Real kappa = Sqrt(1-q*q);			\
  const Complex iepskappa = I*epsilon*kappa;		\
  const Real tau = (epsilon - m*q)/kappa
#define alpha(n)							\
  ({ typeof(nu) macroa_1 = (n) + nu;					\
    typeof(macroa_1) macroa_2 = macroa_1 - 1.0;				\
    iepskappa * macroa_1 * (2.0*macroa_1 - 1.0)				\
      * (macroa_2*macroa_2 + epsilonsq) * (macroa_1 + 1.0 + I*tau); })

#define PREP_GAMMA()					\
  const Real kappa = Sqrt(1-q*q);			\
  const Complex iepskappa = I*epsilon*kappa;		\
  const Real tau = (epsilon - m*q)/kappa
#define gamma(n)							\
  ({ typeof(nu) macrog_1 = (n) + nu;					\
    typeof(macrog_1) macrog_2 = macrog_1 + 2.0;				\
    -iepskappa * (macrog_1 + 1.0) * (2.0*macrog_1 + 3.0)		\
      * (macrog_2*macrog_2 + epsilonsq) * (macrog_1 - I*tau); })

#define PREP_ALPHA_GAMMA()					\
  const Real kappa = Sqrt(1-q*q);				\
  const Complex iepskappa = I*epsilon*kappa;			\
  const Real tau = (epsilon - m*q)/kappa
