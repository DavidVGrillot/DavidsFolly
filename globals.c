// globals.c
#include "globals.h"
 /* following input variables are global to all functions in this file */
char a_type[3];             /* if bg else ei eigenvalue */
char b_type[5];            /* type baryon xi+, sig+, del+, pro+, neu0 */
int l;                     /* Har. osc. ang mom. quantum num of trial sol */
int i_pot;                  /* indicates the type of n-n potential 2 + types ?? */
int n;                      /* Harmonic osc quantum number of trial solution */
int sp;                     /* spin  vector either 0 or 1 */
int tau;                    /* isospin vector either 0 or 1 */
int j;                      /* total angular  momentum j = l + s */
double eng;                /* hbar*omega single particle energy parameter to vary */
/* **************end of input variable names */
/* following are globabl vars */

int nmesh = 80;           /* number of mess points */
int nmax;
double mu[3];              /* reduced masses of three quarks in MeV/c */
double muu;                /* the real red. mass used in calculations */
double tmat;               /* tmatrix <t|v|u> ; */
/* N.B. the reduced mass is always put into mu(1) */
double rho;                /* dimensionless position value analog of r ro x */
double *u;                 /*solu. wave fct*/
double *r;                 /*solu. harmonic oscill. wave fct*/
double *pot;
double mesh=.485;          /* mesh width for calculations */
double gam;                /* converts fentometers to to dimensionless  */
double pi=3.14159;
double hbar;  /* planck constant divided by 2pi in units of mev sec */
