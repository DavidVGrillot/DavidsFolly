// globals.h

/* following input variables are global to all functions in this file */
extern char a_type[3];             /* if bg else ei eigenvalue */
extern char b_type[5];            /* type baryon xi+, sig+, del+, pro+, neu0 */
extern int l;                     /* Har. osc. ang mom. quantum num of trial sol */
extern int i_pot;                  /* indicates the type of n-n potential 2 + types ?? */
extern int n;                      /* Harmonic osc quantum number of trial solution */
extern int sp;                     /* spin  vector either 0 or 1 */
extern int tau;                    /* isospin vector either 0 or 1 */
extern int j;                      /* total angular  momentum j = l + s */
extern double eng;                /* hbar*omega single particle energy parameter to vary */
/* **************end of input variable names */
/* following are globabl vars */

extern int nmesh;                 /* number of mess points */
extern int nmax;
extern double mu[3];              /* reduced masses of three quarks in MeV/c */
extern double muu;                /* the real red. mass used in calculations */
/* N.B. the reduced mass is always put into mu(1) */
extern double rho;                /* dimensionless position value analog of r ro x */
extern double *u;     /*solu. wave fct*/
extern double *r;     /*solu. harmonic oscill. wave fct*/
extern double *pot;  /* N-N potential box or HJ */
extern double mesh;          /* mesh width for calculations */
extern double gam;
extern double pi;
extern double hbar;  /* planck constant divided by 2pi in units of mev sec */
extern double tmat;  /* t matrix  <r|v|u>  */
