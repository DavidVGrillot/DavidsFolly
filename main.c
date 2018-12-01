#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "globals.h"
#include "get_red_mass.h"
#include "get_integral.h"
#include "get_pot.h"
#include "get_ho.h"
#include "bg_solver.h"
#include "read_write.h"
#include "bg_solver_2.h"
// void commons_to_angmom_();
// double gmosh_(int n,int l,int nc, int lc, int n1,int l1,int n2,int l2,int lr,double d);

void main() {

    double c;
    nmax = nmesh - 1;
    u = (double *) calloc(nmesh, sizeof(double));    /*solu. wave fct*/
    r = (double *) calloc(nmesh, sizeof(double));    /*solu. wave fct*/
    pot = (double *) calloc(nmesh, sizeof(double));  /* N-N potential box or HJ */
    hbar = 6.582 * pow(10, -22);
    c = 3 * pow(10, 10);
    gam = muu * eng / (c * c * hbar*hbar);  /* converts x(schroedinger) to rho(dimensionless) */
    gam = 1;

    //  commons_to_angmom_();
    //  double d = gmosh_(0,0,0,0,0,0,0,0,0, 0.);

    /* Read input data file: a_type, q_type_i, b_type, ipot, n, l, s, tau, eng,  */

    FILE *io = fopen("../data_input", "r+w");
    read_file(io);
//    b_type = "pro+"
    j = l + sp;                     /* total angular momentum */


/* CASE 1 m1 and m2  */

/* set the reduced mass to mu(1) */

    get_red_mass();                 /* get the reduced masses and charges: of the quark input particle type */
    muu = mu[0];


//    mesh = sqrt(gam)*mesh;             /* convert mesh to dimensionless */
    mesh = mesh/10;

/* call numerical Bethe-Goldstone or Eigenvalue _solver */

   if(strcmp(a_type, "bg") ==0) {

       bg_solver(0);

   }else {
       bg_solver_2();
   }



/* CASE 2 m1 and m3  */
/* set the reduced mass in m(2) = m(1) */

    muu = mu[1];

//   bg_solver(); */

/* output to file  and exit */
/* need code here ???????????????? */

/* CASE 3 m2 and m3 */
/* set the reduced mass in m(3) = m(1) */
    muu = mu[2];


/* bg_solver(); */

/* output to file  and exit */
/* need code here ???????????????? */

//    const int c_results_size = nmesh;
//    double resultsArray[c_results_size];
//    for(int i = 0; i < c_results_size; ++i)
//    resultsArray[i] = r[i];

//    write_file(resultsArray, c_results_size, io);
//     write_file(r, nmesh, io);
//     write_file(u, nmesh, io);

    fclose(io);

    return;
}