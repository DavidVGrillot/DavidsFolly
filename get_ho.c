//
// Created by David on 7/30/2018.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "globals.h"
#include "get_red_mass.h"
#include "get_integral.h"
// #include "read_write.h"
// #include "get_integral.h"

/* put subroutine to be tested here */


 void get_ho() {

/* based on following article: */
/* https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential */

    int i;
    int k;
    double x;
    double fact1;        /* factorial values of various integer quantities */
    double fact2;
    double fact3;
    double powe;
    double lag;  /* Laguerre polynomial */
    double integrand[nmesh];
    double temp1;
    double temp2;

/* calculate the r values for Laguerre polynomial  but first initialize r[i] to zero */

    for (i = 0; i < nmax; i++) {
        r[i] = 0;
    }

/* get the factorial terms */

    fact1 = 1;
    fact2 = 1;
    fact3 = 1;                     /* zero factorial */

    for (i = 1; i <= (n - l); i++) fact1 = fact1 * i;   /* N - L */
    for (i = 1; i <= (n + l); i++) fact2 = fact2 * i;
    for (i = 1; i <= (n + l + 1); i++) fact3 = fact3 * i;

    lag = 0.0;
    k = (n - l) / 2;
    for (i = 0; i < nmax; i++) {
        x = (i * mesh);

        if (k == 0) {
            lag = 1;
        } else if (k == 1) {
            lag = 1.0 - x;
        } else if (k == 2) {
            lag = (1.0 / 2.0) * (x * x - 4.0 * x + 2, 0);
        } else if (k == 3) {
            lag = (1.0 / 6.0) * (-x * x * x + 9.0 * x * x - 18.0 * x + 6.0);
        } else if (k == 4) {
            lag = (1.0 / 24.0) * (x * x * x * x - 16.0 * x * x * x + 72.0 * x * x - 96.0 * x + 24.0);
        } else if (k == 5) {
            lag = (1.0 / 120.0) *
                  (-x * x * x * x * x + 25.0 * x * x * x * x - 200.0 * x * x * x + 600.0 * x * x - 600.0 * x + 120.0);
        } else if (k == 6) {
            lag = (1.0 / 720.0) *
                  (x * x * x * x * x * x - 36.0 * x * x * x * x * x + 450.0 * x * x * x * x - 2400.0 * x * x * x +
                   5400.0 * x * x -
                   4320.0 * x + 720.0);
        }

        powe = -.5 * x * x;
        temp1 = pow(x, (l + 1.));
        temp2 = pow(2.7183, powe);
        r[i] = temp1 *temp2* lag;
    }

    for (i = 0; i < nmax; i++) {
        integrand[i] = r[i] * r[i];          /* this is to normalize the wavefunction to one */
    }

    get_integral(integrand);                 /* calculates tmat the squae o f the normalization constant */

    for (i = 0; i < nmax; i++) {
        r[i] = r[i]/sqrt(tmat);             /* put in terms of rho */
    }
    return;
}

