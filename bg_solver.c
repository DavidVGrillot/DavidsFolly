//
// Created by David on 7/30/2018.
//

/* ******************************************** */

/* Subroutine BG_solver */
#include <stdio.h>
#include "globals.h"
#include "get_red_mass.h"
#include "get_integral.h"
#include "get_pot.h"
#include "get_ho.h"
// #include "bg_solver.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include "read_write.h"

void bg_solver() {
/*  program to solve either eigenvalue or Bethe-Goldstone equation for
      wave function u */
/* output is the calculated final wave function, u and t-matrix */

//    double *u     =  (double*) calloc((size_t)nmesh, sizeof(double));

/* ******************************************************* */
/* see notes for definition of following terms */

    double *a = (double *) calloc((size_t) nmesh, sizeof(double));
    double *b = (double *) calloc((size_t) nmesh, sizeof(double));
    double *c = (double *) calloc((size_t) nmesh, sizeof(double));
    double *d1 = (double *) calloc((size_t) nmesh, sizeof(double));  /* small d in notes */
    double *d2 = (double *) calloc((size_t) nmesh, sizeof(double));  /* big d  in notes */
    double *d3 = (double *) calloc((size_t) nmesh, sizeof(double));
    /* test thingy */
    double *e = (double *) calloc((size_t) nmesh, sizeof(double));
    double *f = (double *) calloc((size_t) nmesh, sizeof(double));
    double *h = (double *) calloc((size_t) nmesh, sizeof(double));
    double *q = (double *) calloc((size_t) nmesh, sizeof(double));

    double *pot = (double *) calloc((size_t) nmesh, sizeof(double));  /* n-n potential */
    double *delu = (double *) calloc((size_t) nmesh, sizeof(double));  /* wfc change */
    double *pctd1 = (double *) calloc((size_t) nmesh, sizeof(double));
    double *lamda = (double *) calloc((size_t) nmesh, sizeof(double));  /* energy parm  */
    double sumdelu;                                                   /* sum of delus */
    double coef;         /*used in Pauli */
    int i;
    int j;
    double dellam;  /* change in the eigenvalue at end of each iteration */
    double lam;         /* eigen value */
    double integrand[nmesh];
    double ru;
    double integrand1[nmesh];
    double rvu;
    double e0;          /* harmonic energy for ground state */
    double k;           /* constant to fix boundary conditions at large i */
    double *out1 = (double *) calloc((size_t) nmesh, sizeof(double));
    double *out2 = (double *) calloc((size_t) nmesh, sizeof(double));

/* initialize stuff. n.b. nmesh is set in globals.c  */
    FILE *io = fopen("../data_input", "r+w");

    get_ho();                 /* Get trial wave harmonic osc wave funtcion */
    sumdelu = 1.0;
    lam = eng;
    coef = 1.;
    j = 0;
    a[0] = 1. + (mesh * mesh / 12.) * (lam - pot[1]);
    b[0] = -2. + (5.0 / 6.0) * mesh * mesh * (lam - pot[0]);
    c[0] = 1. - (mesh * mesh / 12.) *lam;
    d1[0] = b[0];
    e0 = (n + 3.0 / 2.0) *2* lam;
    d3[0] = -(mesh * mesh / 12.) * (lam - e0) * r[1];
    q[0] = a[0] * r[1] + d3[0];
    f[0] = -(q[0] + d3[0]) / a[0];
    h[0] = -a[0] / b[0];
    u[0] = 0;

/* Basic loop dee--loop */
    while (fabs(sumdelu) > .01) {
        sumdelu = 0;
        for (i = 0; i < nmax; i++) {

            delu[i] = 0;
            if (j == 0) u[i] = r[i];  /*  first time fix the H.O. wfct for summation in Pauli excl. princ */
            integrand1[i] = r[i] * u[i];
            pot[i] = get_pot(i);              /* output is in v[] v[rho] = <Y[jlst]| v[hamada] | Y[jlst]> */
            integrand[i] = r[i] * pot[i] * u[i];
            lamda[i] = (lam - e0)*r[i]*ru + 12.0*coef*(r[i] * rvu);
        }
        get_integral(integrand);
        rvu = tmat;
        get_integral(integrand1);
        ru = tmat;
        j = j + 1;

        for (i = 1; i < nmax; i++) {
            a[i] = 1. + (mesh*mesh / 12.) * (lam - pot[i + 1]);
            b[i] = -2. + (5.0 / 6.0) * (mesh*mesh/12.) * (lam - pot[i]);
            c[i] = 1. + (mesh*mesh/12.) * (lam - pot[i - 1]);
            q[i] = (a[i] * u[i + 1] + b[i] * u[i] + c[i] * u[i - 1]);
            d1[i] = b[i] + c[i] * h[i - 1];

            pctd1[i] = (d1[i]-d1[i-1])/d1[i-1];
//            if(fabs(pctd1[i]) > 0.1) d1[i] = d1[i-1];
            out1[i] = u[i];
            d3[i]=  -(mesh * mesh / 12.)*(lamda[i + 1] + 10.0 * lamda[i] + lamda[i - 1] + 12*r[i]*rvu);
            q[i] = q[i] + d3[i];
            f[i] = -(q[i] + c[i] * f[i - 1]) / d1[i];
            h[i] = -a[i] / d1[i];

            k = pow(2.718, -mesh * mesh * (2.0 * nmax + 1) / 2);
            if (strcmp(a_type, "bg") == 0) {
                if (i == 1) {
                    delu[i] = -q[0] / a[0];
                } else if (i == 2) {
                      delu[2] = (delu[1] - f[1])/h[1];
//                      delu[2] = -(q[1] + b[1]*delu[1])/a[1];
                } else {
                      delu[i] = (delu[i - 1] - f[i - 1]) / h[i - 1];
//                    delu[i] = -(q[i - 1] + b[i - 1] * delu[i - 1] + c[i - 1] * delu[i - 2]) / a[i];
                }
            }
            delu[i] = delu[i] * f[i - 1] / (1 + k * h[i - 1]);
            u[i] = u[i] + delu[i];
            sumdelu = sumdelu + delu[i];
        }
        printf("sumdelu = %g\n", sumdelu);
        fflush(stdout);
        write_file(f, nmesh, io);
        write_file(h, nmesh, io);
        write_file(d1, nmesh, io);


        for (i = 0; i < nmax; i++) {
            integrand1[i] = u[i] * u[i];              /* normalize u */
        }
        get_integral(integrand1);

        for (i = 0; i < nmax; i++) {
            u[i] = u[i] / tmat;                       /* finish normalization */
        }
        if (j > 1) {
            printf("no convergence after tries = %d\n", j);
            fflush(stdout);
            sumdelu = 0.0;
            dellam = 0.0;
        }

        dellam = (f[nmax - 1]) / e[nmax - 1];
//        dellam= (delu[nmax-1]+h[0]*delu[1]-f[nmax-1]-f[0])/(e[nmax-1]-e[0]);
        if (strcmp(a_type, "bg") == 0) dellam = 0.0;
        lam = lam + dellam;
        printf("lam = %g\n", lam);
        fflush(stdout);
        printf("dellam = %g\n", dellam);
        fflush(stdout);
    }

    for (i = 0; i < nmax; i++) {
        integrand1[i] = u[i] * u[i];              /* normalize u */
    }
    get_integral(integrand1);

    for (i = 0; i < nmax; i++) {
        u[i] = u[i] / tmat;
        integrand[i] = r[i] * pot[i] * u[i];
    }
    get_integral(integrand);   /* recomp t-mat if new u */
    return;
}



