//
// Created by David on 11/20/2018.
//
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
// #include <stdlib.h>
#include "read_write.h"

#include "bg_solver_2.h"

void bg_solver_2() {
/*  program to solve either eigenvalue or Bethe-Goldstone equation for
      wave function u */
/* output is the calculated final wave function, u and t-matrix */

//    double *u     =  (double*) calloc((size_t)nmesh, sizeof(double));

/* ******************************************************* */
/* see notes for definition of following terms */

    double *a = (double *) calloc((size_t) nmesh, sizeof(double));
    double *b = (double *) calloc((size_t) nmesh, sizeof(double));
    double *c = (double *) calloc((size_t) nmesh, sizeof(double));
    double *d3 = (double *) calloc((size_t) nmesh, sizeof(double));
    double *q = (double *) calloc((size_t) nmesh, sizeof(double));
    double *pot = (double *) calloc((size_t) nmesh, sizeof(double));  /* n-n potential */
    double *delu = (double *) calloc((size_t) nmesh, sizeof(double));  /* wfc change */
    double *e = (double *) calloc((size_t) nmesh, sizeof(double));
    double *f = (double *) calloc((size_t) nmesh, sizeof(double));
    double *h = (double *) calloc((size_t) nmesh, sizeof(double));
    double *d = (double *) calloc((size_t) nmesh, sizeof(double));
    double sumdelu;                                                   /* sum of delus */
    int i;
    int j;
    double dellam;  /* change in the eigenvalue at end of each iteration */
    double lam;         /* eigen value */
    double integrand[nmesh];
    double integrand1[nmesh];
    double k;           /* constant to fix boundary conditions at large i */
    double *out1 = (double *) calloc((size_t) nmesh, sizeof(double));
    double *out2 = (double *) calloc((size_t) nmesh, sizeof(double));

/* initialize stuff. n.b. nmesh is set in globals.c  */
    FILE *io = fopen("../data_input", "r+w");

    get_ho();                 /* Get trial wave harmonic osc wave funtcion */
    dellam = 1;           /* change in egienvalue; always zero for BG. */
    sumdelu = 1.0;
    lam = -15.95;
    j = 0;
    pot[0] = get_pot(0);
    pot[1] = get_pot(1);
    a[0] = 1. + (mesh * mesh / 12.) * (lam - pot[1]);
    b[0] = -2. + (5.0 / 6.0) * pow(mesh, 2.) * (lam - pot[0]);
    c[0] = 1. - (pow(mesh, 2.) / 12.) * lam;
    d3[0] = (mesh*mesh/12.)*u[1]*dellam;
    u[0] = 0;
    q[0] = a[0] * u[1];
    f[0] = -(a[0]*r[1])/b[0];
    h[0] = -a[0]/b[0];
    e[0] = -(mesh*mesh*r[1]/12.)/b[0];

/* Basic loop dee--loop */
    while (fabs(sumdelu) > .001) {
        sumdelu = 0;
        for (i = 0; i < nmax; i++) {
            delu[i] = 0;
            if (j == 0) u[i] = r[i];  /*  first time fix the H.O. wfct for summation in Pauli excl. princ */
            pot[i] = get_pot(i);              /* output is in v[] v[rho] = <Y[jlst]| v[hamada] | Y[jlst]> */
        }

        j = j + 1;

        for (i = 1; i < nmax; i++) {
            a[i] = 1. + (pow(mesh, 2.) / 12.) * (lam - pot[i + 1]);
            b[i] = -2. + (5.0 / 6.0) * pow(mesh, 2.) * (lam - pot[i]);
            c[i] = 1. + (pow(mesh, 2.) / 12.) * (lam - pot[i - 1]);
            q[i] = (a[i] * u[i + 1] + b[i] * u[i] + c[i] * u[i - 1]);
            d3[i] = (mesh*mesh/12.0)*(u[i + 1] + 10.0*u[i] + u[i - 1]);
            d[i] = b[i]+ c[i]*h[i-1];
            f[i] = -(q[i]+ c[i]*f[i-1])/d[i];
            h[i] = -a[i]/d[i];
            e[i] = -(d3[i] + c[i]*e[i-1]);
            out1[i] = f[i]/e[i];

            k = pow(2.718, -mesh * mesh * (2.0 * nmax + 1) / 2);
            if (i == 1) {
                delu[i] = -(e[0]*dellam + f[0])/h[0];
                delu[i] =0.0;
            } else {
               delu[i] = (delu[i-1] - e[i]*dellam - f[i])/h[i];
            }
            u[i] = u[i] + delu[i];
            sumdelu = sumdelu + delu[i];
        }
        printf("sumdelu = %g\n", sumdelu);
        fflush(stdout);

        write_file(u, nmesh, io);
//        write_file(e, nmesh, io);
//        write_file(out1, nmesh, io);

        if (j > 1) {
            printf("no convergence after tries = %d\n", j);
            fflush(stdout);
            sumdelu = 0.0;
            dellam = 0.0;
        }
        dellam = -f[nmax-1]/e[nmax-1];
        dellam = 1.0;
        lam = lam + dellam;
        printf("lam = %g\n", lam);
        fflush(stdout);
        printf("dellam = %g\n", dellam);
        fflush(stdout);

        for (i = 0; i < nmax; i++) {
            integrand1[i] = u[i] * u[i];              /* normalize u */
        }
        get_integral(integrand1);

        for (i = 0; i < nmax; i++) {
            u[i] = u[i] / tmat;                       /* finish normalization */
        }

        for (i = 0; i < nmax; i++) {
            integrand[i] = r[i] * pot[i] * u[i];
        }
        get_integral(integrand);   /* recomp t-mat if new u */

    }
    return;
}