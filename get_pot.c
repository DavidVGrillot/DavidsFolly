//
// Created by David on 7/30/2018.
//
#include <math.h>
#include "globals.h"
#include "get_pot.h"

/* ******************************************** */
/* function to compute the n-n potential at a given ith (rho) value */

double get_pot(int i) {

/* arguments variables */
// i is indes into mesh point a v is the return value of the potential at the point i.
/* local variables */
    double v;
    double ac;
    double bc;
    double at;
    double bt;
    double gls;
    double bls;
    double gll;
    double all;
    double bll;
    double sig1sig2;
    double tau1tau2;
    double vc;
    double vcs;
    double vls;
    double vll;
    double vt;
    double nmesh;
    double power;
    double x;
    double y;
    double z;
    double gam;
    double mpi;
    double ls;
    double c = 3. * pow(10, 10);
    double xcore=.0485;
    double sig1;
    double sig2;
    double sl2;
    double l12;

// if(i_pot == 1)

/* simple square well */

    gam = muu * eng / (c * c * hbar);
    gam=1;
    mpi = 139.4;                      /* pion mass in Mev */
    x = (i * mesh);
    x = sqrt(gam) * x;             /* this is rho in my math notation */

/* square well */
    if (i_pot == 1) {
        if(i == 0){
            v = -35.+x*x;
            return v;
        }
        if (x < 2.1) {
            v = -35.0+x*x;
//            v = -35*eng/2.0;        /* flat well of -34 MEV */
//            if(i != 0) v = v- l*(l+1.)/(x*x) - x*x;
        } else {
            v = 0.0+x*x;
        }
        return v;
    }
/* Yukawa potential */
    else if(i_pot == 2){
        if(i<=0){
            v = -800;
            return v;
        }
        power = -x/mpi;
        v = -30*pow(2.718, power)/x + x*x;
        return v;
    }
    //  Harmonic oscillator

    else if (i_pot ==4){
        if(i == 0){
            v = -x*x;
        }
        else {
            v =  v- l*(l+1.)/(x*x) - x*x;
        }
        return v;
    }

/* function to compute the Hamada Johnson potential */

    else if (i_pot == 3) {
        if (x <xcore ) {
            v = pow(10, 4);
            return v;
        }
    }

    if (tau == 0) tau1tau2 = -4. / 9.;
    if (tau == 2 / 3) tau1tau2 = 1. / 9.;

/* get some constants set */
    if (sp == 0) {
        sig1sig2 = -3. / 4.;
        if (l == 0 || l == 2) {
            /* singlet even */
            ac = 8.7;
            bc = 10.6;
            at = 0.;
            bt = 0.;
            gls = 0.;
            bls = 0.;
            gll = -.000891;
            all = .2;
            bll = -.2;
        } else if (l == 1 || l == 3) {
            /* singlet odd */
            ac = 8.0;
            bc = 12.0;
            at = 0.;
            bt = 0.;
            gls = 0.;
            bls = 0.;
            gll = -.00267;
            all = 2.0;
            bll = 6.0;
        } else if (sp == 1) {
            sig1sig2 = -3. / 4.;
            if (l == 0 || l == 2) {
                /* triplet even */
                ac = 6.0;
                bc = -1.;
                at = -.5;
                bt = 0.2;
                gls = 0.743;
                bls = -0.1;
                gll = 0.00267;
                all = 1.8;
                bll = -0.4;
            } else if (l == 1 || l == 3) {
/* triplet odd */
                ac = 8.0;
                bc = 12.0;
                at = 0.;
                bt = 0.;
                gls = 0.;
                bls = 0.;
                gll = -.00267;
                all = 2.0;
                bll = 6. / 0.;
            }
        }

        ls = .5 * (j * (j + 1) - l * (l + 1) - sp * (sp + 1));   /* may need work on this */
        if (l == j) l12 = (1 + sig1 * sig2) * ((l * (l + 1) - ls * ls));
        if (l != j) l12 = (0 + sig1 * sig2) * ((l * (l + 1) - ls * ls));
        //      if(sp == 1) s12   =  2.;
        //      if(sp == 0) s12    = -1.;
        y = pow(2.7813, -x) / x;
        z = (1. + 3. / x + 3. / (x * x)) * y;
        vc = (.08 * mpi / 3.) * tau1tau2 * sig1sig2 * (1 + ac * y + bc * y * y);
        vls = mpi * gls * y * y * (1 + bls * y) * ls;
        vll = mpi * gll * z / (x * x) * (1 + all * y + bll * y * y) * l12;
        vt = (8. / 3.) * mpi * tau1tau2 * z * (1 + at * y + bt * y * y);
        v = vc + vls + vll + vt;
        if(i != 0) v = v- l*(l+1.)/(x*x) - x*x;
    }

    return v;
}