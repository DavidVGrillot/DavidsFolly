/*  subroutine to integrate the input integrand and return value in tmat (global) */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "globals.h"
#include "get_integral.h"

void get_integral(double integrand[])    /* function to output the integral of integrand typicall r*v*u */

{
/* input integrand to be integrated...output is put into global tmat */
    double lhs;                     /* left hand sum */
    double rhs;                     /* right hand sum */
    int i;
    int max;



/* integration routine computes (LHS +RHS)/2 return  tmat */
    lhs = 0.0;
    rhs = 0.0;
    max = nmax - 1;
    for(i=0; i<max; i++)  {

        lhs = lhs + integrand[i];
        rhs = rhs +integrand[i + 1];
    }
    lhs = lhs*mesh;
    rhs = rhs*mesh;
    tmat = (lhs + rhs)/2.;

    return;
}
