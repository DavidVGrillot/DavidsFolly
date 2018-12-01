/* *************************************  */
/* Subroutine: Given deuteron or quark and baryon type compute the reduced masses */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "globals.h"
#include "get_red_mass.h"
#include "get_integral.h"

int get_red_mass(){

/* out put is the reduced masses mu  in units MEV/c^2 */

double m1;                    /* mass of first particle */
double m2;                    /* mass of second particle */
double m3;                    /* mass of third particle */
double cq[3];                 /* charges of three quarks */

if(!strcmp(b_type , "deut") ){
/* deuterium nucleus */
    m1 = 938.272;          /* proton mass */
    m2 = 939.0/575.0;          /* neutron mass */
    m3 = 0;
    cq[0] = 1.;
    cq[1] = 0.;

}
else if (!strcmp(b_type , "prot")) {
/*  quarks uud */
    m1 = 2.2;
    m2 = 2.2;
    m3 = 4.7;
    /* set charge */
    cq[0] = 2./3.;
    cq[1] = 2./3.;
    cq[2] = -1./3.;
}
else if (!strcmp(b_type , "neu0")) {
/* DUD */
    m1 = 4.7;
    m2 = 2.2;
    m3 = 4.7;
    /* set charge */
    cq[0] = -1./3.;
    cq[1] = 2./3.;
    cq[2] = -1./3.;

}
else if (!strcmp(b_type , "sig+")) {
/*  UUS */

    m1 = 2.2;
    m2 = 2.2;
    m3 = 96;
    /* set charge */
    cq[0] = 2./3.;
    cq[1] = 2./3.;
    cq[2] = -1./3.;
}
else if (!strcmp(b_type , "sig-")) {
/* DDS */
    m1 = 4.7;
    m2 = 4.7;
    m3 = 96.0;
    /* set charge */
   cq[0] = -1./3.;
   cq[1] = -1./3.;
   cq[2] = -1./3.;

}
else if (!strcmp(b_type , "sig0")) {
/* UDS */
   m1 = 2.2;
   m2 = 4.7;
   m3 = 96.0;
   /* set charge */
   cq[0] = 2./3.;
   cq[1] = -1./3.;
   cq[2] = -1./3.;

}
else if (!strcmp(b_type , "xii-")) {
/*  SDS */
   m1 = 96.0;
   m2 = 4.7;
   m3 = 96.0;
   /* set charge */
   cq[0] = -1./3.;
   cq[1] = -1./3.;
   cq[2] = -1./3.;

}
else if (!strcmp(b_type , "xii0")) {
/* SUS */
   m1 = 96.0;
   m2 = 2.2;
   m3 = 96.0;
   /* set charge */
   cq[0] = -1./3.;
   cq[1] = 2./3.;
   cq[2] = -1./3.;

}
else if (!strcmp(b_type ,  "del-")) {
/* DDD */
    m1 = 4.7;
    m2 = 4.7;
    m3 = 4.7;
    /* set charge */
    cq[0] = -1./3.;
    cq[1] = -1./3.;
    cq[2] = -1./3.;

}
else if(!strcmp(b_type ,  "del0")) {
/* DUD hmmmmmmmmmmmmm same as neutron?? */
    m1 = 4.7;
    m2 = 2.2;
    m3 = 4.7;
    /* set charge */
    cq[0] = -1./3.;
    cq[1] = 2./3.;
    cq[2] = -1./3.;

}
else if(!strcmp(b_type , "del+")) {
       /* same as proton dud */
    m1 = 7.7;
    m2 = 7.7;
    m3 = 2.2;
    /* set charge */
    cq[0] = -1./3.;
    cq[1] = 2./3.;
    cq[2] = -1./3.;

}
else if(!strcmp(b_type , "del++")) {
    m1 = 2.2;
    m2 = 2.2;
    m3 = 2.2;
    /* set charge */
    cq[0] = 2./3.;
    cq[1] = 2./3.;
    cq[2] = 2./3.;

}
else {
    return 1;
}
mu[0] = m1 * m2 / (m1 + m2);
mu[1] = m1 * m3 / (m1 + m3);
mu[2] = m2 * m3 / (m2 + m3);
return 0;
}