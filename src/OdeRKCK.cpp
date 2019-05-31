#include "OdeRKCK.hpp"

OdeRKCK::OdeRKCK (unsigned long long neq_) : OdeBaseARK(neq_, 4) {

    //time stepping arrays
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    k5 = new double[neq];
    k6 = new double[neq];
    //coefficents of tableau
    c2 =  1.0/5; a21 =        1.0/5;
    c3 = 3.0/10; a31 =       3.0/40; a32 =    9.0/40;
    c4 =  3.0/5; a41 =       3.0/10; a42 =   -9.0/10; a43 =        6.0/5;
    c5 =    1.0; a51 =     -11.0/54; a52 =     5.0/2; a53 =     -70.0/27; a54 = 35.0/27;
    c6 =  7.0/8; a61 = 1631.0/55296; a62 = 175.0/512; a63 =  575.0/13824; a64 = 44275.0/110592; a65 = 253.0/4096;
                  b1 =     37.0/378; b2  =       0.0; b3 =     250.0/621; b4  =      125.0/594; b5 =         0.0; b6 = 512.0/1771;
                  d1 = 2825.0/27648; d2 =        0.0; d3 = 18575.0/48384; d4 =   13525.0/55296; d5 = 277.0/14336; d6 = 1.0/4;
                  e1 =      19.0/54; e2 =        0.0; e3 =      -10.0/27; e4 =         55.0/54;
                  f1 =       -3.0/2; f2 =      5.0/2;
                  g1 =          1.0;
}

//destructor
OdeRKCK::~OdeRKCK() {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
}

void OdeRKCK::step (double dt) {

    unsigned long long i;

    //K1
    ode_funk(t, sol, k1);

    //order 1 solution
    //for (i=0; i<neq; i++) sollo[i] = sol[i] + dt*g1*k1[i];

    //K2
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //order 2 solution
    //for (i=0; i<neq; i++) solhi[i] = sol[i] + dt*(f1*k1[i] + f2*k2[i]);

    //K3
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i]
                                                  + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //K4
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a41*k1[i]
                                                  + a42*k2[i]
                                                  + a43*k3[i]);
    ode_funk(t + dt*c4, soltemp, k4);

    //order 3 solution
    //for (i=0; i<neq; i++) {
    //    sollo[i] = solhi[i];
    //    solhi[i] = sol[i] + dt*(e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i]);
    //}

    //K5
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a51*k1[i]
                                                  + a52*k2[i]
                                                  + a53*k3[i]
                                                  + a54*k4[i]);
    ode_funk(t + dt*c5, soltemp, k5);

    //K6
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a61*k1[i]
                                                  + a62*k2[i]
                                                  + a63*k3[i]
                                                  + a64*k4[i]
                                                  + a65*k5[i]);
    ode_funk(t + dt*c6, soltemp, k6);

    //order 4 and 5 solutions
    for (i=0; i<neq; i++) {
        sollo[i] = sol[i] + dt*(d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i]);
        solhi[i] = sol[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b6*k6[i]);
    }

    //use the 5th order solution
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increment the time
    t += dt;
    //update function evaluation counter
    neval += 6;
    //add a step to the counter
    nstep++;
}
