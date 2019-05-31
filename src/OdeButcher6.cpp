#include "OdeButcher6.hpp"

OdeButcher6::OdeButcher6 (unsigned long long neq_) : OdeBaseRK (neq_) {

    //time stepping vectors
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    k5 = new double[neq];
    k6 = new double[neq];
    k7 = new double[neq];
    //tableau of coefficients
    c2 = 1.0/2; a21 =     1.0/2;
    c3 = 2.0/3; a31 =     2.0/9; a32 =    4.0/9;
    c4 = 1.0/3; a41 =    7.0/36; a42 =    2.0/9; a43 =  -1.0/12;
    c5 = 5.0/6; a51 = -35.0/144; a52 = -55.0/36; a53 =  35.0/48; a54 =    15.0/8;
    c6 = 1.0/6; a61 =  -1.0/360; a62 = -11.0/36; a63 =   -1.0/8; a64 =     1.0/2; a65 =   1.0/10;
    c7 =   1.0; a71 = -41.0/260; a72 =  22.0/13; a73 = 43.0/156; a74 = -118.0/39; a75 = 32.0/195; a76 = 80.0/39;
                b1 =   13.0/200; b2 =       0.0; b3 =   11.0/40; b4 =    11.0/40; b5 =    4.0/25; b6 =    4.0/25; b7 = 13.0/200;
}

OdeButcher6::~OdeButcher6 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
    delete [] k7;
}

void OdeButcher6::step (double dt) {

    unsigned long long i;

    //k1
    ode_funk(t, sol, k1);

    //k2
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //k3
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i]
                                                  + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //k4
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a41*k1[i]
                                                  + a42*k2[i]
                                                  + a43*k3[i]);
    ode_funk(t + dt*c4, soltemp, k4);

    //k5
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a51*k1[i]
                                                  + a52*k2[i]
                                                  + a53*k3[i]
                                                  + a54*k4[i]);
    ode_funk(t + dt*c5, soltemp, k5);

    //k6
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a61*k1[i]
                                                  + a62*k2[i]
                                                  + a63*k3[i]
                                                  + a64*k4[i]
                                                  + a65*k5[i]);
    ode_funk(t + dt*c6, soltemp, k6);

    //k7
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a71*k1[i]
                                                  + a72*k2[i]
                                                  + a73*k3[i]
                                                  + a74*k4[i]
                                                  + a75*k5[i]
                                                  + a76*k6[i]);
    ode_funk(t + dt*c7, soltemp, k7);

    //store solution
    for (i=0; i<neq; i++) sol[i] += dt*(b1*k1[i]
                                      + b3*k3[i]
                                      + b4*k4[i]
                                      + b5*k5[i]
                                      + b6*k6[i]
                                      + b7*k7[i]);

    //increase time
    t += dt;
    //update function evaluation counter
    neval += 7;
    //add a step to the counter
    nstep++;
}
