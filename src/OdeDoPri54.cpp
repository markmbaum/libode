#include "OdeDoPri54.hpp"

OdeDoPri54::OdeDoPri54 (unsigned long long neq_) : OdeBaseARK (neq_, 4) {

    //time stepping arrays
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    k5 = new double[neq];
    k6 = new double[neq];
    k7 = new double[neq];
    //turn fsal on
    fsal = true;
    //pointer to final k values for fsal step rejection
    klast = k7;
    //tableau of coefficents
    c2 =  1.0/5; a21 =        1.0/5;
    c3 = 3.0/10; a31 =       3.0/40; a32 =        9.0/40;
    c4 =  4.0/5; a41 =      44.0/45; a42 =      -56.0/15; a43 = 32.0/9;
    c5 =  8.0/9; a51 = 19372.0/6561; a52 = -25360.0/2187; a53 = 64448.0/6561; a54 = -212.0/729;
    c6 =    1.0; a61 =  9017.0/3168; a62 =     -355.0/33; a63 = 46732.0/5247; a64 =   49.0/176; a65 =   -5103.0/18656;
    c7 =    1.0; a71 =     35.0/384; a72 =           0.0; a73 =   500.0/1113; a74 =  125.0/192; a75 =    -2187.0/6784; a76 =    11.0/84;
                  b1 =     35.0/384; b2  =           0.0; b3  =   500.0/1113; b4  =  125.0/192; b5  =    -2187.0/6784;  b6 =    11.0/84;
                  d1 = 5179.0/57600; d2  =           0.0; d3  = 7571.0/16695; d4 =   393.0/640; d5  = -92097.0/339200;  d6 = 187.0/2100; d7 = 1.0/40;
}

OdeDoPri54::~OdeDoPri54 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
    delete [] k7;
}

void OdeDoPri54::step (double dt) {

    unsigned long long i;

    //first rk step, FSAL
    for (i=0; i<neq; i++) k1[i] = klast[i];

    //second step, first using k1 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //third step, first using k2 and k1 to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i]
                                                  + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //fourth step, first using previous k values to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a41*k1[i]
                                                  + a42*k2[i]
                                                  + a43*k3[i]);
    ode_funk(t + dt*c4, soltemp, k4);

    //fifth step, first using previous k values to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a51*k1[i]
                                                  + a52*k2[i]
                                                  + a53*k3[i]
                                                  + a54*k4[i]);
    ode_funk(t + dt*c5, soltemp, k5);

    //sixth step, first using previous k values to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a61*k1[i]
                                                  + a62*k2[i]
                                                  + a63*k3[i]
                                                  + a64*k4[i]
                                                  + a65*k5[i]);
    ode_funk(t + dt*c6, soltemp, k6);

    //seventh step, first using previous k values to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a71*k1[i]
                                                  + a73*k3[i]
                                                  + a74*k4[i]
                                                  + a75*k5[i]
                                                  + a76*k6[i]);
    ode_funk(t + dt*c7, soltemp, k7);

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //fourth order result
        solhi[i] = sol[i] + dt*(b1*k1[i]
                              + b3*k3[i]
                              + b4*k4[i]
                              + b5*k5[i]
                              + b6*k6[i]);
        //third order result
        sollo[i] = sol[i] + dt*(d1*k1[i]
                              + d3*k3[i]
                              + d4*k4[i]
                              + d5*k5[i]
                              + d6*k6[i]
                              + d7*k7[i]);
    }

    //take the 5th order solution
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
