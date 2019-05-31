#include "OdeRK43.hpp"

OdeRK43::OdeRK43 (unsigned long long neq_) : OdeBaseARK(neq_, 3) {

    //time stepping arrays
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    k5 = new double[neq];
    //turn fsal on
    fsal = true;
    //pointer to final k values for fsal step rejection
    klast = k5;
    //tableau of coefficents
    c2 = 1.0/3; a21 = 1.0/3;
    c3 = 2.0/3; a31 = -1./3; a32 =     1;
    c4 =   1.0; a41 =     1; a42 =    -1; a43 = 1;
    c5 =   1.0; a51 = 1.0/8; a52 = 3.0/8; a53 = 3.0/8; a54 = 1.0/8;
                b1 =  1.0/8; b2  = 3.0/8; b3 =  3.0/8; b4  = 1.0/8;
                d1 =  1./12; d2 =   1./2; d3  = 1.0/4; d4 =      0; d5  = 1./6;
}

//destructor
OdeRK43::~OdeRK43 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
}

void OdeRK43::step (double dt) {

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

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //fourth order result
        solhi[i] = sol[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]);
        //third order result
        sollo[i] = sol[i] + dt*(d1*k1[i] + d2*k2[i] + d3*k3[i] + d5*k5[i]);
    }

    //take the 4th order solution
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increment the time
    t += dt;
    //update function evaluation counter
    neval += 4;
    //add a step to the counter
    nstep++;
}
