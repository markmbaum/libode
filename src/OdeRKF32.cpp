#include "OdeRKF32.hpp"

OdeRKF32::OdeRKF32 (unsigned long long neq_) : OdeBaseARK(neq_, 2) {

    //time stepping arrays
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    //tableau of coefficients
    c2 =   1.0; a21 = 1.0;
    c3 = 1.0/2; a31 = 1.0/4; a32 = 1.0/4;
                b1 =  1.0/2; b2  = 1.0/2; b3 =  0.0;
                d1 =  1.0/6; d2 =  1.0/6; d3  = 4.0/6;
}

//destructor
OdeRKF32::~OdeRKF32() {
    delete [] k1;
    delete [] k2;
    delete [] k3;
}

void OdeRKF32::step (double dt) {

    unsigned long long i;

    //first rk step
    ode_funk(t, sol, k1);
    //second step, first using k1 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);
    //third step, first using k2 and k1 to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i] + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //second order result
        sollo[i] = sol[i] + dt*(b1*k1[i] + b2*k2[i]);
        //third order result
        solhi[i] = sol[i] + dt*(d1*k1[i] + d2*k2[i] + d3*k3[i]);
    }

    //take the third order values as the new ones, keeping previous values
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increment the time
    t += dt;
    //update function evaluation counter
    neval += 3;
    //add a step to the counter
    nstep++;
}
