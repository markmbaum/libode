#include "OdeRKF21.hpp"

OdeRKF21::OdeRKF21 (unsigned long long neq_) : OdeBaseARK(neq_, 1) {

    //time stepping vectors
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    //turn on fsal
    fsal = true;
    //pointer to final k values for fsal step rejection
    klast = k3;
    //tableau of coefficients
    c2 = 1.0/2; a21 =   1.0/2;
    c3 =   1.0; a31 = 1.0/256; a32 = 255.0/256;
                b1  = 1.0/256; b2  = 255.0/256;
                d1  = 1.0/512; d2  = 255.0/256; d3 = 1.0/512;
}

//destructor
OdeRKF21::~OdeRKF21() {
    delete [] k1;
    delete [] k2;
    delete [] k3;
}

void OdeRKF21::step (double dt) {

    unsigned long long i;

    //first rk step, FSAL
    for (i=0; i<neq; i++) k1[i] = klast[i];
    //second step, first using k1 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);
    //third step, first using k2 and k1 to get temp y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i] + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //fourth order result
        sollo[i] = sol[i] + dt*(b1*k1[i] + b2*k2[i]);
        //third order result
        solhi[i] = sol[i] + dt*(d1*k1[i] + d2*k2[i] + d3*k3[i]);
    }

    //take the 2nd order solution
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increment the time
    t += dt;
    //update function evaluation counter
    neval += 2;
    //add a step to the counter
    nstep++;
}
