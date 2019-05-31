#include "OdeTrapz.hpp"

OdeTrapz::OdeTrapz (unsigned long long neq_) : OdeBaseARK (neq_, 1) {

    //time stepping vectors
    k1 = new double[neq];
    k2 = new double[neq];
    //tableau of coefficients
    c2 = 1.0; a21 =  1.0;
              b1  = 1.0/2; b2 = 1.0/2;
}

//destructor
OdeTrapz::~OdeTrapz() {
    delete [] k1;
    delete [] k2;
}

void OdeTrapz::step (double dt) {

    unsigned long long i;

    //k1
    ode_funk(t, sol, k1);
    //second step, first using k1 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //second order result
        solhi[i] = sol[i] + dt*(b1*k1[i] + b2*k2[i]);
        //first order result (euler)
        sollo[i] = sol[i] + dt*k1[i];
    }

    //take the third order values as the new ones, keeping previous values
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increase time
    t += dt;
    //update function evaluation counter
    neval += 2;
    //add a step to the counter
    nstep++;
}
