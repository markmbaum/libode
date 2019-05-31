#include "OdeRK4.hpp"

OdeRK4::OdeRK4 (unsigned long long neq_) : OdeBaseRK (neq_) {

    //time stepping vectors
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    //tableau of coefficients
    c2 = 1.0/2; a21 = 1.0/2;
    c3 = 1.0/2;              a32 = 1.0/2;
    c4 =   1.0;                           a43 =  1.0;
                 b1 = 1.0/6; b2  = 1.0/3; b3 = 1.0/3; b4 = 1.0/6;
}

OdeRK4::~OdeRK4 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
}

void OdeRK4::step (double dt) {

    unsigned long long i;

    //k1
    ode_funk(t, sol, k1);
    //k2
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);
    //k3
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a32*k2[i];
    ode_funk(t + dt*c3, soltemp, k3);
    //k4
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a43*k3[i];
    ode_funk(t + dt*c4, soltemp, k4);

    //store solution
    for (i=0; i<neq; i++) sol[i] += dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]);

    //increase time
    t += dt;
    //update function evaluation counter
    neval += 4;
    //add a step to the counter
    nstep++;
}
