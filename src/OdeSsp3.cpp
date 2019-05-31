
#include "OdeSsp3.hpp"

OdeSsp3::OdeSsp3 (unsigned long long neq_) : OdeBaseRK (neq_) {
    //time stepping vectors
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    //tableau of coefficients
    c2 =   1.0; a21 =   1.0;
    c3 = 1.0/2; a31 = 1.0/4; a32 = 1.0/4;
                 b1 = 1.0/6; b2 =  1.0/6; b3 = 2.0/3;

}

OdeSsp3::~OdeSsp3 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
}

void OdeSsp3::step (double dt) {

    unsigned long long i;

    //k1
    ode_funk(t, sol, k1);

    //second step, first using k1 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //third step, first using k1 and k2 to get temporary y values
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i] + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //store solution
    for (i=0; i<neq; i++) sol[i] += dt*(b1*k1[i] + b2*k2[i] + b3*k3[i]);

    //increase time
    t += dt;
    //update function evaluation counter
    neval += 3;
    //add a step to the counter
    nstep++;
}
