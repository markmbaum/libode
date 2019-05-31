#include "OdeEuler.hpp"

OdeEuler::OdeEuler (unsigned long long neq_) : OdeBaseRK (neq_) {
    //time stepping vectors
    k1 = new double[neq];
}

OdeEuler::~OdeEuler () {
    delete [] k1;
}

void OdeEuler::step (double dt) {

    unsigned long long i;

    //compute slope
    ode_funk(t, sol, k1);

    //compute solution
    for (i=0; i<neq; i++) sol[i] += dt*k1[i];

    //increase time
    t += dt;
    //update function evaluation counter
    neval += 1;
    //add a step to the counter
    nstep++;
}
