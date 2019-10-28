//! \file ode_euler.cc

#include "ode_euler.h"

OdeEuler::OdeEuler (unsigned long neq) :
    OdeBase(neq, false),
    OdeRK (neq, 1) {

    method_ = "Euler";
}

void OdeEuler::step_ (double dt) {

    //index
    unsigned long i;
    //compute slope
    ode_fun_(sol_, k_[0]);
    //compute solution
    for (i=0; i<neq_; i++) sol_[i] += dt*k_[0][i];
}
