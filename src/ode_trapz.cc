//! \file ode_trapz.cc

#include "ode_trapz.h"

OdeTrapz::OdeTrapz (unsigned long neq) :
    OdeBase (neq, false),
    OdeRK (neq, 2),
    OdeERK (neq) {

    method_ = "Trapz";
    //tableau of coefficients
    c2 = 1.0; a21 =  1.0;
              b1  = 1.0/2; b2 = 1.0/2;
}

void OdeTrapz::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    ode_fun_(sol_, k_[0]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a21*k_[0][i];
    ode_fun_(soltemp_, k_[1]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) sol_[i] = sol_[i] + dt*(b1*k_[0][i] + b2*k_[1][i]);
}
