//! \file ode_rk_4.cc

#include "ode_rk_4.h"

OdeRK4::OdeRK4 (unsigned long neq) :
    OdeAdaptive (neq, false),
    OdeRK (neq, 4),
    OdeERK (neq) {

    method_ = "RK4";
    //tableau of coefficients
    c2 = 1.0/2; a21 =    c2;
    c3 = 1.0/2;              a32 = 1.0/2;
    c4 =   1.0;                           a43 =  1.0;
                 b1 = 1.0/6; b2  = 1.0/3; b3 = 1.0/3; b4 = 1.0/6;
}

void OdeRK4::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    ode_fun_(sol_, k_[0]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a21*k_[0][i];
    ode_fun_(soltemp_, k_[1]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a32*k_[1][i];
    ode_fun_(soltemp_, k_[2]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a43*k_[2][i];
    ode_fun_(soltemp_, k_[3]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) sol_[i] += dt*(b1*k_[0][i]
                                        + b2*k_[1][i]
                                        + b3*k_[2][i]
                                        + b4*k_[3][i]);
}
