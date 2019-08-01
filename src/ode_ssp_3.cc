#include "ode_ssp_3.h"

OdeSsp3::OdeSsp3 (unsigned long neq) :
    OdeBase (neq, false),
    OdeRK (neq, 3),
    OdeERK (neq) {

    method_ = "Ssp3";
    //tableau of coefficients
    c2 =   1.0; a21 =    c2;
    c3 = 1.0/2; a31 = 1.0/4; a32 = 1.0/4;
                 b1 = 1.0/6; b2 =  1.0/6; b3 = 2.0/3;
}

void OdeSsp3::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    ode_fun_(sol_, k_[0]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a21*k_[0][i];
    ode_fun_(soltemp_, k_[1]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a31*k_[0][i] + a32*k_[1][i]);
    ode_fun_(soltemp_, k_[2]);

    //------------------------------------------------------------------
    //store solution
    for (i=0; i<neq_; i++) sol_[i] += dt*(b1*k_[0][i] + b2*k_[1][i] + b3*k_[2][i]);
}
