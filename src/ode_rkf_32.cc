#include "ode_rkf_32.h"

OdeRKF32::OdeRKF32 (unsigned long neq) :
    OdeEmbedded (neq, false, 2),
    OdeRK (neq, 3),
    OdeERK (neq) {

    method_ = "RKF32";
    //tableau of coefficients
    c2 =   1.0; a21 = 1.0;
    c3 = 1.0/2; a31 = 1.0/4; a32 = 1.0/4;
                b1 =  1.0/2; b2  = 1.0/2;
                d1 =  1.0/6; d2 =  1.0/6; d3  = 4.0/6;
}

void OdeRKF32::step_ (double dt) {

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
    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + dt*(b1*k_[0][i] + b2*k_[1][i]);
        sol_[i] = sol_[i] + dt*(d1*k_[0][i] + d2*k_[1][i] + d3*k_[2][i]);
    }
}
