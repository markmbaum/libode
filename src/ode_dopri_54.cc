#include "ode_dopri_54.h"

OdeDoPri54::OdeDoPri54 (unsigned long neq) :
    OdeEmbedded (neq, false, 4),
    OdeRK (neq, 7),
    OdeERK (neq) {

    method_ = "DoPri54";
    //tableau of coefficents
    c2 =  1.0/5; a21 =        1.0/5;
    c3 = 3.0/10; a31 =       3.0/40; a32 =        9.0/40;
    c4 =  4.0/5; a41 =      44.0/45; a42 =      -56.0/15; a43 = 32.0/9;
    c5 =  8.0/9; a51 = 19372.0/6561; a52 = -25360.0/2187; a53 = 64448.0/6561; a54 = -212.0/729;
    c6 =    1.0; a61 =  9017.0/3168; a62 =     -355.0/33; a63 = 46732.0/5247; a64 =   49.0/176; a65 =   -5103.0/18656;
    c7 =    1.0; a71 =     35.0/384; a72 =           0.0; a73 =   500.0/1113; a74 =  125.0/192; a75 =    -2187.0/6784; a76 =    11.0/84;
                  b1 =     35.0/384; b2  =           0.0; b3  =   500.0/1113; b4  =  125.0/192; b5  =    -2187.0/6784;  b6 =    11.0/84;
                  d1 = 5179.0/57600; d2  =           0.0; d3  = 7571.0/16695; d4 =   393.0/640; d5  = -92097.0/339200;  d6 = 187.0/2100; d7 = 1.0/40;
}

void OdeDoPri54::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    ode_fun_(sol_, k_[0]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a21*k_[0][i];
    ode_fun_(soltemp_, k_[1]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a31*k_[0][i]
                                                     + a32*k_[1][i]);
    ode_fun_(soltemp_, k_[2]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a41*k_[0][i]
                                                     + a42*k_[1][i]
                                                     + a43*k_[2][i]);
    ode_fun_(soltemp_, k_[3]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a51*k_[0][i]
                                                     + a52*k_[1][i]
                                                     + a53*k_[2][i]
                                                     + a54*k_[3][i]);
    ode_fun_(soltemp_, k_[4]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a61*k_[0][i]
                                                     + a62*k_[1][i]
                                                     + a63*k_[2][i]
                                                     + a64*k_[3][i]
                                                     + a65*k_[4][i]);
    ode_fun_(soltemp_, k_[5]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a71*k_[0][i]
                                                     + a72*k_[1][i]
                                                     + a73*k_[2][i]
                                                     + a74*k_[3][i]
                                                     + a75*k_[4][i]
                                                     + a76*k_[5][i]);
    ode_fun_(soltemp_, k_[6]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + dt*(d1*k_[0][i]
                                 + d2*k_[1][i]
                                 + d3*k_[2][i]
                                 + d4*k_[3][i]
                                 + d5*k_[4][i]
                                 + d6*k_[5][i]
                                 + d7*k_[6][i]);
        sol_[i] = sol_[i] + dt*(b1*k_[0][i]
                              + b2*k_[1][i]
                              + b3*k_[2][i]
                              + b4*k_[3][i]
                              + b5*k_[4][i]
                              + b6*k_[5][i]);
    }
}