//! \file ode_sdirk_43.cc

#include "ode_sdirk_43.h"

namespace ode {

void NewtonSDIRK43::f_Newton (double *x, double *y) {

    (void)x;
    unsigned long i;
    int m;
    double dt = *dt_;

    //evaluate solution based on current k values
    for (i=0; i<neq_; i++) {
        //diagonal term
        soltemp_[i] = sol_[i] + dt*gam*k_[ik_][i];
        //terms below the diagonal
        for (m=0; m<ik_; m++) soltemp_[i] += dt*a[ik_][m]*k_[m][i];
    }
    fun(soltemp_, ftemp_);
    //evaluate the Newton function
    for (i=0; i<neq_; i++) y[i] = k_[ik_][i] - ftemp_[i];
}

void NewtonSDIRK43::J_Newton (double *x, double **J) {

    (void)x;
    unsigned long i,j;
    int m;
    double dt = *dt_;

    //evaluate Jacobian with current k values if needed
    if ( !get_modified() ) {
        for (i=0; i<neq_; i++) {
            //diagonal term
            soltemp_[i] = sol_[i] + dt*gam*k_[ik_][i];
            //terms below the diagonal
            for (m=0; m<ik_; m++) soltemp_[i] += dt*a[ik_][m]*k_[m][i];
        }
        jac(soltemp_, Jac_);
    }
    //evaluate the Newton system's Jacobian
    for (i=0; i<neq_; i++) {
        for (j=0; j<neq_; j++) {
            J[i][j] = -dt*gam*Jac_[i][j];
        }
        J[i][i] += 1.0;
    }
}

//------------------------------------------------------------------------------

OdeSDIRK43::OdeSDIRK43 (unsigned long neq) :
    OdeEmbedded (neq, true, 3),
    OdeRK (neq, 5) {

    method_ = "SDIRK43";

    int nk = 5;
    a = new double*[nk];
    for (int i=0; i<nk; i++) a[i] = new double[nk];
    b = new double[nk];
    d = new double[nk];

    gam = 1.0/4;

    a[1][0] =      1.0/2;
    a[2][0] =    17.0/50; a[2][1] =     -1.0/25;
    a[3][0] = 371.0/1360; a[3][1] = -137.0/2720; a[3][2] = 15.0/544;
    a[4][0] =    25.0/24; a[4][1] =    -49.0/48; a[4][2] = 125.0/16; a[4][3] = -85.0/12;

    b[0] =       a[4][0]; b[1] =        a[4][1]; b[2] =     a[4][2]; b[3] =     a[4][3]; b[4] = gam;
    d[0] =       59.0/48; d[1] =       -17.0/96; d[2] =    225.0/32; d[3] =    -85.0/12;

    newton_ = new NewtonSDIRK43(neq, this);
    newton_->set_modified(true);
}

OdeSDIRK43::~OdeSDIRK43 () {
    for (int i=0; i<nk_; i++) delete [] a[i];
    delete [] a;
    delete [] b;
    delete [] d;
}

void OdeSDIRK43::step_ (double dt) {

    unsigned long i;
    //compute the ode system's Jacobian if using modified Newton iterations
    if ( newton_->get_modified() ) ode_jac_(sol_, Jac_);
    //solve for first k vector
    for (i=0; i<neq_; i++) k_[0][i] = 0.0;
    newton_->set_ignore_JLU(false);
    newton_->set_ik(0);
    newton_->solve_Newton(k_[0]);
    //solve for subsequent stage k vectors
    newton_->set_ignore_JLU(true);
    for (int j=1; j<nk_; j++) {
        //use previous stage values as initial guess
        for (i=0; i<neq_; i++) k_[j][i] = k_[j-1][i];
        //solve
        newton_->set_ik(j);
        newton_->solve_Newton(k_[j]);
    }
    //compute solution
    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + dt*(d[0]*k_[0][i]
                                 + d[1]*k_[1][i]
                                 + d[2]*k_[2][i]
                                 + d[3]*k_[3][i]);
        sol_[i] += dt*(b[0]*k_[0][i]
                     + b[1]*k_[1][i]
                     + b[2]*k_[2][i]
                     + b[3]*k_[3][i]
                     + b[4]*k_[4][i]);
    }
}

} // namespace ode
