//! \file ode_backward_euler.cc

#include "ode_backward_euler.h"

void NewtonBackwardEuler::f_Newton (double *x, double *y) {

    (void)x; //supress unused variable warning
    unsigned long i;
    double dt = *dt_; //get time step from solver object

    //evaluate temporary solution based on current k values
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*k_[0][i];
    //compute k1 values
    fun(soltemp_, ftemp_);
    //evaluate the Newton function's k1 component
    for (i=0; i<neq_; i++) y[i] = k_[0][i] - ftemp_[i];
}

void NewtonBackwardEuler::J_Newton (double *x, double **J) {

    (void)x; //supress unused variable warning
    unsigned long i,j;
    double dt = *dt_; //get time step from solver object

    if ( get_modified() ) jac(sol_, Jac_);

    if ( !get_modified() ) {
        //evaluate temporary solution based on current k values for k1
        for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*k_[0][i];
        //evaluate J at soltemp
        jac(soltemp_, Jac_);
    }
    //evaluate the Newton jacobian
    for (i=0; i<neq_; i++) {
        for (j=0; j<neq_; j++) {
            J[i][j] = -dt*Jac_[i][j];
        }
        J[i][i] += 1.0; //subtract from the identity matrix
    }
}

//------------------------------------------------------------------------------

OdeBackwardEuler::OdeBackwardEuler (unsigned long neq) :
    OdeAdaptive (neq, true),
    OdeIRK (neq, 1) {

    method_ = "BackwardEuler";

    newton_ = new NewtonBackwardEuler(neq, 1, this);
    newton_->set_modified(true);
}

OdeBackwardEuler::~OdeBackwardEuler () {}

void OdeBackwardEuler::step_ (double dt) {

    unsigned long i;
    //set k values to zero
    for (i=0; i<neq_; i++) k_[0][i] = 0.0;
    //solve for slopes
    newton_->solve_Newton(k_[0]);
    //compute new solution
    for (i=0; i<neq_; i++) sol_[i] += dt*k_[0][i];
}
