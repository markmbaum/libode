//! \file ode_radau_iia_5.cc

#include "ode_radau_iia_5.h"

void NewtonRadauIIA5::f_Newton (double *x, double *y) {

    (void)x; //supress unused variable warning
    unsigned long i;
    int j,n;
    double dt = *dt_; //get time step from solver object

    for (n=0; n<nk_; n++) {
        //evaluate temporary solution based on current k values
        for (i=0; i<neq_; i++) {
            soltemp_[i] = sol_[i];
            for (j=0; j<nk_; j++) soltemp_[i] += dt*a[n][j]*k_[j][i];
        }
        //compute k1 values
        fun(soltemp_, ftemp_);
        //evaluate the Newton function's k1 component
        for (i=0; i<neq_; i++) y[i+n*neq_] = k_[n][i] - ftemp_[i];
    }
}

void NewtonRadauIIA5::J_Newton (double *x, double **J) {

    (void)x; //supress unused variable warning
    unsigned long i,j;
    int m,n;
    double dt = *dt_; //get time step from solver object

    if ( get_modified() ) jac(sol_, Jac_);

    for (n=0; n<nk_; n++) {
        if ( !get_modified() ) {
            //evaluate temporary solution based on current k values for k1
            for (i=0; i<neq_; i++) {
                soltemp_[i] = sol_[i];
                for (m=0; m<nk_; m++) soltemp_[i] += dt*a[n][m]*k_[m][i];
            }
            //evaluate J at soltemp
            jac(soltemp_, Jac_);
        }
        //evaluate the Newton jacobian
        for (i=0; i<neq_; i++) {
            for (j=0; j<neq_; j++) {
                for (m=0; m<nk_; m++) {
                    J[i+n*neq_][j+m*neq_] = -dt*a[n][m]*Jac_[i][j];
                }
            }
            J[i+n*neq_][i+n*neq_] += 1.0; //subtract from the identity matrix
        }
    }
}

//------------------------------------------------------------------------------

OdeRadauIIA5::OdeRadauIIA5 (unsigned long neq) :
    OdeBase (neq, true),
    OdeIRK (neq, 3) {

    method_ = "RadauIIA5";

    int nk = 3;
    a = new double*[nk];
    for (int i=0; i<nk; i++) a[i] = new double[nk];
    b = new double[nk];

    double r = sqrt(6.0);

    a[0][0] =     (88 - 7*r)/360; a[0][1] = (296 - 169*r)/1800; a[0][2] = (-2 + 3*r)/225;
    a[1][0] = (296 + 169*r)/1800; a[1][1] =     (88 + 7*r)/360; a[1][2] = (-2 - 3*r)/225;
    a[2][0] =        (16 - r)/36; a[2][1] =        (16 + r)/36; a[2][2] =          1.0/9;

    b[0] =           (16 - r)/36; b[1] =           (16 + r)/36; b[2] =             1.0/9;

    newton_ = new NewtonRadauIIA5(neq, nk, this);
    newton_->set_modified(true);
}

OdeRadauIIA5::~OdeRadauIIA5 () {
    for (int i=0; i<nk_; i++) delete [] a[i];
    delete [] a;
    delete [] b;
}

void OdeRadauIIA5::step_ (double dt) {

    unsigned long i;
    //set k values to zero
    for (i=0; i<neq_*nk_; i++) kall_[i] = 0.0;
    //solve for slopes
    newton_->solve_Newton(kall_);
    //compute new solution
    for (i=0; i<neq_; i++) sol_[i] += dt*(b[0]*k_[0][i]
                                        + b[1]*k_[1][i]
                                        + b[2]*k_[2][i]);
}
