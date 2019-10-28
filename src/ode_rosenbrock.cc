//! \file ode_rosenbrock.cc

#include "ode_rosenbrock.h"

OdeRosenbrock::OdeRosenbrock (unsigned long neq, int nk) {

    //stage count
    nk_ = nk;
    //set up stage derivatives
    k_ = new double*[nk];
    for (int i=0; i<nk; i++) k_[i] = new double[neq];
    //allocate permutation array
    p_ = new int[neq];
    //right hand side of matrix equations
    rhs_ = new double[neq];
    //temporary sol vector
    soltemp_ = new double[neq];
}

OdeRosenbrock::~OdeRosenbrock () {
    for (int i=0; i<nk_; i++) delete [] k_[i];
    delete [] k_;
    delete [] p_;
    delete [] rhs_;
    delete [] soltemp_;
}

void OdeRosenbrock::prep_jac (double **Jac, unsigned long n, double dt, int *p) {

    unsigned long i,j;

    for (i=0; i<n; i++) {
        //flip the sign of the Jacobian across the row and multiply by gam*dt
        for (j=0; j<n; j++) Jac[i][j] = -Jac[i][j]*gam*dt;
        //subtract it from the identity matrix
        Jac[i][i] += 1.0;
    }
    ode_crout_LU(Jac, n, p);
}
