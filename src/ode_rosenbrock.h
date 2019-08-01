#ifndef ODE_ROSENBROCK_H_
#define ODE_ROSENBROCK_H_

#include "ode_linalg.h"

class OdeRosenbrock {

    public:
        //constructor
        OdeRosenbrock (unsigned long neq, int nk);
        //destructor
        ~OdeRosenbrock ();

    protected:
        //parameter multipying Jacobian or single diagonal gamma
        double gam;
        //permutation array for LU factorization
        int *p_;
        //right hand side of matrix equations
        double *rhs_;
        //temporary sol vector
        double *soltemp_;
        //stage derivatives
        double **k_;
        //do necessary arithmetic with the Jacobian then lu factor it
        void prep_jac (double **Jac, unsigned long n, double dt, int *p);

    private:
        //number of stages
        int nk_;
};

#endif
