#include <cmath>

#include "ode_newton.h"

using namespace ode;

//function to find root of
void f_f (double *x, double *y) {
    y[0] = 3*x[0] - cos(x[1]*x[2]) - 0.5;
    y[1] = x[0]*x[0] - 81.0*(x[1] + 0.1)*(x[1] + 0.1) + sin(x[2]) + 1.06;
    y[2] = exp(-x[0]*x[1]) + 20.0*x[2] + (10*M_PI - 3.0)/3.0;
}

//Jacobian of f
void f_J (double *x, double **J) {
    J[0][0] = 3.0;
    J[0][1] = sin(x[1]*x[2])*x[1];
    J[0][2] = sin(x[1]*x[2])*x[2];

    J[1][0] = 2*x[0];
    J[1][1] = -162.0*(x[1] + 0.1);
    J[1][2] = -cos(x[2]);

    J[2][0] = -x[1]*exp(-x[0]*x[1]);
    J[2][1] = -x[0]*exp(-x[0]*x[1]);
    J[2][2] = 20.0;
}

class TestNewton : public OdeNewton {
public:
    TestNewton () : OdeNewton (3) {}
    void f_Newton (double *x, double *y) { f_f(x, y); }
    void J_Newton (double *x, double **J) { f_J(x, J); }
};

int main () {

    //system and solver
    TestNewton sys;
    sys.set_modified(true);
    //initial values
    double *x = new double[3];

    //print the solution
    printf("Exact solution:\n");
    printf("x0 = %16g\nx1 = %16g\nx2 = %16g\n", 0.5, 0.0, -0.52359877);

    //set initial values
    x[0] = 1.0; x[1] = 1.0; x[2] = 1.0;
    //solve with class functions
    sys.solve_Newton(x);
    //print results
    printf("Newton solution:\n");
    printf("x0 = %16g\nx1 = %16g\nx2 = %16g\n", x[0], x[1], x[2]);
    printf("nJLU = %lu\n", sys.get_nJLU());

    //finish
    delete [] x;

    return(0);
}
