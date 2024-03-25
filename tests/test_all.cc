/*
This program tests all the integrators on a simple test problem with known solution to see if any of them are wildly incorrect. Each method should also be tested for the expected order of convergence.
*/

#include <cstdio>
#include <cmath>

#include "test_systems.h"

using namespace ode;

template<class T>
void test_sys (T sys) {

    //integration time
    double tint = 3;
    //small time step
    double dt = 0.001;
    //error "tolerance" sort of
    double tol = 0.1;
    //exact solutions for oscilator 1
    double y0ex = exp(sin(tint*tint)),
           y1ex = exp(cos(tint*tint));

    double y0er, y1er;
    //integrate the system
    sys.solve_fixed(tint, dt);
    //calculate errors
    y0er = fabs(sys.get_sol(0) - y0ex);
    y1er = fabs(sys.get_sol(1) - y1ex);
    //check the errors
    printf("%s:\n", sys.get_method());
    printf("  y0 error = %12e", y0er);
    if ( y0er > tol ) {
        printf(" <---- DANGER!!!\n");
        exit(EXIT_FAILURE);
    } else {
        printf("\n");
    }
    printf("  y1 error = %12e", y1er);
    if ( y1er > tol ) {
        printf(" <---- DANGER!!!\n");
        exit(EXIT_FAILURE);
    } else {
        printf("\n");
    }
}

int main () {

    //integrate all the systems and check their errors
    test_sys(Osc2<OdeEuler>());
    test_sys(Osc2<OdeTrapz>());
    test_sys(Osc2<OdeSsp3>());
    test_sys(Osc2<OdeRKF32>());
    test_sys(Osc2<OdeRK4>());
    test_sys(Osc2<OdeRK43>());
    test_sys(Osc2<OdeRKCK>());
    test_sys(Osc2<OdeDoPri54>());
    test_sys(Osc2<OdeVern65>());
    test_sys(Osc2<OdeVern76>());
    test_sys(Osc2<OdeDoPri87>());
    test_sys(Osc2<OdeVern98>());
    test_sys(Osc2<OdeGRK4A>());
    test_sys(Osc2<OdeROW6A>());
    test_sys(Osc2<OdeBackwardEuler>());
    test_sys(Osc2<OdeGauss6>());
    test_sys(Osc2<OdeLobattoIIIC6>());
    test_sys(Osc2<OdeRadauIIA5>());
    test_sys(Osc2<OdeGeng5>());
    test_sys(Osc2<OdeSDIRK43>());

}
