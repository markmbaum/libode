#include <cstdio>

#include "test_systems.h"

int main () {

    //choose an integration time
    double tint = 50;
    //time step
    double dt = 0.01;
    //choose which system and solver to use
    Brus<OdeSDIRK43> sys;
    //sys.set_tol(1e-3);

    printf("Solving system '%s' with method '%s'\n",
        sys.get_name(), sys.get_method());
    sys.solve_adaptive(tint, dt, "out", 1);
    printf("%llu function evaluations\n", sys.get_neval());
    printf("%llu steps\n", sys.get_nstep());
    printf("%llu rejected steps\n", sys.get_nrej());
    printf("%llu ODE Jacobian evaluations\n", sys.get_nJac());
    printf("%lu Newton JLU updates (LU decompositions)\n", sys.get_newton()->get_nJLU());
    printf("%lu Newton LU solves (matrix equation solves)\n", sys.get_newton()->get_n_solve_LU());

    return(0);
}
