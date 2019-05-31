#include <cstdio>
#include <cmath>
#include <vector>

#include "ode_explicit_test_systems.hpp"

int main () {

    //choose an integration time
    double tend = 7;
    //time step
    double dt = 1e-3;
    //output directory
    char dirout[4] = "out";
    //choose which system and solver to use
    RKCK_Osc1 sys;

    sys.solve_fixed(tend, dt, dirout);
    printf("%llu function evaluations, %llu steps\n", sys.neval, sys.nstep);

    printf("\ntest complete\n");
    return(0);
}
