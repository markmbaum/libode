#include <cstdio>

#include "ode_explicit_test_systems.hpp"

int main () {

    //choose an integration time
    double tend = 6.0;
    //choose error tolerance
    double abstol = 1e-9;
    double reltol = abstol;
    //output directory
    char dirout[4] = "out";
    //choose which system and solver to use
    DoPri54_Osc2 sys;

    printf("starting test\n");
    sys.abstol = abstol;
    sys.reltol = reltol;
    sys.solve_adaptive(tend, tend*1e-3, dirout);
    printf("%llu function evaluations, %llu steps, %llu rejected steps\n",
        sys.neval, sys.nstep, sys.nrej);

    printf("\ntest complete\n");
    return(0);
}
