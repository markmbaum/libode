#include <cstdio>

#include "test_systems.h"

using namespace ode;

int main () {

    //choose an integration time
    double tint = 5;
    //time step
    double dt = 0.05;
    //choose which system and solver to use
    Osc1<OdeRadauIIA5> sys;

    printf("Solving system '%s' with method '%s'\n",
        sys.get_name(),
        sys.get_method());
    sys.solve_fixed(tint, dt, "out", 1);
    printf("%lu function evaluations, %lu steps\n",
        sys.get_neval(),
        sys.get_nstep());

    return(0);
}
