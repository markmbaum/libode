#include <cstdio>

#include "test_systems.h"

int main () {

    //choose an integration time
    double tint = 5;
    //choose which system and solver to use
    Osc2<OdeDoPri54> sys;
    //choose error tolerance
    sys.set_tol(1e-6);

    printf("Solving system '%s' with method '%s'\n",
        sys.get_name(), sys.get_method());
    sys.solve_adaptive(tint, tint/10.0, "out", 1);
    printf("%d function evaluations, %d steps, %d rejected steps\n",
        int(sys.get_neval()),
        int(sys.get_nstep()),
        int(sys.get_nrej()));

    return(0);
}
