#include <cstdio>

#include "test_systems.h"

int main () {

    //choose an integration time
    double tint = 5;
    //time step
    double dt = tint/1000.0;
    //choose the number of snaps
    int nsnap = 5;
    //choose which system and solver to use
    Osc1<OdeVern65> sys;

    double *tsnap = new double[5];
    tsnap[0] = 1.1;
    tsnap[1] = 2.22;
    tsnap[2] = 3;
    tsnap[3] = 3.14;
    tsnap[4] = 4.9;

    printf("Solving system '%s' with method '%s'\n",
        sys.get_name(),
        sys.get_method());

    /* solve with evenly spaced snaps */
    //sys.solve_adaptive(tint, dt, nsnap, "out");

    /* solve with snap times from tsnap */
    sys.solve_adaptive(dt, tsnap, nsnap, "out");

    printf("%lu function evaluations, %lu steps\n",
        sys.get_neval(),
        sys.get_nstep());

    delete [] tsnap;
    return(0);
}
