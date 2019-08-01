#include <cstdio>

#include "ode_dopri_54.h"
#include "lorentz.h"

int main () {

    //desired solution parameters
    double tend = 20;
    double dt = tend/10000;

    //construct the system and integrator
    Lorentz<OdeDoPri54> sys;
    sys.set_name("lorentz");
    sys.sigma = 10.0;
    sys.beta = 8.0/3;
    sys.rho = 28.0;
    //initial conditions
    sys.set_sol(0, 1.0);
    sys.set_sol(1, 1.0);
    sys.set_sol(2, 1.0);

    //integrate and write results into the "out" directory
    printf("solving...\n");
    sys.solve_fixed(tend, dt, "out", 1);

    printf("%llu function evaluations, %llu steps\n", sys.get_neval(), sys.get_nstep());

    return(0);
}
