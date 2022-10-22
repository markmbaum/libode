#include <cstdio>

#include "ode_ssp_3.h"
#include "swarmalator.h"

int main () {

    //number of agents
    int nagent = 128;
    //choose an integration time
    double tend = 100;
    //choose time step
    double dt = tend/5000;

    //create integrator
    Swarmalator<OdeSsp3> sys(nagent);
    sys.set_name("swarmalator");

    //coupling parameters
    //uniform blob
    //sys.A = 1.0; sys.B = 1.0; sys.J = 0.5; sys.K = 0.5;
    //static
    //sys.A = 1.0; sys.B = 1.0; sys.J = 0.3; sys.K = -0.2;
    //circle
    sys.A = 1.0; sys.B = 1.0; sys.J = 1.0; sys.K = -0.2;

    printf("integrating...\n");
    sys.solve_fixed(tend, dt, "out", 15);
    printf("%lu function evaluations, %lu steps\n", sys.get_neval(), sys.get_nstep());

    return(0);
}
