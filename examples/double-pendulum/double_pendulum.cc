#include <cstdio>

#include "ode_gauss_6.h"
#include "double_pendulum.h"

using namespace ode;

int main () {

    //create integrator
    DoublePendulum<OdeGauss6> sys;
    sys.set_name("double_pendulum");

    //physical parameters
    sys.g = 9.8, //gravity (m/s^2)
    sys.m1 = 1.0, //mass 1 (kg)
    sys.m2 = 1.0, //mass 2 (kg)
    sys.L1 = 1.0, //length 1 (m)
    sys.L2 = 1.5; //length 2 (m)

    //initial conditions
    sys.set_sol(0,  2); //initial angle of first mass
    sys.set_sol(1,  0); //initial angular momentum of first mass
    sys.set_sol(2, -2); //initial angle of second mass
    sys.set_sol(3,  0); //initial angular momentum of second mass

    //integration parameters
    double tend = 25;
    double dt = 1e-4;

    //integrate
    printf("integrating system '%s' with method '%s'...\n", sys.get_name(), sys.get_method());
    sys.solve_fixed(tend, dt, "out", 125);
    printf("finished.\n%lu steps, %lu function evaluations\n",
        sys.get_nstep(), sys.get_neval());

    return(0);
}
