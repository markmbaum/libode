#include <cstdio>

#include "ode_geng_5.h"
#include "star.h"

int main () {

    //create integrator
    Star<OdeGeng5> sys;
    //Star<OdeGeng5> sys;
    sys.set_name("star");

    //set parameters
    sys.a = 1.25;
    sys.b = 1.0;
    sys.c = 0.75;
    sys.A = 1.0;
    sys.C = 1.0;
    sys.Omega = 0.25;

    //set initial conditions
    sys.q[0] = 2.5;
    sys.q[1] = 0.0;
    sys.q[2] = 0.0;
    sys.p[0] = 0.0;
    sys.p[1] = 1.6888370059044742;
    sys.p[2] = 0.2;

    //integration parameters
    double tend = 1000;
    double dt = 1.0/50;

    //integrate
    printf("integrating...\n");
    sys.solve_fixed(tend, dt, "out", 1);
    ode_write("out/star_H", sys.H.data(), sys.H.size());
    printf("finished.\n%lu steps, %lu function evaluations\n",
        sys.get_nstep(), sys.get_neval());

    return(0);
}
