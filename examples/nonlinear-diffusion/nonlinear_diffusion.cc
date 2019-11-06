#include <cstdio>
#include <cmath>

#include "ode_trapz.h"
#include "nonlinear_diffusion.h"

int main () {

    //number of spatial cells
    int nx = 1000;
    //domain length
    double L = 0.5;
    //integration time
    double tint = 1e-1;
    //time step
    double dt = tint/1e5;

    //create system
    NonlinearDiffusion<OdeTrapz> sys(nx, L);
    sys.set_name("nonlinear_diffusion");
    sys.gam = 3.1;
    //initial conditions
    for (int i=0; i<nx; i++) sys.set_sol(i, 0.1 + 0.9*exp(-sys.xc[i]*sys.xc[i]*50));

    //integrate
    sys.set_prescribe_adapt(true);
    sys.solve_fixed(tint, dt, 15, "out");
    printf("%llu steps\n", sys.get_nstep());

    //write the grid to file
    ode_write("out/x", sys.xc, nx);

    return(0);
}
