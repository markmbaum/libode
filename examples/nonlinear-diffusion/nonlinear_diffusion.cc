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
    double tint = 2e-1;

    //create system
    NonlinearDiffusion<OdeTrapz> sys(nx, L);
    sys.set_name("nonlinear_diffusion");
    sys.gam = 3.1;
    //initial conditions
    for (int i=0; i<nx; i++) sys.set_sol(i, 0.1 + 0.9*exp(-sys.xc[i]*sys.xc[i]*50));

    //integrate
    sys.solve_adaptive(tint, tint*1e-9, 15, "out");
    printf("%lu steps\n", sys.get_nstep());

    //write the grid to file
    ode_write("out/x", sys.xc, nx);

    return(0);
}
