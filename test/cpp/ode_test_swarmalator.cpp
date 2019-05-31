#include <cstdio>

#include "ode_explicit_test_systems.hpp"

int main () {

    //number of agents
    int N = 100;
    //choose an integration time
    double tend = 50;
    //number of snaps
    int snaps = int(tend)*10;
    //error tolerance
    double tol = 1e-4;
    //coupling parameters
    //double J=0.5, K=0.5;  // i
    double J=0.3, K=-0.2; // ii
    //double J=1.0, K=-0.2; // iii


    printf("starting swarmalators\n");
    Swarmalator sys(N);
    sys.abstol = tol; sys.reltol = tol;
    sys.J = J; sys.K = K;
    sys.dtmax = 0.49*tend/double(snaps);

    sys.solve_adaptive(tend, sys.dtmax, snaps, "out");

    printf("%llu function evaluations, %llu steps, %llu rejected steps\n",
        sys.neval, sys.nstep, sys.nrej);

    printf("\nswarmalators complete\n\n");
    return(0);
}
