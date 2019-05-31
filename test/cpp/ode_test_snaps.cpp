#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include "ode_explicit_test_systems.hpp"

int main () {

    //choose an integration time
    double tend = 10;
    //time step
    double dt = tend/1000;
    //output directory
    std::string dirout = "out";
    //choose the number of snaps
    int snaps = 20;
    //choose which system and solver to use
    RK4_Osc1 sys;

    sys.solve_fixed(tend, dt, snaps, dirout.c_str());
    printf("%llu function evaluations, %llu steps\n", sys.neval, sys.nstep);

    //write the number of snaps
    std::string fn = dirout + "/nsnap.txt";
    ode_check_file_write(fn.c_str());
    FILE *ofile;
    ofile = fopen(fn.c_str(), "w");
    fprintf(ofile, "%d", snaps);
    fclose(ofile);
    printf("file written: %s\n", fn.c_str());

    printf("\ntest complete\n");
    return(0);
}
