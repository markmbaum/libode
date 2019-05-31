#include <cstdio>
#include <cmath>
#include <vector>

#include "ode_explicit_test_systems.hpp"
#include "ode_misc.hpp"

//oscillator 1 equations
double y1exact (double t) {
    //oscillator 1
    return(cos(t*t/2));
    //oscillator 2
    //return(exp(sin(t*t)));
}
double y2exact (double t) {
    //oscillator 1
    return(sin(t*t/2));
    //oscillator 2
    //return(exp(cos(t*t)));
}

int main () {

    double h;
    double hmax = 1e-1;
    double hmin = 1e-5;
    double tend = 5;
    int iters = 0;
    double frac = 1.25;
    std::vector<double> y1err;
    std::vector<double> y2err;
    std::vector<double> y1;
    std::vector<double> y2;
    std::vector<double> hstore;

    printf("%14s %14s %14s %14s %14s %14s\n",
        "h", "tend", "y1", "y2", "y1_error", "y2_error");
        printf("%14s %14s %14s %14s %14s %14s\n",
            "----------", "----------", "----------",
            "----------", "----------", "----------");
    h = hmax;
    while (h >= hmin) {
        Ssp3_Osc1 sys; //choose system and method
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        hstore.push_back(h);
        printf("%14g %14g %14g %14g %14g %14g\n",
                h,
                sys.t,
                sys.sol[0],
                sys.sol[1],
                fabs(y1err.back()),
                fabs(y2err.back())
        );
        h = h/frac;
        iters++;
    }

    ode_write_double("out/sol1err", y1err.data(), iters);
    ode_write_double("out/sol2err", y2err.data(), iters);
    ode_write_double("out/h", hstore.data(), iters);

    printf("\ntest complete\n");
    return(0);
}
