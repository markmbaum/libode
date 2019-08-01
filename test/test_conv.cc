#include <cstdio>
#include <cmath>
#include <vector>

#include "test_systems.h"

//oscillator 1 equations
double y1exact (double t) {
    //Dahl
    //return(exp(-t));
    //oscillator 1
    return(cos(t*t/2.0));
    //oscillator 2
    //return(exp(sin(t*t)));
}
double y2exact (double t) {
    //Dahl
    //return(exp(-2*t));
    //oscillator 1
    return(sin(t*t/2.0));
    //oscillator 2
    //return(exp(cos(t*t)));
}

int main () {

    double h;
    double hmax = 1e-1;
    double hmin = 1e-4;
    double tint = 7;
    int iters = 0;
    double frac = 1.2;
    std::vector<double> y1err;
    std::vector<double> y2err;
    std::vector<double> y1;
    std::vector<double> y2;
    std::vector<double> hstore;

    printf("\n%14s %14s %14s %14s %14s %14s\n",
        "h", "tint", "y1", "y2", "y1_error", "y2_error");
        printf("%14s %14s %14s %14s %14s %14s\n",
            "----------", "----------", "----------",
            "----------", "----------", "----------");
    h = hmax;
    while (h >= hmin) {
        Osc1<OdeLobattoIIIC6> sys; //choose system and method
        sys.solve_fixed(tint, h);
        y1err.push_back( fabs(sys.get_sol(0) - y1exact(sys.get_t())) );
        y2err.push_back( fabs(sys.get_sol(1) - y2exact(sys.get_t())) );
        hstore.push_back(h);
        printf("%14g %14g %14g %14g %14g %14g\n",
                h,
                sys.get_t(),
                sys.get_sol(0),
                sys.get_sol(1),
                fabs(y1err.back()),
                fabs(y2err.back()) );
        h = h/frac;
        iters++;
    }

    ode_write("out/sol1err", y1err.data(), iters);
    ode_write("out/sol2err", y2err.data(), iters);
    ode_write("out/h", hstore.data(), iters);

    printf("\ntest complete\n");
    return(0);
}
