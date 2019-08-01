#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include "test_systems.h"

double y1exact (double t) {

    //oscillator 1
    //return(cos(t*t/2));

    //oscillator 2
    return(exp(sin(t*t)));
}

double y2exact (double t) {

    //oscillator 1
    //return(sin(t*t/2));

    //oscillator 2
    return(exp(cos(t*t)));
}

template<class T>
void test_work (T sys, double tint, double frac, double dtmax, double dtmin, double *ic, const char *name) {

    int iters = 0;
    double dt = dtmax;
    std::vector<double> y1err, y2err, neval;

    while ( dt >= dtmin ) {
        sys.reset(0.0, ic);
        sys.solve_fixed(tint, dt);
        y1err.push_back( fabs(sys.get_sol(0) - y1exact(sys.get_t())) );
        y2err.push_back( fabs(sys.get_sol(1) - y2exact(sys.get_t())) );
        neval.push_back( double(sys.get_neval()) );
        dt *= frac;
        iters++;
    }

    std::string name_ = name;
    ode_write(("out/sol1err_" + name_).data(), y1err.data(), iters);
    ode_write(("out/sol2err_" + name_).data(), y2err.data(), iters);
    ode_write(("out/neval_" + name_).data(), neval.data(), iters);
}

int main () {

    //total solution time
    double tint = 12.0;
    //fraction time step reduction
    double frac = 0.75;
    //initial conditions for resetting
    double ic[2] = {1.0, exp(1.0)};

    Osc2<OdeEuler> euler;
    test_work(euler, tint, frac, 5e-4, 1e-6, ic, "Euler");
    printf("Euler\n");

    Osc2<OdeTrapz> trapz;
    test_work(trapz, tint, frac, 3e-2, 1e-5, ic, "Trapz");
    printf("Trapz\n");

    Osc2<OdeSsp3> ssp3;
    test_work(ssp3, tint, frac, 3e-2, 3e-5, ic, "Ssp3");
    printf("Ssp3\n");

    Osc2<OdeRK4> rk4;
    test_work(rk4, tint, frac, 5e-2, 5e-4, ic, "RK4");
    printf("RK4\n");

    Osc2<OdeDoPri54> dopri54;
    test_work(dopri54, tint, frac, 3e-2, 8e-4, ic, "DoPri54");
    printf("DoPri54\n");

    Osc2<OdeVern65> vern65;
    test_work(vern65, tint, frac, 5e-2, 1e-3, ic, "Vern65");
    printf("Vern65\n");

    Osc2<OdeVern76> vern76;
    test_work(vern76, tint, frac, 5e-2, 2e-3, ic, "Vern76");
    printf("Vern76\n");

    Osc2<OdeDoPri87> dopri87;
    test_work(dopri87, tint, frac, 1e-1, 5e-3, ic, "DoPri87");
    printf("DoPri87\n");

    Osc2<OdeVern98> vern98;
    test_work(vern98, tint, frac, 1e-1, 5e-3, ic, "Vern98");
    printf("Vern98\n");

    return(0);
}
