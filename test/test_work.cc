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
void test_work (T sys, double tint, double frac, double dtmax, double dtmin) {

    int iters = 0;
    double dt = dtmax;
    std::vector<double> y1err, y2err, neval;
    double *ic = new double[sys.get_neq()];

    printf("  %s\n", sys.get_method());

    for (unsigned long i=0; i<sys.get_neq(); i++) ic[i] = sys.get_sol(i);

    while ( dt >= dtmin ) {
        sys.reset(0, ic);
        sys.solve_fixed(tint, dt);
        y1err.push_back( fabs(sys.get_sol(0) - y1exact(sys.get_t())) );
        y2err.push_back( fabs(sys.get_sol(1) - y2exact(sys.get_t())) );
        neval.push_back( double(sys.get_neval()) );
        dt *= frac;
        iters++;
    }

    std::string method = sys.get_method();
    ode_write(("out/sol1err_" + method).data(), y1err.data(), iters);
    ode_write(("out/sol2err_" + method).data(), y2err.data(), iters);
    ode_write(("out/neval_" + method).data(), neval.data(), iters);
}

int main () {

    //total solution time
    double tint = 12.0;
    //fraction time step reduction
    double frac = 0.9;

    printf("finished:\n");
    test_work(Osc2<OdeEuler>(), tint, frac, 5e-4, 1e-6);
    test_work(Osc2<OdeTrapz>(), tint, frac, 3e-2, 1e-5);
    test_work(Osc2<OdeSsp3>(), tint, frac, 3e-2, 3e-5);
    test_work(Osc2<OdeRK4>(), tint, frac, 5e-2, 5e-4);
    test_work(Osc2<OdeDoPri54>(), tint, frac, 3e-2, 8e-4);
    test_work(Osc2<OdeVern65>(), tint, frac, 5e-2, 1e-3);
    test_work(Osc2<OdeVern76>(), tint, frac, 5e-2, 2e-3);
    test_work(Osc2<OdeDoPri87>(), tint, frac, 1e-1, 5e-3);
    test_work(Osc2<OdeVern98>(), tint, frac, 1e-1, 5e-3);

    return(0);
}
