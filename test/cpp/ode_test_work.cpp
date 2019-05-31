#include <cstdio>
#include <cmath>
#include <vector>

#include "ode_explicit_test_systems.hpp"
#include "ode_misc.hpp"

//oscillator 1 equations
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

int main () {

    int iters;
    double tend, frac, hmax, h, hmin, tolmax, tol, tolmin;
    std::vector<double> y1err, y2err, neval;

    iters = 0;
    tend = 10, frac = 1.25, hmax = 5e-4, h = hmax, hmin = 1e-5;
    y1err.clear(); y2err.clear(); neval.clear();
    while (h >= hmin) {
        Euler_Osc2 sys;
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        h = h/frac;
        iters++;
    }
    ode_write_double("out/sol1err_euler", y1err.data(), iters);
    ode_write_double("out/sol2err_euler", y2err.data(), iters);
    ode_write_double("out/neval_euler", neval.data(), iters);

    iters = 0;
    tend = 10, frac = 1.25, hmax = 3e-2, h = hmax, hmin = 1e-5;
    y1err.clear(); y2err.clear(); neval.clear();
    while (h >= hmin) {
        Midpnt_Osc2 sys;
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        h = h/frac;
        iters++;
    }
    ode_write_double("out/sol1err_midpnt", y1err.data(), iters);
    ode_write_double("out/sol2err_midpnt", y2err.data(), iters);
    ode_write_double("out/neval_midpnt", neval.data(), iters);

    iters = 0;
    tend = 10, frac = 1.25, hmax = 3e-2, h = hmax, hmin = 1e-5;
    y1err.clear(); y2err.clear(); neval.clear();
    while (h >= hmin) {
        Trapz_Osc2 sys;
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        h = h/frac;
        iters++;
    }
    ode_write_double("out/sol1err_trapz", y1err.data(), iters);
    ode_write_double("out/sol2err_trapz", y2err.data(), iters);
    ode_write_double("out/neval_trapz", neval.data(), iters);

    iters = 0;
    tolmax = 1e-3, tol = tolmax, tolmin = 1e-10, tend = 10, frac = 1.25, hmax = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (tol >= tolmin) {
        RKF32_Osc2 sys;
        sys.abstol = tol; sys.reltol = tol;
        sys.solve_adaptive(tend, hmax);
        printf("rkf23 tol=%g, nrej=%llu\n", tol, sys.nrej);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        tol = tol/frac;
        iters++;
    }
    ode_write_double("out/sol1err_rkf32", y1err.data(), iters);
    ode_write_double("out/sol2err_rkf32", y2err.data(), iters);
    ode_write_double("out/neval_rkf32", neval.data(), iters);

    iters = 0;
    tend = 10, frac = 1.15, hmax = 5e-2, h = hmax, hmin = 5e-4;
    y1err.clear(); y2err.clear(); neval.clear();
    while (h >= hmin) {
        RK4_Osc2 sys;
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        h = h/frac;
        iters++;
    }
    ode_write_double("out/sol1err_rk4", y1err.data(), iters);
    ode_write_double("out/sol2err_rk4", y2err.data(), iters);
    ode_write_double("out/neval_rk4", neval.data(), iters);

    iters = 0;
    tolmax = 1e-3, tol = tolmax, tolmin = 1e-11, tend = 10, frac = 1.25, hmax = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (tol >= tolmin) {
        RK43_Osc2 sys;
        sys.abstol = tol; sys.reltol = tol;
        sys.solve_adaptive(tend, hmax);
        printf("rk34 tol=%g, nrej=%llu\n", tol, sys.nrej);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        tol = tol/frac;
        iters++;
    }
    ode_write_double("out/sol1err_rk43", y1err.data(), iters);
    ode_write_double("out/sol2err_rk43", y2err.data(), iters);
    ode_write_double("out/neval_rk43", neval.data(), iters);

    iters = 0;
    tolmax = 1e-3, tol = tolmax, tolmin = 1e-11, tend = 10, frac = 1.25, hmax = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (tol >= tolmin) {
        RKCK_Osc2 sys;
        sys.abstol = tol; sys.reltol = tol;
        sys.solve_adaptive(tend, hmax);
        printf("rkck tol=%g, nrej=%llu\n", tol, sys.nrej);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        tol = tol/frac;
        iters++;
    }
    ode_write_double("out/sol1err_rkck", y1err.data(), iters);
    ode_write_double("out/sol2err_rkck", y2err.data(), iters);
    ode_write_double("out/neval_rkck", neval.data(), iters);

    iters = 0;
    tolmax = 1e-3, tol = tolmax, tolmin = 1e-12, tend = 10, frac = 1.25, hmax = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (tol >= tolmin) {
        DoPri54_Osc2 sys;
        sys.abstol = tol; sys.reltol = tol;
        sys.solve_adaptive(tend, hmax);
        printf("dopri45 tol=%g, nrej=%llu\n", tol, sys.nrej);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        tol = tol/frac;
        iters++;
    }
    ode_write_double("out/sol1err_dopri54", y1err.data(), iters);
    ode_write_double("out/sol2err_dopri54", y2err.data(), iters);
    ode_write_double("out/neval_dopri54", neval.data(), iters);

    iters = 0;
    tend = 10, frac = 1.15, hmax = 5e-2, h = hmax, hmin = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (h >= hmin) {
        Butcher6_Osc2 sys;
        sys.solve_fixed(tend, h);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        h = h/frac;
        iters++;
    }
    ode_write_double("out/sol1err_butcher6", y1err.data(), iters);
    ode_write_double("out/sol2err_butcher6", y2err.data(), iters);
    ode_write_double("out/neval_butcher6", neval.data(), iters);

    iters = 0;
    tolmax = 1e-3, tol = tolmax, tolmin = 2e-11, tend = 10, frac = 1.25, hmax = 1e-3;
    y1err.clear(); y2err.clear(); neval.clear();
    while (tol >= tolmin) {
        DoPri87_Osc2 sys;
        sys.abstol = tol; sys.reltol = tol;
        sys.solve_adaptive(tend, hmax);
        printf("dopri78 tol=%g, nrej=%llu\n", tol, sys.nrej);
        y1err.push_back(fabs(sys.sol[0] - y1exact(sys.t)));
        y2err.push_back(fabs(sys.sol[1] - y2exact(sys.t)));
        neval.push_back(double(sys.neval));
        tol = tol/frac;
        iters++;
    }
    ode_write_double("out/sol1err_dopri87", y1err.data(), iters);
    ode_write_double("out/sol2err_dopri87", y2err.data(), iters);
    ode_write_double("out/neval_dopri87", neval.data(), iters);

    return(0);
}
