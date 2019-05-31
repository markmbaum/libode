/*
This file defines test systems of equations for each of the explicit, single-step Runge-Kutta methods. It uses three different systems, two coupled oscillators and the _Dahlquist test system (y' = a*t). The first two systems, osc1 and osc2, have known analytical solutions.
*/

#ifndef ODE_TEST_SYSTEMS_HPP
#define ODE_TEST_SYSTEMS_HPP

#include <cmath>
#include <cstdlib>

#include "ode_misc.hpp"
#include "OdeEuler.hpp"
#include "OdeMidpnt.hpp"
#include "OdeTrapz.hpp"
#include "OdeRKF21.hpp"
#include "OdeSsp3.hpp"
#include "OdeRKF32.hpp"
#include "OdeRK4.hpp"
#include "OdeRK43.hpp"
#include "OdeRKCK.hpp"
#include "OdeDoPri54.hpp"
#include "OdeButcher6.hpp"
#include "OdeDoPri87.hpp"

//------------------------------------------------------------------------------
//test systems of ODEs

//This is a system of two oscillating equations, which become increasingly
//steep with time. The solutions are
//  y1 = cos(t^2/2),   y1(0) = 1
//  y2 = sin(t^2/2),   y2(0) = 0
void osc1 (double t_, double *solin, double *fout) {
    //coupled oscillators
    fout[0] = -t_*solin[1];
    fout[1] = t_*solin[0];
}
//This is a second system of two oscillating equations, taken from Hairer II.4
//The solutions are:
//  y1 = cos(t^2/2),   y1(0) = 1
//  y2 = sin(t^2/2),   y2(0) = e
void osc2 (double t_, double *solin, double *fout) {
    //other coupled oscillators
    solin[1] > 1e-3 ? fout[0] = 2*t_*solin[0]*log(solin[1])
                  : fout[0] = 2*t_*solin[0]*log(1e-3);
    solin[0] > 1e-3 ? fout[1] = -2*t_*solin[1]*log(solin[0])
                  : fout[1] = -2*t_*solin[1]*log(1e-3);
}

//_Dahlquist test equation, the solution of which is just an exponential
void dahl (double t_, double *solin, double *fout) {
    (void)t_; //supress unused variable warning
    fout[0] = -solin[0];
    fout[1] = -2*solin[1];
}

//This is a classic test system of two coupled, nonlinear equations called
//the _Brusselator (_Brussels oscillator).
void brusselator (double t_, double *solin, double *fout) {
    (void)t_; //supress unused variable warning
    fout[0] = 1.0 + solin[0]*solin[0]*solin[1] - 4.0*solin[0];
    fout[1] = 3.0*solin[0] - solin[0]*solin[0]*solin[1];
}

//------------------------------------------------------------------------------
//other stuff

//quick function for getting a pseudorandom number between 1 and 0
double rand_ () { return(double(rand())/double(RAND_MAX)); }

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//classes defining test systems for many different methods

//------------------------------------------------------------------------------
//Euler's method

class Euler_Dahl : public OdeEuler {
public:
    Euler_Dahl () : OdeEuler(2) { sol[0] = 1.0; sol[1] = 1.0; };
    ~Euler_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl(t_, solin, fout);
    }
};
class Euler_Osc1 : public OdeEuler {
public:
    Euler_Osc1 () : OdeEuler(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~Euler_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class Euler_Osc2 : public OdeEuler {
public:
    Euler_Osc2 () : OdeEuler(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~Euler_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//midpoint rule

class Midpnt_Dahl : public OdeMidpnt {
public:
    Midpnt_Dahl () : OdeMidpnt(2) { sol[0] = 1.0; sol[1] = 1.0; };
    ~Midpnt_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl(t_, solin, fout);
    }
};
class Midpnt_Osc1 : public OdeMidpnt {
public:
    Midpnt_Osc1 () : OdeMidpnt(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~Midpnt_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class Midpnt_Osc2 : public OdeMidpnt {
public:
    Midpnt_Osc2 () : OdeMidpnt(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~Midpnt_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//trapzedoial rule

class Trapz_Dahl : public OdeTrapz {
public:
    Trapz_Dahl () : OdeTrapz(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~Trapz_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class Trapz_Osc1 : public OdeTrapz {
public:
    Trapz_Osc1 () : OdeTrapz(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~Trapz_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class Trapz_Osc2 : public OdeTrapz {
public:
    Trapz_Osc2 () : OdeTrapz(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~Trapz_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//Fehlberg order 1 and 2

class RKF21_Dahl : public OdeRKF21 {
public:
    RKF21_Dahl () : OdeRKF21(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~RKF21_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class RKF21_Osc1 : public OdeRKF21 {
public:
    RKF21_Osc1 () : OdeRKF21(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~RKF21_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class RKF21_Osc2 : public OdeRKF21 {
public:
    RKF21_Osc2 () : OdeRKF21(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~RKF21_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
/*strong stability preserving rule of order 3, from Shu and Osher

    Shu, C.W. and Osher, S., 1988. Efficient implementation of essentially non-oscillatory shock-capturing schemes. Journal of computational physics, 77(2), pp.439-471.
*/

class Ssp3_Dahl : public OdeSsp3 {
public:
    Ssp3_Dahl () : OdeSsp3(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~Ssp3_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class Ssp3_Osc1 : public OdeSsp3 {
public:
    Ssp3_Osc1 () : OdeSsp3(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~Ssp3_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class Ssp3_Osc2 : public OdeSsp3 {
public:
    Ssp3_Osc2 () : OdeSsp3(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~Ssp3_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//third and fourth order embedded with fsal

class RK43_Dahl : public OdeRK43 {
public:
    RK43_Dahl () : OdeRK43(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~RK43_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class RK43_Osc1 : public OdeRK43 {
public:
    RK43_Osc1 () : OdeRK43(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~RK43_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class RK43_Osc2 : public OdeRK43 {
public:
    RK43_Osc2 () : OdeRK43(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~RK43_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//the classic fourth order Runge-Kutta method

class RK4_Dahl : public OdeRK4 {
public:
    RK4_Dahl () : OdeRK4(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~RK4_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class RK4_Osc1 : public OdeRK4 {
public:
    RK4_Osc1 () : OdeRK4(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~RK4_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1(t_, solin, fout);
    }
};
class RK4_Osc2 : public OdeRK4 {
public:
    RK4_Osc2 () : OdeRK4(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~RK4_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//second and third order embedded from Fehlberg

class RKF32_Dahl : public OdeRKF32 {
public:
    RKF32_Dahl () : OdeRKF32(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~RKF32_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class RKF32_Osc1 : public OdeRKF32 {
public:
    RKF32_Osc1 () : OdeRKF32(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~RKF32_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1 (t_, solin, fout);
    }
};
class RKF32_Osc2 : public OdeRKF32 {
public:
    RKF32_Osc2 () : OdeRKF32(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~RKF32_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//Cash-Karp with 5th, 4th, 3rd, 2nd, and 1st order solutions, a full family

class RKCK_Dahl : public OdeRKCK {
public:
    RKCK_Dahl () : OdeRKCK(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~RKCK_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class RKCK_Osc1 : public OdeRKCK {
public:
    RKCK_Osc1 () : OdeRKCK(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~RKCK_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1 (t_, solin, fout);
    }
};
class RKCK_Osc2 : public OdeRKCK {
public:
    RKCK_Osc2 () : OdeRKCK(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~RKCK_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//The much used fourth and fifth order solver from Dormand and Prince

class DoPri54_Dahl : public OdeDoPri54 {
public:
    DoPri54_Dahl () : OdeDoPri54(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~DoPri54_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class DoPri54_Osc1 : public OdeDoPri54 {
public:
    DoPri54_Osc1 () : OdeDoPri54(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~DoPri54_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1 (t_, solin, fout);
    }
};
class DoPri54_Osc2 : public OdeDoPri54 {
public:
    DoPri54_Osc2 () : OdeDoPri54(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~DoPri54_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//A sixth order method from Butcher

class Butcher6_Dahl : public OdeButcher6 {
public:
    Butcher6_Dahl () : OdeButcher6(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~Butcher6_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class Butcher6_Osc1 : public OdeButcher6 {
public:
    Butcher6_Osc1 () : OdeButcher6(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~Butcher6_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1 (t_, solin, fout);
    }
};
class Butcher6_Osc2 : public OdeButcher6 {
public:
    Butcher6_Osc2 () : OdeButcher6(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~Butcher6_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
//The eighth and seventh order solver from Dormand and Prince, as recommended
// and reprinted by Hairer-Norsett-Wanner

class DoPri87_Dahl : public OdeDoPri87 {
public:
    DoPri87_Dahl () : OdeDoPri87(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~DoPri87_Dahl () {};
    void ode_funk (double t_, double *solin, double *fout) {
        dahl (t_, solin, fout);
    }
};
class DoPri87_Osc1 : public OdeDoPri87 {
public:
    DoPri87_Osc1 () : OdeDoPri87(2) { sol[0] = 1.0; sol[1] = 0.0; };
    ~DoPri87_Osc1 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc1 (t_, solin, fout);
    }
};
class DoPri87_Osc2 : public OdeDoPri87 {
public:
    DoPri87_Osc2 () : OdeDoPri87(2) { sol[0] = 1.0; sol[1] = exp(1.0); };
    ~DoPri87_Osc2 () {};
    void ode_funk (double t_, double *solin, double *fout) {
        osc2 (t_, solin, fout);
    }
};
class DoPri87_Brus : public OdeDoPri87 {
public:
    DoPri87_Brus () : OdeDoPri87(2) { sol[0] = 1.5; sol[1] = 3.0; };
    ~DoPri87_Brus () {};
    void ode_funk (double t_, double *solin, double *fout) {
        brusselator (t_, solin, fout);
    }
};

//------------------------------------------------------------------------------
/*
Oscillators that sync and swarm, or SWARMALATORS, from:

    3K. P. O’Keeffe, H. Hong, and S. H. Strogatz, Oscillators that sync and swarm, Nat. Commun. 8, 1504 (2017). doi:10.1038/s41467-017-01190-3

Evaluation of the system could be sped up by pre-calculating distances to avoid duplication, but oh well...
*/

class Swarmalator : public OdeRKF32 {
public:
    //number of agents
    int n;
    //constants
    double A, B, J, K, N;

    //each agent gets three equations, x, y, and phase
    Swarmalator (int n_) : OdeRKF32 (3*n_) {
        //number of agents
        n = n_;
        //constants
        A = 1.0;
        B = 1.0;
        N = double(n);
        J = 1.0;
        K = -0.2;
        //initial x coord, y coord, and phase (theta)
        double r1, r2;
        for (int i=0; i<n; i++) {
            //get random x,y coords in the unit circle
            r1 = 2*rand_() - 1; r2 = 2*rand_() - 1;
            while (sqrt(r1*r1 + r2*r2) > 1.0) {
                r1 = 2*rand_() - 1; r2 = 2*rand_() - 1;
            }
            sol[i] = r1;
            sol[i+n] = r2;
            //get a random phase in [0,2*pi)
            sol[i+2*n] = 2*3.1415926535897932384626433832795*rand_();
        }
    }

    ~Swarmalator () {};

    void ode_funk (double t_, double *solin, double *fout) {
        (void)t_; //suppress warning
        //zero out the time derivatives before incrementing
        for (int i=0; i<3*n; i++) fout[i] = 0.0;
        //variable for distances
        double d, dx, dy, dtheta;
        //inefficient calculation of time derivatives
        //can be sped up by avoiding dulplicate calculations
        for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
                //distances
                dx = solin[j] - solin[i];
                dy = solin[j+n] - solin[i+n];
                d = sqrt(dx*dx + dy*dy);
                dtheta = solin[j+2*n] - solin[i+2*n];
                //time derivaties
                fout[i] += ((dx/d)*(A + J*cos(dtheta)) - B*(dx/(d*d)))/N; // dx_i/dt
                fout[j] += ((-dx/d)*(A + J*cos(-dtheta)) - B*(-dx/(d*d)))/N; // dx_j/dt
                fout[i+n] += ((dy/d)*(A + J*cos(dtheta)) - B*(dy/(d*d)))/N; // dy_i/dt
                fout[j+n] += ((-dy/d)*(A + J*cos(-dtheta)) - B*(-dy/(d*d)))/N; // dy_j/dt
                fout[i+2*n] += (K*sin(dtheta)/d)/N; //d theta_i/dt
                fout[j+2*n] -= fout[i+2*n]; //d theta_j/dt
            }
        }
    }
};

//------------------------------------------------------------------------------

#endif
