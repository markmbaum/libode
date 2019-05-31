/*
A 6th order method from Butcher with coefficients with reasonable numbers of digits

From table 6.1 of:
    E. Hairer, S. P. Nørsett, and G. Wanner. 1993. Solving Ordinary Differential Equations I (2nd Revised. Ed.): Nonstiff Problems. Springer-Verlag, Berlin, Heidelberg.
*/

#ifndef ODE_BUTCHER6_HPP
#define ODE_BUTCHER6_HPP

#include "OdeBaseRK.hpp"

class OdeButcher6 : public OdeBaseRK {
public:
    //constructor
    OdeButcher6 (unsigned long long neq_);
    //destructor
    ~OdeButcher6 ();
    //take a time step
    void step (double dt);
private:
    //arrays for time stepping routine
    double *k1, *k2, *k3, *k4, *k5, *k6, *k7;
    //coefficents of tableau
    double c2, a21,
           c3, a31, a32,
           c4, a41, a42, a43,
           c5, a51, a52, a53, a54,
           c6, a61, a62, a63, a64, a65,
           c7, a71, a72, a73, a74, a75, a76,
                b1,  b2,  b3,  b4,  b5,  b6, b7;

};


#endif
