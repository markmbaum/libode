/*
This class implements a 3rd and 4th order method with the FSAL (first same as last) property
*/

#ifndef ODE_RK43_HPP
#define ODE_RK43_HPP

#include "OdeBaseARK.hpp"

class OdeRK43 : public OdeBaseARK {
public:
    //constructor
    OdeRK43 (unsigned long long neq_);
    //destructor
    ~OdeRK43 ();
    //function for taking a single time step
    void step (double dt);
private:
    //arrays for time stepping routine
    double *k1, *k2, *k3, *k4, *k5;
    //coefficents of tableau
    double c2, a21,
           c3, a31, a32,
           c4, a41, a42, a43,
           c5, a51, a52, a53, a54,
                b1,  b2,  b3,  b4,
                d1,  d2,  d3,  d4, d5;

};

#endif
