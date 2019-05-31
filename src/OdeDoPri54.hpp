/*
This class implements a 4th and 5th order method with the FSAL (first same as last) property
*/

#ifndef ODE_DOPRI54_HPP
#define ODE_DOPRI54_HPP

#include "OdeBaseARK.hpp"

class OdeDoPri54 : public OdeBaseARK {
public:
    //constructor
    OdeDoPri54 (unsigned long long neq_);
    //destructor
    ~OdeDoPri54 ();
    //function for taking a single time step
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
                b1,  b2,  b3,  b4,  b5,  b6,
                d1,  d2,  d3,  d4,  d5,  d6, d7;

};

#endif
