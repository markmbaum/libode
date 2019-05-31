/*
A second order adaptive solver from Fehlberg with the FSAL property
*/

#ifndef ODE_RKF21_HPP
#define ODE_RKF21_HPP

#include "OdeBaseARK.hpp"

class OdeRKF21 : public OdeBaseARK {
public:
    //constructor
    OdeRKF21 (unsigned long long neq_);
    //destructor
    ~OdeRKF21 ();
    //function for taking a single time step
    void step (double dt);
private:
    //vectors for stepping
    double *k1, *k2, *k3;
    //coefficents of tableau
    double c2, a21,
           c3, a31, a32,
                b1,  b2,
                d1,  d2, d3;
};

#endif
