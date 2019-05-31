/*
This class implements a 2nd and 3rd solver developed by Fehlberg.
*/

#ifndef ODE_RKF32_HPP
#define ODE_RKF32_HPP

#include "OdeBaseARK.hpp"

class OdeRKF32 : public OdeBaseARK {
public:
    //constructor
    OdeRKF32 (unsigned long long neq_);
    //destructor
    ~OdeRKF32 ();
    //function for taking a single time step
    void step (double dt);
private:
    //vectors for stepping
    double *k1, *k2, *k3;
    //coefficents of tableau
    double c2, a21,
           c3, a31, a32,
                b1,  b2, b3,
                d1,  d2, d3;

};

#endif
