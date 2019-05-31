//the oldest runge kutta method: Euler's explicit method

#ifndef ODE_EULER_HPP
#define ODE_EULER_HPP

#include "OdeBaseRK.hpp"

class OdeEuler : public OdeBaseRK {
public:
    //constructor
    OdeEuler (unsigned long long neq_);
    //destructor
    ~OdeEuler ();
    //temporary y values for stepping
    double *ytemp;
    //take a time step
    void step (double dt);
private:
    //vectors for time stepping
    double *k1;
};

#endif
