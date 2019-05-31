/*
THE classic Runge-Kutta 4th order method
*/

#ifndef ODE_RK4_HPP
#define ODE_RK4_HPP

#include "OdeBaseRK.hpp"

class OdeRK4 : public OdeBaseRK {
public:
    //constructor
    OdeRK4 (unsigned long long neq_);
    //destructor
    ~OdeRK4 ();
    //take a time step
    void step (double dt);
private:
    //arrays for time stepping routine
    double *k1, *k2, *k3, *k4;
    //coefficents of tableau
    double c2, a21,
           c3,      a32,
           c4,           a43,
                b1,  b2,  b3,  b4;

};


#endif
