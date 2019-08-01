//the simplest runge kutta method, forward Euler's method

#ifndef ODE_EULER_H_
#define ODE_EULER_H_

#include "ode_base.h"
#include "ode_rk.h"

class OdeEuler : public OdeBase, public OdeRK {
    public:
        //constructor
        OdeEuler (unsigned long neq);
    private:
        //take a time step
        void step_ (double dt);
};

#endif
