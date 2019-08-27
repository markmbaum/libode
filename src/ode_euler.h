#ifndef ODE_EULER_H_
#define ODE_EULER_H_

#include "ode_base.h"
#include "ode_rk.h"

//!The simplest runge kutta method, forward Euler's
class OdeEuler : public OdeBase, public OdeRK {
    public:
        //!constructs
        /*!
        \param[in] neq size of ODE system
        */
        OdeEuler (unsigned long neq);
    private:
        //take a time step
        void step_ (double dt);
};

#endif
