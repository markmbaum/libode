#ifndef ODE_EULER_H_
#define ODE_EULER_H_

//! \file ode_euler.h

#include "ode_adaptive.h"
#include "ode_rk.h"

//!The simplest runge kutta method, forward Euler's
/*!
    + https://en.wikipedia.org/wiki/Euler_method
*/
class OdeEuler : public OdeAdaptive, public OdeRK {
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
