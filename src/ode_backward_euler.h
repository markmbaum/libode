/*
Backward Euler's method, unconditionally stable but relatively inaccurate
*/

#ifndef ODE_BACKWARD_EULER_H_
#define ODE_BACKWARD_EULER_H_

#include "ode_base.h"
#include "ode_irk.h"
#include "ode_newton_bridge.h"

class OdeBackwardEuler;

class BackwardEulerNewton : public OdeNewtonIRK<OdeBackwardEuler> {
    public:
        BackwardEulerNewton (unsigned long neq, unsigned long nnew, OdeBackwardEuler *integrator) : OdeNewtonIRK (neq, nnew, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeBackwardEuler : public OdeBase, private OdeIRK {
    //friends!
    friend class OdeNewtonBridge<OdeBackwardEuler>;
    friend class OdeNewtonIRK<OdeBackwardEuler>;

    public:
        OdeBackwardEuler (unsigned long neq);
        ~OdeBackwardEuler ();
        BackwardEulerNewton get_newton () { return(*newton_); }
    private:
        double **a;
        double *b;
        BackwardEulerNewton *newton_;
        void step_ (double dt);
};

#endif
