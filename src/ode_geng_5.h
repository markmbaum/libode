/*
The fifth-order, symplectic, fully-implicit Geng integrator with 3 stages
*/

#ifndef ODE_GENG_5_H_
#define ODE_GENG_5_H_

#include "ode_base.h"
#include "ode_irk.h"
#include "ode_newton_bridge.h"

class OdeGeng5;

class Geng5Newton : public OdeNewtonIRK<OdeGeng5> {
    public:
        Geng5Newton (unsigned long neq, unsigned long nnew, OdeGeng5 *integrator) : OdeNewtonIRK (neq, nnew, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeGeng5 : public OdeBase, private OdeIRK {
    //friends!
    friend class OdeNewtonBridge<OdeGeng5>;
    friend class OdeNewtonIRK<OdeGeng5>;

    public:
        OdeGeng5 (unsigned long neq);
        ~OdeGeng5 ();
        Geng5Newton get_newton () { return(*newton_); }
    private:
        double **a;
        double *b;
        Geng5Newton *newton_;
        void step_ (double dt);
};

#endif
